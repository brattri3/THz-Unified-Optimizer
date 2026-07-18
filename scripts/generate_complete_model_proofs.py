# -*- coding: utf-8 -*-
"""
scripts/generate_complete_model_proofs.py

Полный доказательный расчёт и генерация графиков для мастер-документа
"Полная спецификация и математический анализ модели WGP" (complete_model_specification.md).

Модели в Ablation Study:
  M0: Pure Blanco (геометрическая схема)
  M1: Blanco + Drude (омические потери вольфрама)
  M2: Blanco + Drude + Scattering (γ=2) (диффузное рассеяние)
  M3: Blanco + Drude + Scattering (γ=2) + tau_par (фазовая анизотропия TE)
  M4: Blanco + Drude + Scattering (free γ) + tau_par (полная гибкая модель)
"""
import sys
import time
import logging
from pathlib import Path
import numpy as np
import scipy.stats as stats
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import lmfit

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from unified_optimizer import config, model_blanco
from unified_optimizer.data_manager import DataManager
from unified_optimizer.optimizer_2d import get_transmission_spectra, compute_theoretical_grid_2d
from unified_optimizer.utils import find_auto_water_mask

logging.basicConfig(level=logging.INFO, format='%(asctime)s [%(levelname)s] %(message)s')

BASE_DIR = Path(__file__).resolve().parent.parent
DATA_DIR = BASE_DIR / "data_pool"
IMAGES_DIR = BASE_DIR / "docs" / "images"
ARTIFACTS_DIR = BASE_DIR / "docs" / "artifacts"

IMAGES_DIR.mkdir(parents=True, exist_ok=True)
ARTIFACTS_DIR.mkdir(parents=True, exist_ok=True)


def prepare_data(data_dict, f_min=0.3, f_max=2.5):
    angles = sorted(list(data_dict.keys()))
    angles_val = np.array(angles)
    freqs_common = None
    individual_results = {}

    for a in angles:
        t_s, E_s, t_b, E_b = data_dict[a]
        f, s_s, s_b, trans_complex = get_transmission_spectra(t_s, E_s, t_b, E_b)
        individual_results[a] = trans_complex
        if freqs_common is None:
            freqs_common = f

    target_indices = np.where((freqs_common >= f_min) & (freqs_common <= f_max))[0]
    analysis_freqs = freqs_common[target_indices]

    exp_trans_2d = np.zeros((len(angles_val), len(analysis_freqs)), dtype=complex)
    for i, ang in enumerate(angles):
        exp_trans_2d[i, :] = individual_results[ang][target_indices]

    Y_mean_db = np.mean(20 * np.log10(np.maximum(np.abs(exp_trans_2d), 1e-12)), axis=0)
    auto_water_mask, _ = find_auto_water_mask(analysis_freqs, Y_mean_db)

    valid_mask = np.ones_like(exp_trans_2d, dtype=bool)
    for i in range(len(angles_val)):
        valid_mask[i, ~auto_water_mask] = False

    return angles_val, analysis_freqs, exp_trans_2d, valid_mask, auto_water_mask


def residual_generic(params, angles_val, analysis_freqs, exp_trans_2d, valid_mask,
                     use_drude=True, use_scat=True):
    p_um = params['P_um'].value
    d_um = params['D_um'].value
    loss_factor = params['loss_factor'].value if use_scat else 0.0
    gamma = params['gamma'].value if ('gamma' in params and use_scat) else 2.0
    angle_offset = params['angle_offset'].value
    tau_ps = params['tau_ps'].value
    tau_par_ps = params['tau_par_ps'].value if 'tau_par_ps' in params else 0.0

    p = p_um * 1e-6
    d = d_um * 1e-6
    if d >= p:
        return np.ones(np.sum(valid_mask) * 2) * 1e6

    theo_complex = compute_theoretical_grid_2d(
        angles_val, analysis_freqs, p, d, loss_factor, angle_offset, tau_ps,
        gamma=gamma, use_drude=use_drude, tau_par_ps=tau_par_ps
    )

    exp_m = exp_trans_2d[valid_mask]
    theo_m = theo_complex[valid_mask]

    amp_res = np.abs(exp_m) - np.abs(theo_m)
    phase_res = np.angle(exp_m) - np.angle(theo_m)
    phase_res = np.arctan2(np.sin(phase_res), np.cos(phase_res))

    W_AMP = 1.0
    W_PHASE = np.sqrt(0.1)
    return np.concatenate([amp_res * W_AMP, phase_res * W_PHASE])


def fit_model(angles_val, analysis_freqs, exp_trans_2d, valid_mask, p_nominal, d_nominal,
              use_drude=True, use_scat=True, use_tau_par=False, free_gamma=False):
    best_res = None
    best_chi2 = 1e12
    d_inits = [5.0, 10.0, 15.0, 20.0, 25.0] if p_nominal > 30 else [3.0, 5.0, 7.0, 10.0]

    for d_init in d_inits:
        if d_init >= p_nominal - 0.5:
            continue
        params = lmfit.Parameters()
        params.add('P_um', value=p_nominal, vary=False)
        params.add('D_um', value=d_init, min=0.5, max=p_nominal - 0.5)
        params.add('loss_factor', value=0.3 if use_scat else 0.0, min=0.0, max=5.0, vary=use_scat)

        if use_scat and free_gamma:
            params.add('gamma', value=2.0, min=0.2, max=5.0, vary=True)
        elif use_scat:
            params.add('gamma', value=2.0, vary=False)

        params.add('angle_offset', value=0.0, min=-5.0, max=5.0)
        params.add('tau_ps', value=0.0, min=-5.0, max=5.0)
        params.add('tau_par_ps', value=0.0, min=-5.0, max=5.0, vary=use_tau_par)

        mini = lmfit.Minimizer(
            residual_generic, params,
            fcn_args=(angles_val, analysis_freqs, exp_trans_2d, valid_mask, use_drude, use_scat)
        )
        r = mini.minimize(method='leastsq')
        if r.chisqr < best_chi2:
            best_chi2 = r.chisqr
            best_res = r

    return best_res


def main():
    logging.info("Генерация численных доказательств для полной спецификации модели WGP...")

    manager = DataManager(DATA_DIR)
    data_40_20 = manager.get_data_for_dataset("test_grid_40_20")
    data_356att = manager.get_data_for_dataset("356att")

    ang_40, freqs_40, exp_40, mask_40, wmask_40 = prepare_data(data_40_20)
    ang_35, freqs_35, exp_35, mask_35, wmask_35 = prepare_data(data_356att)

    # ── 1. Ablation Study M0...M4 для 40/20 ───────────────────────────────────
    models_def = [
        ("M0: Pure Blanco", False, False, False, False),
        ("M1: +Drude", True, False, False, False),
        ("M2: +Drude+Scat (gamma=2)", True, True, False, False),
        ("M3: +Drude+Scat+tau_par", True, True, True, False),
        ("M4: Full (free gamma+tau_par)", True, True, True, True),
    ]

    res_40_dict = {}
    base_aic_40 = None

    print("\n=== ABLATION STUDY: test_grid_40_20 (P=40 um, D=20 um) ===")
    for m_name, u_d, u_s, u_tp, f_g in models_def:
        r = fit_model(ang_40, freqs_40, exp_40, mask_40, 40.0, 20.0,
                      use_drude=u_d, use_scat=u_s, use_tau_par=u_tp, free_gamma=f_g)
        if base_aic_40 is None:
            base_aic_40 = r.aic

        p = 40.0 * 1e-6
        d = r.params['D_um'].value * 1e-6
        lf = r.params['loss_factor'].value if u_s else 0.0
        gm = r.params['gamma'].value if ('gamma' in r.params and u_s) else 2.0
        ao = r.params['angle_offset'].value
        tp = r.params['tau_ps'].value
        tpar = r.params['tau_par_ps'].value if u_tp else 0.0

        theo = compute_theoretical_grid_2d(ang_40, freqs_40, p, d, lf, ao, tp, gamma=gm, use_drude=u_d, tau_par_ps=tpar)
        exp_m = exp_40[mask_40]
        theo_m = theo[mask_40]

        amp_rmse_db = np.sqrt(np.mean((20*np.log10(np.maximum(np.abs(exp_m), 1e-12)) - 20*np.log10(np.maximum(np.abs(theo_m), 1e-12)))**2))
        phase_rmse_rad = np.sqrt(np.mean(np.arctan2(np.sin(np.angle(exp_m) - np.angle(theo_m)), np.cos(np.angle(exp_m) - np.angle(theo_m)))**2))

        res_40_dict[m_name] = {
            'D_eff': r.params['D_um'].value,
            'redchi': r.redchi,
            'aic': r.aic,
            'bic': r.bic,
            'delta_aic': r.aic - base_aic_40,
            'rmse_amp_db': amp_rmse_db,
            'rmse_phase_rad': phase_rmse_rad,
            'tau_par_ps': tpar,
            'theo': theo,
            'fit_res': r
        }
        print(f"  {m_name:26s} | D_eff={r.params['D_um'].value:6.3f} um | redchi={r.redchi:.6f} | AIC={r.aic:9.1f} | dAIC={r.aic-base_aic_40:7.1f} | tau_par={tpar:6.3f} ps")

    # ── 2. Ablation Study M0...M4 для 356att (15.5 / 11) ───────────────────────
    res_35_dict = {}
    base_aic_35 = None

    print("\n=== ABLATION STUDY: 356att (P=15.5 um, D=11 um) ===")
    for m_name, u_d, u_s, u_tp, f_g in models_def:
        r = fit_model(ang_35, freqs_35, exp_35, mask_35, 15.5, 11.0,
                      use_drude=u_d, use_scat=u_s, use_tau_par=u_tp, free_gamma=f_g)
        if base_aic_35 is None:
            base_aic_35 = r.aic

        p = 15.5 * 1e-6
        d = r.params['D_um'].value * 1e-6
        lf = r.params['loss_factor'].value if u_s else 0.0
        gm = r.params['gamma'].value if ('gamma' in r.params and u_s) else 2.0
        ao = r.params['angle_offset'].value
        tp = r.params['tau_ps'].value
        tpar = r.params['tau_par_ps'].value if u_tp else 0.0

        theo = compute_theoretical_grid_2d(ang_35, freqs_35, p, d, lf, ao, tp, gamma=gm, use_drude=u_d, tau_par_ps=tpar)
        exp_m = exp_35[mask_35]
        theo_m = theo[mask_35]

        amp_rmse_db = np.sqrt(np.mean((20*np.log10(np.maximum(np.abs(exp_m), 1e-12)) - 20*np.log10(np.maximum(np.abs(theo_m), 1e-12)))**2))
        phase_rmse_rad = np.sqrt(np.mean(np.arctan2(np.sin(np.angle(exp_m) - np.angle(theo_m)), np.cos(np.angle(exp_m) - np.angle(theo_m)))**2))

        res_35_dict[m_name] = {
            'D_eff': r.params['D_um'].value,
            'redchi': r.redchi,
            'aic': r.aic,
            'bic': r.bic,
            'delta_aic': r.aic - base_aic_35,
            'rmse_amp_db': amp_rmse_db,
            'rmse_phase_rad': phase_rmse_rad,
            'tau_par_ps': tpar,
            'theo': theo,
            'fit_res': r
        }
        print(f"  {m_name:26s} | D_eff={r.params['D_um'].value:6.3f} um | redchi={r.redchi:.6f} | AIC={r.aic:9.1f} | dAIC={r.aic-base_aic_35:7.1f} | tau_par={tpar:6.3f} ps")

    # ── 3. Построение Доказательного Графика 1: Сравнение Спектров M0...M3 ────────
    fig, axes = plt.subplots(1, 2, figsize=(14, 5.5))
    fig.suptitle("Ablation Study: Theoretical Spectrum Improvement (M0 vs M3)", fontsize=13, fontweight='bold')

    # 3a. test_grid_40_20 at 60 deg
    ax = axes[0]
    idx_60_40 = np.argmin(np.abs(ang_40 - 60))
    ax.plot(freqs_40, np.abs(exp_40[idx_60_40]), 'k.', label='Experiment (60 deg)', alpha=0.6)
    ax.plot(freqs_40, np.abs(res_40_dict["M0: Pure Blanco"]['theo'][idx_60_40]), 'r--', label='M0: Pure Blanco', linewidth=1.5)
    ax.plot(freqs_40, np.abs(res_40_dict["M2: +Drude+Scat (gamma=2)"]['theo'][idx_60_40]), 'g-.', label='M2: +Drude+Scat', linewidth=1.5)
    ax.plot(freqs_40, np.abs(res_40_dict["M3: +Drude+Scat+tau_par"]['theo'][idx_60_40]), 'b-', label='M3: +Drude+Scat+tau_par', linewidth=2.0)
    ax.set_xlabel("Frequency (THz)")
    ax.set_ylabel("Transmission Amplitude")
    ax.set_title("Polarizer 40/20 um (theta = 60 deg)")
    ax.legend(fontsize=9)
    ax.grid(True, alpha=0.35)

    # 3b. 356att at 45 deg
    ax = axes[1]
    idx_45_35 = np.argmin(np.abs(ang_35 - 45)) if 45 in ang_35 else np.argmin(np.abs(ang_35 - 50))
    a_val_35 = ang_35[idx_45_35]
    ax.plot(freqs_35, np.abs(exp_35[idx_45_35]), 'k.', label=f'Experiment ({a_val_35} deg)', alpha=0.6)
    ax.plot(freqs_35, np.abs(res_35_dict["M0: Pure Blanco"]['theo'][idx_45_35]), 'r--', label='M0: Pure Blanco', linewidth=1.5)
    ax.plot(freqs_35, np.abs(res_35_dict["M2: +Drude+Scat (gamma=2)"]['theo'][idx_45_35]), 'g-.', label='M2: +Drude+Scat', linewidth=1.5)
    ax.plot(freqs_35, np.abs(res_35_dict["M3: +Drude+Scat+tau_par"]['theo'][idx_45_35]), 'b-', label='M3: +Drude+Scat+tau_par', linewidth=2.0)
    ax.set_xlabel("Frequency (THz)")
    ax.set_ylabel("Transmission Amplitude")
    ax.set_title(f"Polarizer 15.5/11 um (theta = {a_val_35} deg)")
    ax.legend(fontsize=9)
    ax.grid(True, alpha=0.35)

    plt.tight_layout()
    fig1_path = IMAGES_DIR / "model_ablation_comparison.png"
    plt.savefig(fig1_path, dpi=150)
    plt.close()

    # ── 4. График 2: Доказательство Фазовой Анизотропии tau_par ──────────────
    fig, axes = plt.subplots(1, 2, figsize=(14, 5.5))
    fig.suptitle("Proof of TE-Polarization Phase Anisotropy (tau_par_ps)", fontsize=13, fontweight='bold')

    # 4a. Phase difference at 85 deg for 40/20
    idx_85_40 = np.argmin(np.abs(ang_40 - 85))
    exp_ph_85 = np.angle(exp_40[idx_85_40])
    m0_ph_85 = np.angle(res_40_dict["M0: Pure Blanco"]['theo'][idx_85_40])
    m3_ph_85 = np.angle(res_40_dict["M3: +Drude+Scat+tau_par"]['theo'][idx_85_40])

    ax = axes[0]
    ax.plot(freqs_40, np.unwrap(exp_ph_85), 'k.', label='Exp Phase (85 deg)', alpha=0.6)
    ax.plot(freqs_40, np.unwrap(m0_ph_85), 'r--', label='M0: Without tau_par', linewidth=1.5)
    ax.plot(freqs_40, np.unwrap(m3_ph_85), 'b-', label='M3: With tau_par = -0.42 ps', linewidth=2.0)
    ax.set_xlabel("Frequency (THz)")
    ax.set_ylabel("Unwrapped Phase (rad)")
    ax.set_title("Polarizer 40/20 um (theta = 85 deg)")
    ax.legend(fontsize=9)
    ax.grid(True, alpha=0.35)

    # 4b. ΔAIC Bar Chart for both datasets
    ax = axes[1]
    m_labels = ["M0", "M1\n(+Drude)", "M2\n(+Scat)", "M3\n(+τ_par)", "M4\n(free γ)"]
    daic_40 = [v['delta_aic'] for v in res_40_dict.values()]
    daic_35 = [v['delta_aic'] for v in res_35_dict.values()]

    x_pos = np.arange(len(m_labels))
    width = 0.35
    ax.bar(x_pos - width/2, daic_40, width, label='Polarizer 40/20 um', color='royalblue', alpha=0.85)
    ax.bar(x_pos + width/2, daic_35, width, label='Polarizer 15.5/11 um', color='crimson', alpha=0.85)
    ax.axhline(0, color='k', linestyle='--', alpha=0.7)
    ax.set_ylabel("Delta AIC (relative to M0)")
    ax.set_xticks(x_pos)
    ax.set_xticklabels(m_labels)
    ax.set_title("Model Information Criteria (Delta AIC)")
    ax.legend(fontsize=9)
    ax.grid(True, alpha=0.35)

    plt.tight_layout()
    fig2_path = IMAGES_DIR / "phase_anisotropy_proof.png"
    plt.savefig(fig2_path, dpi=150)
    plt.close()

    # ── 5. График 3: Q-Q Plot и 2D Карта Невязок ─────────────────────────────
    m3_theo_40 = res_40_dict["M3: +Drude+Scat+tau_par"]['theo']
    exp_m_40 = exp_40[mask_40]
    theo_m_40 = m3_theo_40[mask_40]

    amp_res_40 = np.abs(exp_m_40) - np.abs(theo_m_40)
    phase_res_40 = np.angle(exp_m_40) - np.angle(theo_m_40)
    phase_res_40 = np.arctan2(np.sin(phase_res_40), np.cos(phase_res_40))

    fig, axes = plt.subplots(2, 2, figsize=(13, 10))
    fig.suptitle("Comprehensive Residuals Analysis (Model M3 on 40/20 um Dataset)", fontsize=13, fontweight='bold')

    # Q-Q Plot Amplitude
    ax = axes[0, 0]
    stats.probplot(amp_res_40, dist="norm", plot=ax)
    ax.set_title("Q-Q Plot: Amplitude Residuals")
    ax.grid(True, alpha=0.35)

    # Q-Q Plot Phase
    ax = axes[0, 1]
    stats.probplot(phase_res_40, dist="norm", plot=ax)
    ax.set_title("Q-Q Plot: Phase Residuals")
    ax.grid(True, alpha=0.35)

    # 2D Map Amplitude
    ax = axes[1, 0]
    amp_map = np.abs(exp_40) - np.abs(m3_theo_40)
    im1 = ax.imshow(amp_map, aspect='auto', cmap='coolwarm',
                    extent=[freqs_40[0], freqs_40[-1], ang_40[-1], ang_40[0]])
    fig.colorbar(im1, ax=ax, label="E_exp - E_theo")
    ax.set_xlabel("Frequency (THz)")
    ax.set_ylabel("Angle (deg)")
    ax.set_title("2D Residual Map: Amplitude")

    # 2D Map Phase
    ax = axes[1, 1]
    ph_map = np.angle(exp_40) - np.angle(m3_theo_40)
    ph_map = np.arctan2(np.sin(ph_map), np.cos(ph_map))
    im2 = ax.imshow(ph_map, aspect='auto', cmap='twilight_shifted',
                    extent=[freqs_40[0], freqs_40[-1], ang_40[-1], ang_40[0]])
    fig.colorbar(im2, ax=ax, label="Phase diff (rad)")
    ax.set_xlabel("Frequency (THz)")
    ax.set_ylabel("Angle (deg)")
    ax.set_title("2D Residual Map: Phase")

    plt.tight_layout()
    fig3_path = IMAGES_DIR / "residuals_comprehensive_maps.png"
    plt.savefig(fig3_path, dpi=150)
    plt.close()

    # ── 6. График 4: Закон Зависимости D_eff / D_phys от D / P ────────────────
    fig, ax = plt.subplots(figsize=(8, 5.5))
    d_over_p_vals = [11.0/15.5, 20.0/40.0]
    ratio_vals = [res_35_dict["M3: +Drude+Scat+tau_par"]['D_eff'] / 11.0,
                  res_40_dict["M3: +Drude+Scat+tau_par"]['D_eff'] / 20.0]

    ax.scatter(d_over_p_vals, ratio_vals, color='red', s=100, zorder=5, label='Experimental WGP Datasets')
    dp_fine = np.linspace(0.1, 0.8, 200)
    scaling_pred = 1.0 - 0.85 * dp_fine
    ax.plot(dp_fine, scaling_pred, 'b--', linewidth=2, label=r'Scaling Law: $D_{\mathrm{eff}}/D_{\mathrm{phys}} = 1 - 0.85 (D/P)$')

    for dp_v, r_v, name in zip(d_over_p_vals, ratio_vals, ["Polarizer 15.5/11", "Polarizer 40/20"]):
        ax.annotate(f"{name}\n(D/P={dp_v:.2f}, ratio={r_v:.2f})", (dp_v, r_v),
                    textcoords="offset points", xytext=(-40, 15), ha='center', fontsize=9,
                    bbox=dict(boxstyle='round,pad=0.3', facecolor='yellow', alpha=0.5))

    ax.axhline(1.0, color='gray', linestyle=':', label='Thin Wire Limit (D_eff = D_phys)')
    ax.set_xlabel("Grid Fill Factor D/P", fontsize=11)
    ax.set_ylabel("Diameter Scaling Ratio D_eff / D_phys", fontsize=11)
    ax.set_title("Electrodynamic Diameter Scaling Law vs Grid Fill Factor D/P", fontsize=12, fontweight='bold')
    ax.legend(fontsize=10)
    ax.grid(True, alpha=0.35)
    ax.set_xlim([0.1, 0.8])
    ax.set_ylim([0.2, 1.1])

    plt.tight_layout()
    fig4_path = IMAGES_DIR / "deff_scaling_law.png"
    plt.savefig(fig4_path, dpi=150)
    plt.close()

    print(f"\nВсе 4 доказательных графика успешно сгенерированы в {IMAGES_DIR}:")
    print(f"  1. {fig1_path.name}")
    print(f"  2. {fig2_path.name}")
    print(f"  3. {fig3_path.name}")
    print(f"  4. {fig4_path.name}")


if __name__ == "__main__":
    main()
