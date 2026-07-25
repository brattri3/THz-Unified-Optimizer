"""
Reusable fitting/metrics library for the autonomous research (T1 baseline + ablations).
Uses ONLY public functions of the original package (does not edit it).

Model toggles let us build the ablation ladder M0..M4 and future variants.
"""
from __future__ import annotations
import sys, os, json
from pathlib import Path
import numpy as np
import lmfit
from scipy import stats

REPO = Path(__file__).resolve().parents[2]   # THz-Unified-Optimizer
sys.path.insert(0, str(REPO))

from unified_optimizer import config
from unified_optimizer.data_manager import DataManager
from unified_optimizer import model_blanco
from unified_optimizer.optimizer_2d import get_transmission_spectra, compute_theoretical_grid_2d
from unified_optimizer.utils import find_auto_water_mask

# nominal geometry per dataset (P_um, D_phys_um)
GEOMETRY = {
    "356att":          {"P_um": 15.5, "D_phys_um": 11.0},
    "test_grid_40_20": {"P_um": 38.8, "D_phys_um": 18.0},  # micrograph GT (was nominal 40/20)
    # T15 external validation — micrograph ground-truth (2026-07-24, REVISED diam).
    # Micrograph fixes PERIOD tightly; DIAMETER is a lower-ish prior (mirror-blik
    # underestimates; corrected by measuring full trace width @10%). PureWave:
    # P=25.5+-0.9, D=10.6+-0.5 (D/P~0.42). THz measurements NOT yet in data_pool.
    "purewave":        {"P_um": 25.5, "D_phys_um": 10.6},
    # test_grid_40_20 micrograph GT (2026-07-24): P=38.8+-0.3, D=18+-2 (meas 17
    #   @10% + specular-underread correction toward nominal 20), D/P~0.46. P is a
    #   HARD pin; treat D as a wide prior. NOTE: wires span ~200 um along the beam
    #   (non-coplanar) -> candidate angular term (HN12) / delayed 90deg leakage.
    #   Legacy T1-T14 used nominal 40/20 -> re-run those to compare if needed.
    # Specac: micrograph (PureWave-cal transfer, 2026-07-25) P=24.9 (midmag FFT+lattice,
    #   reliable), D glint-underread ~8-12 -> ~14+-3 corrected (wide prior), D/P~0.55.
    #   THz DATA READY in data_pool/specac/ (13 angles 0-100deg incl 84/90/96 near shadow;
    #   drift-aware bg1..bg7). Model-independent leakage floor eta=0.0357, SNR 68x (best).
    "specac":          {"P_um": 24.9, "D_phys_um": 14.0},
    # series1..series5: THz data present but geometry UNKNOWN (blocks full M5 fit).
}

W_AMP = 1.0
W_PHASE = np.sqrt(0.1)


def build_experiment(data_dict):
    """Replicate the 2D experimental grid + water mask (from optimizer_lmfit)."""
    angles = sorted(data_dict.keys())
    if config.ANGLES_LIMIT_2D is not None:
        lo, hi = config.ANGLES_LIMIT_2D
        angles = [a for a in angles if lo <= a <= hi]
    freqs_common, individual = None, {}
    for a in angles:
        t_s, E_s, t_b, E_b = data_dict[a]
        f, s_s, s_b, trans = get_transmission_spectra(t_s, E_s, t_b, E_b)
        individual[a] = trans
        if freqs_common is None:
            freqs_common = f
    f_max = getattr(config, "F_MAX_2D", config.F_MAX)
    idx = np.where((freqs_common >= config.F_MIN) & (freqs_common <= f_max))[0]
    analysis_freqs = freqs_common[idx]
    angles_val = np.array(angles)
    exp = np.zeros((len(angles_val), len(analysis_freqs)), dtype=complex)
    for i, a in enumerate(angles):
        exp[i, :] = individual[a][idx]
    Ymean_db = np.mean(20*np.log10(np.maximum(np.abs(exp), 1e-12)), axis=0)
    water_mask, _ = find_auto_water_mask(analysis_freqs, Ymean_db)
    valid = np.ones_like(exp, dtype=bool)
    for i in range(len(angles_val)):
        valid[i, ~water_mask] = False
    return angles_val, analysis_freqs, exp, valid


def residual(params, angles_val, freqs, exp, valid, use_drude):
    p = params["P_um"].value * 1e-6
    d = params["D_um"].value * 1e-6
    if d >= p:
        return np.ones(int(np.sum(valid))*2) * 1e6
    gamma = params["gamma"].value if "gamma" in params else 2.0
    theo = compute_theoretical_grid_2d(
        angles_val, freqs, p, d,
        params["loss_factor"].value, params["angle_offset"].value,
        params["tau_ps"].value, gamma=gamma, use_drude=use_drude,
        tau_par_ps=params["tau_par_ps"].value if "tau_par_ps" in params else 0.0)
    em, tm = exp[valid], theo[valid]
    amp = np.abs(em) - np.abs(tm)
    ph = np.angle(em) - np.angle(tm)
    ph = np.arctan2(np.sin(ph), np.cos(ph))
    return np.concatenate([amp*W_AMP, ph*W_PHASE])


def fit_model(dataset, use_drude=True, use_scattering=True, free_gamma=True,
              use_tau_par=True, label="M?"):
    """Fit one model variant to one dataset; return dict of results/metrics/residuals."""
    geo = GEOMETRY[dataset]
    dm = DataManager(config.DATA_DIR)
    data = dm.get_data_for_dataset(dataset)
    angles_val, freqs, exp, valid = build_experiment(data)

    p_um = geo["P_um"]
    params = lmfit.Parameters()
    params.add("P_um", value=p_um, min=1.0, max=100.0, vary=False)  # fixed by optical control
    d_init = model_blanco.estimate_deff_initial(p_um, geo["D_phys_um"])
    params.add("D_um", value=d_init, min=0.1, max=p_um-0.5)
    params.add("loss_factor", value=0.3 if use_scattering else 0.0, min=0.0, max=5.0, vary=use_scattering)
    params.add("gamma", value=2.0, min=-5.0, max=10.0, vary=use_scattering and free_gamma)
    params.add("angle_offset", value=0.0, min=-10.0, max=10.0)
    params.add("tau_ps", value=0.0, min=-10.0, max=10.0)
    params.add("tau_par_ps", value=0.0, min=-5.0, max=5.0, vary=use_tau_par)

    mini = lmfit.Minimizer(residual, params,
                           fcn_args=(angles_val, freqs, exp, valid, use_drude))
    res = mini.minimize(method="leastsq")

    # residual decomposition (unweighted, physical units)
    p = res.params["P_um"].value*1e-6; d = res.params["D_um"].value*1e-6
    gamma = res.params["gamma"].value
    theo = compute_theoretical_grid_2d(angles_val, freqs, p, d,
              res.params["loss_factor"].value, res.params["angle_offset"].value,
              res.params["tau_ps"].value, gamma=gamma, use_drude=use_drude,
              tau_par_ps=res.params["tau_par_ps"].value if use_tau_par else 0.0)
    em, tm = exp[valid], theo[valid]
    amp_res = np.abs(em) - np.abs(tm)
    amp_res_db = 20*np.log10(np.maximum(np.abs(em),1e-12)) - 20*np.log10(np.maximum(np.abs(tm),1e-12))
    ph_res = np.arctan2(np.sin(np.angle(em)-np.angle(tm)), np.cos(np.angle(em)-np.angle(tm)))

    def safe_norm(x):
        x = x[np.isfinite(x)]
        if len(x) < 8: return {"n": int(len(x))}
        sw = stats.shapiro(x[:5000]); jb = stats.jarque_bera(x)
        dw = float(np.sum(np.diff(x)**2)/np.sum(x**2)) if np.sum(x**2)>0 else None
        return {"n": int(len(x)), "shapiro_W": float(sw.statistic), "shapiro_p": float(sw.pvalue),
                "jb_stat": float(jb.statistic), "jb_p": float(jb.pvalue), "durbin_watson": dw}

    metrics = {
        "label": label, "dataset": dataset, "success": bool(res.success),
        "nfev": int(res.nfev), "ndata": int(res.ndata), "nvarys": int(res.nvarys),
        "chisqr": float(res.chisqr), "redchi": float(res.redchi),
        "aic": float(res.aic), "bic": float(res.bic),
        "rmse_amp_lin": float(np.sqrt(np.mean(amp_res**2))),
        "rmse_amp_db": float(np.sqrt(np.mean(amp_res_db[np.isfinite(amp_res_db)]**2))),
        "rmse_phase_rad": float(np.sqrt(np.mean(ph_res**2))),
        "params": {k: {"value": float(v.value), "stderr": (float(v.stderr) if v.stderr else None), "vary": bool(v.vary)}
                   for k, v in res.params.items()},
        "residual_normality": {"amp": safe_norm(amp_res), "phase": safe_norm(ph_res)},
    }
    arrays = {"angles": angles_val, "freqs": freqs, "exp": exp, "valid": valid,
              "theo": theo, "amp_res_db": amp_res_db}
    return metrics, arrays


if __name__ == "__main__":
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("dataset")
    a = ap.parse_args()
    m, _ = fit_model(a.dataset, label="M3", free_gamma=False)
    print(json.dumps(m, ensure_ascii=False, indent=2))
