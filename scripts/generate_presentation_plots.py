# -*- coding: utf-8 -*-
"""
scripts/generate_presentation_plots.py
Generirovanie prezentacionnyh grafikov s russkimi podpisyami,
diapazon 0.2-1.5 TGc, dva grafika vertikal'no.
"""
import sys
import os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from unified_optimizer.data_manager import DataManager
from unified_optimizer.optimizer_2d import get_transmission_spectra, compute_theoretical_grid_2d

BASE_DIR = Path(__file__).resolve().parent.parent
DATA_DIR = BASE_DIR / "data_pool"
OUT_DIR  = BASE_DIR / "docs" / "images"
OUT_DIR.mkdir(parents=True, exist_ok=True)

plt.rcParams.update({
    'font.family': 'DejaVu Sans',
    'axes.titlesize': 13,
    'axes.labelsize': 11,
    'xtick.labelsize': 10,
    'ytick.labelsize': 10,
    'legend.fontsize': 9,
    'figure.dpi': 150,
})

F_MIN, F_MAX = 0.2, 1.5

COLORS = {
    'exp':  '#555555',
    'M0':   '#e74c3c',
    'M2':   '#2ecc71',
    'M3':   '#2980b9',
}

def load_dataset(name):
    dm = DataManager(str(DATA_DIR))
    return dm.get_data_for_dataset(name)

def get_spectra(data_dict, f_min=0.2, f_max=1.5):
    angles = sorted(list(data_dict.keys()))
    freqs_common = None
    exp_by_angle = {}
    for a in angles:
        t_s, E_s, t_b, E_b = data_dict[a]
        f, s_s, s_b, trans = get_transmission_spectra(t_s, E_s, t_b, E_b)
        if freqs_common is None:
            freqs_common = f
        exp_by_angle[a] = trans
    mask = (freqs_common >= f_min) & (freqs_common <= f_max)
    freqs = freqs_common[mask]
    exp = {a: exp_by_angle[a][mask] for a in angles}
    return angles, freqs, exp

def theo_grid(angles_calc, freqs, params, use_drude=True, use_scat=True, tau_par=0.0):
    p = params['P_um'] * 1e-6
    d = params['D_um'] * 1e-6
    lf = params.get('loss_factor', 0.0) if use_scat else 0.0
    th_off = params.get('angle_offset', 0.0)
    tau_ps = params.get('tau_ps', 0.0)
    gamma = params.get('gamma', 2.0)
    angles_arr = np.array(angles_calc)
    grid = compute_theoretical_grid_2d(
        angles_arr, freqs, p, d, lf, th_off, tau_ps,
        gamma=gamma, use_drude=use_drude, tau_par_ps=tau_par
    )
    return {a: grid[i] for i, a in enumerate(angles_calc)}


# ============================================================
# SLIDE 3: Ablation M0 vs M3 - dva grafika vertikal'no
# ============================================================
def plot_slide3_ablation():
    print("[Slide3] Ablation M0 vs M3 (vertikal'no)...")
    try:
        data = load_dataset("356att")
    except Exception as e:
        print("  OSHIBKA: " + str(e))
        return

    angles, freqs, exp = get_spectra(data, F_MIN, F_MAX)

    params_M0 = dict(P_um=15.5, D_um=3.962, loss_factor=0.0, angle_offset=0.0, tau_ps=0.0, gamma=2.0)
    params_M3 = dict(P_um=15.5, D_um=4.398, loss_factor=0.089, angle_offset=0.80, tau_ps=-0.28, gamma=2.0)
    tau_par_M3 = -0.099

    angle_0  = min(angles, key=lambda a: abs(a - 0))
    angle_60 = min(angles, key=lambda a: abs(a - 60))
    angles_sel = [angle_0, angle_60]

    theo_M0 = theo_grid(angles_sel, freqs, params_M0, use_drude=False, use_scat=False)
    theo_M3 = theo_grid(angles_sel, freqs, params_M3, tau_par=tau_par_M3)

    fig, axes = plt.subplots(2, 1, figsize=(6.5, 7.0), sharex=True)
    fig.subplots_adjust(hspace=0.38, left=0.14, right=0.96, top=0.91, bottom=0.09)
    fig.suptitle("Поляризатор 15.5/11 мкм: улучшение аппроксимации M0→M3",
                 fontsize=12, fontweight='bold')

    titles = {angle_0: "Угол поворота θ = 0° (максим. пропускание)",
              angle_60: "Угол поворота θ = 60° (промежуточное ослабление)"}

    for ax, angle in zip(axes, angles_sel):
        E_exp = np.abs(exp[angle])
        E_M0  = np.abs(theo_M0[angle])
        E_M3  = np.abs(theo_M3[angle])
        ax.scatter(freqs, E_exp, s=8, color=COLORS['exp'], label='Эксперимент', zorder=3, alpha=0.7)
        ax.plot(freqs, E_M0, '--', color=COLORS['M0'], lw=1.8, label='M0: Бланко (без потерь)')
        ax.plot(freqs, E_M3, '-',  color=COLORS['M3'], lw=2.0, label='M3: +Друде+Рассеяние+τ_par')
        ax.set_title(titles[angle], fontsize=11)
        ax.set_ylabel('|E_пр / E_оп|', fontsize=10)
        ax.set_xlim(F_MIN, F_MAX)
        ax.set_ylim(bottom=0)
        ax.grid(True, alpha=0.3)
        ax.legend(loc='upper right', fontsize=8)

    axes[-1].set_xlabel('Частота (ТГц)', fontsize=10)

    out = OUT_DIR / "pres_slide3_ablation.png"
    fig.savefig(out, dpi=150, bbox_inches='tight', facecolor='white')
    plt.close(fig)
    print("  -> Sokhraneno: " + str(out))


# ============================================================
# SLIDE 4: Fazovaya anizotropiya - 2 grafika vertikal'no
# ============================================================
def plot_slide4_phase_anisotropy():
    print("[Slide4] Fazovaya anizotropiya (40/20, vertikal'no)...")
    try:
        data = load_dataset("test_grid_40_20")
    except Exception as e:
        print("  OSHIBKA: " + str(e))
        return

    angles, freqs, exp = get_spectra(data, F_MIN, F_MAX)

    params_M0 = dict(P_um=40.0, D_um=11.396, loss_factor=0.0, angle_offset=0.0, tau_ps=0.0, gamma=2.0)
    params_M3 = dict(P_um=40.0, D_um=11.370, loss_factor=0.163, angle_offset=0.475, tau_ps=-0.15, gamma=2.0)
    tau_par_M3 = -0.022

    angle_sel = min(angles, key=lambda a: abs(a - 80))
    angles_calc = [angle_sel]

    theo_M0 = theo_grid(angles_calc, freqs, params_M0, use_drude=False, use_scat=False)
    theo_M3 = theo_grid(angles_calc, freqs, params_M3, tau_par=tau_par_M3)

    E_exp = exp[angle_sel]
    E_M0  = theo_M0[angle_sel]
    E_M3  = theo_M3[angle_sel]

    fig, axes = plt.subplots(2, 1, figsize=(6.5, 7.0), sharex=True)
    fig.subplots_adjust(hspace=0.38, left=0.14, right=0.96, top=0.91, bottom=0.09)
    fig.suptitle("Поляризатор 40/20 мкм: фазовая анизотропия (τ_par) при θ ≈ " + str(angle_sel) + "°",
                 fontsize=12, fontweight='bold')

    # Верхний: амплитуда
    ax0 = axes[0]
    ax0.scatter(freqs, np.abs(E_exp), s=8, color=COLORS['exp'], label='Эксперимент', zorder=3, alpha=0.7)
    ax0.plot(freqs, np.abs(E_M0), '--', color=COLORS['M0'], lw=1.8, label='M0: без τ_par')
    ax0.plot(freqs, np.abs(E_M3), '-', color=COLORS['M3'], lw=2.0, label='M3: с τ_par = ' + str(tau_par_M3) + ' пс')
    ax0.set_title('Амплитуда пропускания', fontsize=11)
    ax0.set_ylabel('|E_пр / E_оп|', fontsize=10)
    ax0.set_xlim(F_MIN, F_MAX)
    ax0.set_ylim(bottom=0)
    ax0.grid(True, alpha=0.3)
    ax0.legend(loc='upper right', fontsize=8)

    # Нижний: фаза
    ax1 = axes[1]
    phase_exp = np.unwrap(np.angle(E_exp))
    phase_M0  = np.unwrap(np.angle(E_M0))
    phase_M3  = np.unwrap(np.angle(E_M3))
    ax1.scatter(freqs, phase_exp, s=8, color=COLORS['exp'], label='Эксперимент', zorder=3, alpha=0.7)
    ax1.plot(freqs, phase_M0, '--', color=COLORS['M0'], lw=1.8, label='M0: без τ_par')
    ax1.plot(freqs, phase_M3, '-', color=COLORS['M3'], lw=2.0, label='M3: с τ_par')
    ax1.set_title('Разматанная фаза пропускания', fontsize=11)
    ax1.set_ylabel('Фаза (рад)', fontsize=10)
    ax1.set_xlabel('Частота (ТГц)', fontsize=10)
    ax1.set_xlim(F_MIN, F_MAX)
    ax1.grid(True, alpha=0.3)
    ax1.legend(loc='upper left', fontsize=8)

    out = OUT_DIR / "pres_slide4_tau_par.png"
    fig.savefig(out, dpi=150, bbox_inches='tight', facecolor='white')
    plt.close(fig)
    print("  -> Sokhraneno: " + str(out))


# ============================================================
# SLIDE 5: Validaciya 40/20 - 2 grafika vertikal'no
# ============================================================
def plot_slide5_validation():
    print("[Slide5] Validaciya 40/20 mkm (vertikal'no)...")
    try:
        data = load_dataset("test_grid_40_20")
    except Exception as e:
        print("  OSHIBKA: " + str(e))
        return

    angles, freqs, exp = get_spectra(data, F_MIN, F_MAX)

    params_M3 = dict(P_um=40.0, D_um=11.370, loss_factor=0.163, angle_offset=0.475, tau_ps=-0.15, gamma=2.0)
    tau_par_M3 = -0.022

    angle_0  = min(angles, key=lambda a: abs(a - 0))
    angle_60 = min(angles, key=lambda a: abs(a - 60))
    angles_sel = [angle_0, angle_60]

    theo_M3 = theo_grid(angles_sel, freqs, params_M3, tau_par=tau_par_M3)

    fig, axes = plt.subplots(2, 1, figsize=(6.5, 7.0), sharex=True)
    fig.subplots_adjust(hspace=0.38, left=0.14, right=0.96, top=0.91, bottom=0.09)
    fig.suptitle("Поляризатор 40/20 мкм: подгонка модели M3 к эксперименту",
                 fontsize=12, fontweight='bold')

    titles = {angle_0: "θ ≈ " + str(angle_0) + "° (максимальное пропускание)",
              angle_60: "θ ≈ " + str(angle_60) + "° (промежуточное ослабление)"}

    for ax, angle in zip(axes, angles_sel):
        E_exp = np.abs(exp[angle])
        E_M3  = np.abs(theo_M3[angle])
        ax.scatter(freqs, E_exp, s=8, color=COLORS['exp'], label='Эксперимент', zorder=3, alpha=0.7)
        ax.plot(freqs, E_M3, '-', color=COLORS['M3'], lw=2.0, label='Модель M3')
        ax.fill_between(freqs, E_exp, E_M3, alpha=0.15, color=COLORS['M3'], label='Остатки')
        ax.set_title(titles[angle], fontsize=11)
        ax.set_ylabel('|E_пр / E_оп|', fontsize=10)
        ax.set_xlim(F_MIN, F_MAX)
        ax.set_ylim(bottom=0)
        ax.grid(True, alpha=0.3)
        ax.legend(loc='upper right', fontsize=8)

    axes[-1].set_xlabel('Частота (ТГц)', fontsize=10)

    out = OUT_DIR / "pres_slide5_validation.png"
    fig.savefig(out, dpi=150, bbox_inches='tight', facecolor='white')
    plt.close(fig)
    print("  -> Sokhraneno: " + str(out))


# ============================================================
# SLIDE 6: Zakon masshtabirovaniya D_eff (russkie podpisi)
# ============================================================
def plot_slide6_scaling():
    print("[Slide6] Zakon masshtabirovaniya D_eff...")

    d_over_p = np.linspace(0.1, 0.85, 200)
    ratio_law = 1.0 - 0.85 * d_over_p

    exp_points = [
        (0.50, 0.57, "40/20 мкм\n(D/P=0.50)"),
        (0.71, 0.40, "15.5/11 мкм\n(D/P=0.71)"),
    ]

    fig, ax = plt.subplots(figsize=(6.5, 5.5))
    fig.subplots_adjust(left=0.13, right=0.97, top=0.90, bottom=0.12)

    ax.plot(d_over_p, ratio_law, '--', color='#2980b9', lw=2.0,
            label='Закон: D_эфф/D_физ = 1 - 0.85·(D/P)')
    ax.axhline(1.0, color='gray', lw=1.0, ls=':',
               label='Предел тонкой нити (D_эфф = D_физ)')

    for dp, ratio, label in exp_points:
        ax.scatter([dp], [ratio], s=140, color='#e74c3c', zorder=5)
        ax.annotate(label, xy=(dp, ratio),
                    xytext=(dp + 0.04, ratio + 0.05),
                    fontsize=9, color='#222222',
                    arrowprops=dict(arrowstyle='->', color='#555', lw=1.2),
                    bbox=dict(boxstyle='round,pad=0.3', facecolor='#fffde7',
                              edgecolor='#ccc', alpha=0.9))

    ax.scatter([], [], s=80, color='#e74c3c', label='Экспериментальные WGP')
    ax.set_xlabel('Фактор заполнения D/P', fontsize=11)
    ax.set_ylabel('Отношение D_эфф / D_физ', fontsize=11)
    ax.set_title('Закон масштабирования электродинамического диаметра нити',
                 fontsize=12, fontweight='bold')
    ax.set_xlim(0.1, 0.85)
    ax.set_ylim(0.20, 1.08)
    ax.grid(True, alpha=0.3)
    ax.legend(loc='upper right', fontsize=9)

    out = OUT_DIR / "pres_slide6_deff_law.png"
    fig.savefig(out, dpi=150, bbox_inches='tight', facecolor='white')
    plt.close(fig)
    print("  -> Sokhraneno: " + str(out))


if __name__ == "__main__":
    print("=" * 60)
    print("Generaciya prezentacionnyh grafikov")
    print("=" * 60)
    plot_slide3_ablation()
    plot_slide4_phase_anisotropy()
    plot_slide5_validation()
    plot_slide6_scaling()
    print("\nGotovo!")
