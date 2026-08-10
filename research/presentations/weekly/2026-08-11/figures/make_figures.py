#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Графики к недельному отчёту P2 за 04–10 августа 2026.
Источник чисел — research/results/a18_microscopy/SUMMARY.md (A18/A19, 2026-08-10).

⚠ ЭТОТ СКРИПТ НИЧЕГО НЕ ПЕРЕСЧИТЫВАЕТ. Все числа — литералы,
  перенесённые из опубликованного артефакта A18 (SUMMARY.md).
  data_pool не читается, модель не импортируется, фиты не выполняются.
"""
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from pstyle import SURFACE, INK, INK_2, S1, S2, S3, REF, GRID, finish  # noqa

OUT = HERE

# ============================================================
# Рис. 1 — Новые геометрические данные: сравнение с GEOMETRY
# Источник: SUMMARY.md §«Сводная таблица»
# ============================================================

# (образец, P_was, P_new, P_se, D_was, D_new)
# P_was = значение, стоявшее в GEOMETRY до A18
# P_new = результат A18 (граница f=0.50)
GEOM_P = [
    ("purewave",       25.50, 24.44, 0.19),
    ("specac",         24.90, 24.73, 0.22),
    ("test_grid_33_11", 31.82, 31.71, 0.47),
    ("test_grid_40_20", 38.80, 37.18, 1.83),
]

GEOM_D = [
    ("purewave",       10.6,  10.28, 0.014),
    ("specac",         14.0,   9.72, 0.051),  # грубая ошибка — 14.0 → 9.72
    ("test_grid_33_11", 11.2,  10.47, 0.039),
    ("test_grid_40_20", 18.0,  20.30, 0.106),
]


def fig_geometry():
    fig, (ax_p, ax_d) = plt.subplots(1, 2, figsize=(11, 4.5))

    # --- период ---
    idx = np.arange(len(GEOM_P))
    labels = [r[0] for r in GEOM_P]
    was_p = [r[1] for r in GEOM_P]
    new_p = [r[2] for r in GEOM_P]
    err_p = [r[3] for r in GEOM_P]

    ax_p.bar(idx - 0.2, was_p, 0.35, color=GRID, label="было в GEOMETRY", zorder=2)
    ax_p.bar(idx + 0.2, new_p, 0.35, color=S1,   label="A18 (микрофото)", zorder=2)
    ax_p.errorbar(idx + 0.2, new_p, yerr=err_p,
                  fmt="none", color=INK_2, capsize=4, linewidth=1.5, zorder=3)
    ax_p.set_xticks(idx)
    ax_p.set_xticklabels(labels, fontsize=9)
    ax_p.set_ylim(20, 45)
    ax_p.legend(fontsize=9)
    finish(ax_p, "Период P (мкм): было vs стало", "", "P, мкм")

    # --- диаметр ---
    was_d = [r[1] for r in GEOM_D]
    new_d = [r[2] for r in GEOM_D]
    err_d = [r[3] for r in GEOM_D]
    labels_d = [r[0] for r in GEOM_D]

    bars_was = ax_d.bar(idx - 0.2, was_d, 0.35, color=GRID, label="было в GEOMETRY", zorder=2)
    bars_new = ax_d.bar(idx + 0.2, new_d, 0.35, color=S1,   label="A18 (микрофото)", zorder=2)
    ax_d.errorbar(idx + 0.2, new_d, yerr=err_d,
                  fmt="none", color=INK_2, capsize=4, linewidth=1.5, zorder=3)
    # стрелка на ошибку specac
    ax_d.annotate("specac: 14.0 → 9.72\n(ошибка кода исправлена)",
                  xy=(1 + 0.2, 9.72), xytext=(2.3, 13.0),
                  arrowprops=dict(arrowstyle="->", color=S2, lw=1.5),
                  fontsize=8.5, color=S2)
    ax_d.set_xticks(idx)
    ax_d.set_xticklabels(labels_d, fontsize=9)
    ax_d.legend(fontsize=9)
    finish(ax_d, "Диаметр D, мкм (граница на полувысоте f = 0.5)", "", "D, мкм")

    fig.tight_layout(pad=1.5)
    fig.savefig(OUT / "fig1_geometry_update.png", dpi=170)
    plt.close(fig)
    print("  fig1_geometry_update.png  ✓")


# ============================================================
# Рис. 2 — Нерегулярность решётки: σ_pos / P по образцам
# Источник: SUMMARY.md §«Однородность и нерегулярность»
# ============================================================

IRREG = [
    ("purewave",        3.21, 24.44),
    ("specac",          4.69, 24.73),
    ("test_grid_33_11", 5.12, 31.71),
    ("test_grid_40_20", 4.25, 37.18),
]


def fig_irregularity():
    fig, ax = plt.subplots(figsize=(8, 4.2))
    labels = [r[0] for r in IRREG]
    ratio  = [r[1] / r[2] * 100 for r in IRREG]   # σ_pos/P в %
    idx    = np.arange(len(IRREG))

    bars = ax.bar(idx, ratio, 0.5, color=S1, zorder=2)
    for bar, val in zip(bars, ratio):
        ax.text(bar.get_x() + bar.get_width() / 2,
                bar.get_height() + 0.2,
                f"{val:.1f} %",
                ha="center", color=INK, fontsize=10, fontweight="bold")

    ax.axhline(20, color=S2, linewidth=1.5, linestyle="--",
               label="ожидание случайного джиттера (Монте-Карло, 95%)")
    ax.set_xticks(idx)
    ax.set_xticklabels(labels, fontsize=10)
    ax.set_ylim(0, 28)
    ax.legend(fontsize=9)
    finish(ax,
           "Нерегулярность решётки (σ_поз / P)",
           "",
           "стандартное отклонение положения, % от периода")
    fig.tight_layout()
    fig.savefig(OUT / "fig2_irregularity.png", dpi=170)
    plt.close(fig)
    print("  fig2_irregularity.png  ✓")


# ============================================================
# Рис. 3 — D_eff/D_phys: загадка по-прежнему открыта
# Источник чисел: SUMMARY.md + результаты M5 из предыдущей A7
# (BOARD.md строка A, раздел A7, purewave день 1 + 6 образцов A3)
# Значения D_eff/D_phys из артефакта weekly/2026-08-04/SOURCES.md строка Fig.1
# ============================================================

# (образец, D/P физическое, D_eff/D_phys по M5)
# Используем обновлённые D_phys из A18 (f=0.50) и D/P пересчитанные
SAMPLES = [
    ("att-11-16-356", 0.710,  0.373),    # D/P паспорт (att нет микрофото)
    ("att-11-16-s1",  0.688,  0.528),
    ("att-11-16-s2",  0.688,  0.452),
    ("att-11-16-s3",  0.688,  0.371),
    ("specac",        0.393,  0.443),    # D/P = D_new/P_new = 9.72/24.73 ≈ 0.393
    ("test_grid_40_20", 0.546, 0.586),   # D_new/P_new = 20.30/37.18
    ("purewave (д.1)", 0.421, 0.659),   # D_new/P_new = 10.28/24.44
]
# Закон H5: D_eff/D_phys ≈ 1 - 0.85·(D/P)  — из HYPOTHESES.md строка H5


def fig_deff():
    xs_emp = np.linspace(0.30, 0.75, 200)
    ys_h5  = 1 - 0.85 * xs_emp

    fig, ax = plt.subplots(figsize=(8.5, 4.8))
    ax.plot(xs_emp, ys_h5, color=REF, linewidth=2.0, linestyle="--",
            label="эмпирический закон H5:  1 − 0.85·D/P")
    ax.axhline(1.0, color=GRID, linewidth=1.5, linestyle=":",
               label="точный расчёт по решётке (эталон ≈ 1.0)")

    ax.plot([s[1] for s in SAMPLES], [s[2] for s in SAMPLES],
            "o", markersize=9, color=S1, markeredgecolor="white",
            markeredgewidth=1.5, linestyle="none",
            label="извлечено из данных (фит M5)")

    LABEL_OFFSETS = {
        "att-11-16-356": (0,  13, "center"),
        "att-11-16-s1":  (-12, 10, "right"),
        "att-11-16-s2":  (-12, -8, "right"),
        "att-11-16-s3":  (-12, -20, "right"),
        "specac":         (14,  -6, "left"),
        "test_grid_40_20": (0,  13, "center"),
        "purewave (д.1)": (0,  13, "center"),
    }
    for name, xv, yv in SAMPLES:
        dx, dy, ha = LABEL_OFFSETS[name]
        ax.annotate(name, (xv, yv), textcoords="offset points",
                    xytext=(dx, dy), ha=ha, color=INK_2, fontsize=8.5)

    ax.set_xlim(0.28, 0.78)
    ax.set_ylim(0.20, 1.15)
    ax.legend(fontsize=9)
    finish(ax,
           "Загадка D_eff: данные требуют «тонких» проволок, расчёт — нет",
           "D/P физическое (по микрофото)",
           "D_eff / D_phys")
    fig.tight_layout()
    fig.savefig(OUT / "fig3_deff_puzzle.png", dpi=170)
    plt.close(fig)
    print("  fig3_deff_puzzle.png  ✓")


if __name__ == "__main__":
    print("Рисуем графики для недельного отчёта P2 (2026-08-11)…")
    fig_geometry()
    fig_irregularity()
    fig_deff()
    print("Готово: fig1, fig2, fig3")
