#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Графики-сводки к недельному отчёту зоны P за 2026-07-28…2026-08-04.

ЭТОТ СКРИПТ НИЧЕГО НЕ СЧИТАЕТ. Все числа ниже — литералы, перенесённые вручную из уже
опубликованных артефактов сессии A; у каждого блока указан источник. `data_pool/` не
читается, модель не импортируется. Графики, построенные ПО ДАННЫМ (наложение теории и
эксперимента, микрофото дифракции), вынесены в отдельные скрипты `make_overlay.py` и
`make_micrographs.py` — они работают по прямому разрешению владельца от 2026-08-04.

Имена наборов — по `data_pool/DATASETS.md` (переименование 2026-08-04).

Запуск из корня репозитория:
    .venv\\Scripts\\python.exe research/presentations/weekly/2026-08-04/figures/make_figures.py
"""

import sys
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt

sys.path.insert(0, str(Path(__file__).resolve().parent))
from pstyle import SURFACE, INK, INK_2, S1, S2, S3, REF, finish  # noqa: E402

OUT = Path(__file__).resolve().parent


# ---------------------------------------------------------------------------
# Рис. 1 — загадка D_eff: семь наборов против точного эталона
# Источник: BRIEF_FOR_P §3.1 (он же LOG, запись «[A] 2026-08-04 — A7», раздел 4).
# Значение эталона ≈1.0 — SYNTHESIS.md §3 (T5b).
# ---------------------------------------------------------------------------
SAMPLES = [
    # (имя набора по DATASETS.md, D/P, D_eff/D_phys, подпись на графике — паспортные D/P)
    # Источник паспортов: two_wgp_attenuator/README.md, SYNTHESIS.md
    ("att-11-16-356", 0.710, 0.373, "att-11-16\n(D=11, P=16 пасп.)"),
    ("att-11-16-s1",  0.688, 0.528, "att-11-16-s1\n(D=11, P=16 пасп.)"),
    ("att-11-16-s2",  0.688, 0.452, "att-11-16-s2\n(D=11, P=16 пасп.)"),
    ("att-11-16-s3",  0.688, 0.371, "att-11-16-s3\n(D=11, P=16 пасп.)"),
    ("specac",        0.562, 0.443, "specac\n(D≈14, P=24.9)"),
    ("test_grid_40_20", 0.464, 0.586, "test_grid\n(D=20, P=40 пасп.)"),
    ("purewave",      0.416, 0.659, "purewave\n(D=10.6, P=25.5)"),
]
LABEL_OFFSETS = {          # подобраны вручную, чтобы подписи не сталкивались
    "att-11-16-356":   (  0,  14, "center"),
    "att-11-16-s1":    (-14,  10, "right"),
    "att-11-16-s2":    (-14,  -6, "right"),
    "att-11-16-s3":    (-14, -24, "right"),
    "specac":          ( 16,  -5, "left"),
    "test_grid_40_20": (  0,  14, "center"),
    "purewave":        (  0,  14, "center"),
}


def fig_deff():
    fig, ax = plt.subplots(figsize=(9.0, 5.0))

    # Опорный эталон: точный lattice-sum даёт ≈1.0
    ax.axhline(1.0, color=REF, linewidth=2.0, linestyle="--",
               label=r"точный расчёт по решётке (lattice-sum) $\approx$ 1.0")

    # Экспериментальные точки
    ax.plot([s[1] for s in SAMPLES], [s[2] for s in SAMPLES], "o", markersize=11,
            color=S1, markeredgecolor=SURFACE, markeredgewidth=1.5, linestyle="none",
            label="извлечено из данных (фит M5, 7 наборов)")

    # Подписи: «имя (паспортный D/P)»
    for name, xv, yv, lbl in SAMPLES:
        dx, dy, ha = LABEL_OFFSETS[name]
        ax.annotate(lbl, (xv, yv), textcoords="offset points",
                    xytext=(dx, dy), ha=ha, color=INK_2, fontsize=8.5,
                    linespacing=1.25)


    ax.set_ylim(0.25, 1.18)
    ax.set_xlim(0.38, 0.78)
    finish(ax, r"Загадка $D_{eff}$: данные требуют в 1.5–2.7× тоньше реальных проволок",
           r"плотность решётки  $D_{{phys}}/P$  (физическая геометрия)",
           r"$D_{eff}/D_{phys}$")
    ax.legend(loc="upper right", frameon=False, fontsize=10)
    fig.tight_layout()
    fig.savefig(OUT / "fig1_deff_puzzle.png", dpi=170)
    plt.close(fig)


# ---------------------------------------------------------------------------
# Рис. 2 — что воспроизводимо: суммарная утечка η, а не её разложение
# Источник: BRIEF_FOR_P §3.2 (он же LOG, «2026-08-04 · [A] A3_refit_series», «Числа» п.2).
# ---------------------------------------------------------------------------
LEAK = [
    # (имя набора, η(A0) модельно-незав., η₀(M5), |ε_cross|)
    # Источник: BRIEF_FOR_P §3.2 + a7_purewave.json (harmonics.day1_-95_80, fit.day1_-95_80.best)
    ("att-11-16-356", 0.01266, 0.00115, 0.00664),
    ("att-11-16-s1",  0.01210, 0.06593, 0.04746),
    ("att-11-16-s2",  0.01298, 0.04235, 0.04175),
    ("att-11-16-s3",  0.01537, 0.00277, 0.00911),
    ("specac",        0.03594, 0.00039, 0.00330),
    ("test_grid_40_20", 0.03290, 0.00971, 0.01023),
    ("purewave",      0.01847, 0.01119, 0.01710),  # a7, день 1 (-95..80°)
]


def fig_leakage():
    fig, ax = plt.subplots(figsize=(9.0, 4.8))
    idx = np.arange(len(LEAK))

    ax.plot(idx, [r[1] for r in LEAK], "o-", color=S1, markersize=10,
            markeredgecolor=SURFACE, markeredgewidth=1.5,
            label=r"$\eta$ — измерено без модели")
    ax.plot(idx, [r[2] for r in LEAK], "s--", color=S2, markersize=9,
            markeredgecolor=SURFACE, markeredgewidth=1.5,
            label=r"$\eta_0$ — параметр модели M5")
    ax.plot(idx, [r[3] for r in LEAK], "^--", color=S3, markersize=9,
            markeredgecolor=SURFACE, markeredgewidth=1.5,
            label=r"$|\varepsilon_{cross}|$ — параметр модели M5")

    ax.axvline(3.5, color=REF, linewidth=1.2, linestyle=":")
    ax.text(1.5, 1.35e-4, "один прибор ATT-11-16\n(четыре съёмки)",
            ha="center", va="bottom", color=INK_2, fontsize=9.0)
    ax.text(5.0, 1.35e-4, "одиночный WGP\n(три образца)",
            ha="center", va="bottom", color=INK_2, fontsize=9.0)

    ax.set_yscale("log")
    ax.set_ylim(1e-4, 0.5)
    ax.set_xlim(-0.4, len(LEAK) - 0.6)
    ax.set_xticks(idx)
    ax.set_xticklabels([r[0] for r in LEAK], rotation=15, ha="right", fontsize=9)
    finish(ax, "Что воспроизводится от съёмки к съёмке, а что нет",
           "", "утечка (доля амплитуды, логарифмическая шкала)")
    ax.legend(loc="upper left", frameon=False, fontsize=10)
    fig.tight_layout()
    fig.savefig(OUT / "fig2_leakage_reproducibility.png", dpi=170)
    plt.close(fig)


# ---------------------------------------------------------------------------
# Рис. 3 — цена смешивания двух дней съёмки purewave
# Источник: BRIEF_FOR_P §3.4 (он же LOG, «[A] 2026-08-04 — A7», разделы 2 и 3).
# ---------------------------------------------------------------------------
RANGES = [
    ("день 1\n(−95…80°)", 0.490, S1),
    ("день 2\n(82…110°)", 0.801, S1),
    ("оба дня\n(−95…110°)", 1.510, S2),
    ("продакшн-обрезка\n(0…90°)", 2.211, S2),
]


def fig_purewave_days():
    fig, ax = plt.subplots(figsize=(8.4, 4.2))
    idx = np.arange(len(RANGES))
    vals = [r[1] for r in RANGES]
    bars = ax.bar(idx, vals, width=0.58, color=[r[2] for r in RANGES],
                  edgecolor=SURFACE, linewidth=1.5)
    for b, v in zip(bars, vals):
        ax.text(b.get_x() + b.get_width() / 2, v + 0.07, f"{v:.3f}",
                ha="center", color=INK, fontsize=11, fontweight="bold")

    ax.set_xticks(idx)
    ax.set_xticklabels([r[0] for r in RANGES], fontsize=9.5)
    ax.set_ylim(0, 2.6)
    finish(ax, "Невязка модели в тени: внутри дня и при смешивании дней",
           "", "DS_rmse (дБ)")
    fig.tight_layout()
    fig.savefig(OUT / "fig3_purewave_two_days.png", dpi=170)
    plt.close(fig)


if __name__ == "__main__":
    fig_deff()
    fig_leakage()
    fig_purewave_days()
    print("готово: fig1, fig2, fig3")
