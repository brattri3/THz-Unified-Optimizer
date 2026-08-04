#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Графики к недельному отчёту зоны P за 2026-07-28…2026-08-04.

ВАЖНО (гардрейл зоны P, `research/presentations/STATE.md` п.3):
этот скрипт НИЧЕГО НЕ СЧИТАЕТ. Все числа ниже — литералы, перенесённые вручную из
уже опубликованных артефактов сессии A; у каждого блока указан источник. Скрипт не
читает `data_pool/`, не импортирует модель и не выполняет фитов. Единственная
арифметика — эмпирический закон H5 `1 − 0.85·D/P`, который рисуется как опорная
кривая (его формула опубликована в `research/hypotheses/HYPOTHESES.md`, строка H5).

Запуск из корня репозитория:
    .venv\\Scripts\\python.exe research/presentations/weekly/2026-08-04/figures/make_figures.py
"""

from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

OUT = Path(__file__).resolve().parent

# --- Палитра (research/presentations/STYLE_MANIFEST.md §4, прогнана валидатором dataviz) ---
SURFACE = "#0f172a"
INK = "#f8fafc"
INK_2 = "#94a3b8"
GRID = "#334155"
S1 = "#3987e5"   # слот 1 — основная измеренная величина
S2 = "#d95926"   # слот 2 — конкурирующая величина
S3 = "#199e70"   # слот 3 — третья серия
REF = "#94a3b8"  # опорные линии/законы/эталоны — НЕ серия данных

plt.rcParams.update({
    "figure.facecolor": SURFACE,
    "axes.facecolor": SURFACE,
    "savefig.facecolor": SURFACE,
    "text.color": INK,
    "axes.labelcolor": INK,
    "axes.edgecolor": GRID,
    "xtick.color": INK_2,
    "ytick.color": INK_2,
    "grid.color": GRID,
    "grid.alpha": 0.45,
    "font.family": "DejaVu Sans",
    "font.size": 11,
    "axes.titlesize": 13,
    "lines.linewidth": 2.0,
})


def _finish(ax, title, xlabel, ylabel):
    ax.set_title(title, color="#38bdf8", pad=14, loc="left", fontweight="bold")
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.grid(True, linewidth=0.8)
    ax.set_axisbelow(True)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)


# ---------------------------------------------------------------------------
# Рис. 1 — загадка D_eff: 7 образцов против точного эталона
# Источник чисел: RESEARCH_LOG.md, запись «[A] 2026-08-04 — A7», раздел 4 (таблица
#   D/P · D_eff/D_phys(полн.) · H5). Значение эталона ≈1.0 — SYNTHESIS.md §3 (T5b).
# ---------------------------------------------------------------------------
SAMPLES = [
    # (имя, D/P, D_eff/D_phys на ПОЛНОМ угловом диапазоне)
    ("att-11-16-356", 0.710, 0.373),
    ("att-11-16-s1",  0.688, 0.528),
    ("att-11-16-s2",  0.688, 0.452),
    ("att-11-16-s3",  0.688, 0.371),
    ("specac",        0.562, 0.443),
    ("test_grid_40_20", 0.464, 0.586),
    ("purewave (день 1)", 0.416, 0.659),
]


def fig_deff():
    fig, ax = plt.subplots(figsize=(9.0, 5.0))

    x = np.linspace(0.38, 0.75, 100)
    ax.plot(x, 1.0 - 0.85 * x, color=REF, linewidth=2.0, linestyle="-",
            label="эмпирический закон H5:  1 − 0.85·D/P")
    ax.axhline(1.0, color=REF, linewidth=2.0, linestyle="--",
               label="точный lattice-sum (эталон):  ≈ 1.0")

    xs = [s[1] for s in SAMPLES]
    ys = [s[2] for s in SAMPLES]
    ax.plot(xs, ys, "o", markersize=10, color=S1, markeredgecolor=SURFACE,
            markeredgewidth=2.0, linestyle="none", label="фит M5, полный угловой диапазон")

    # смещения подписей подобраны вручную, чтобы не сталкивались у плотной группы
    label_offsets = {
        "att-11-16-356": (0, 14, "center"),
        "att-11-16-s1": (-10, 12, "right"),
        "att-11-16-s2": (-10, -6, "right"),
        "att-11-16-s3": (-10, -20, "right"),
        "specac": (14, -6, "left"),
        "test_grid_40_20": (0, 13, "center"),
        "purewave (день 1)": (0, 13, "center"),
    }
    for name, xv, yv in SAMPLES:
        dx, dy, ha = label_offsets[name]
        ax.annotate(name, (xv, yv), textcoords="offset points",
                    xytext=(dx, dy), ha=ha, color=INK_2, fontsize=9)

    ax.annotate("разрыв ≈ 0.34…0.63\n— это и есть загадка $D_{eff}$",
                xy=(0.485, 0.83), color=INK, fontsize=11, ha="center",
                bbox=dict(boxstyle="round,pad=0.45", facecolor="#1e293b",
                          edgecolor=GRID))
    ax.annotate("", xy=(0.562, 1.0), xytext=(0.562, 0.443),
                arrowprops=dict(arrowstyle="<->", color=INK_2, linewidth=1.5))

    ax.set_ylim(0.25, 1.12)
    ax.set_xlim(0.38, 0.76)
    _finish(ax,
            "Загадка $D_{eff}$: данные требуют «тонких» проволок, точный расчёт — нет",
            "плотность решётки  D/P  (физическая геометрия)",
            "$D_{eff}/D_{phys}$")
    ax.legend(loc="upper right", bbox_to_anchor=(1.0, 0.93), frameon=False,
              labelcolor=INK, fontsize=10)
    fig.tight_layout()
    fig.savefig(OUT / "fig1_deff_puzzle.png", dpi=170)
    plt.close(fig)


# ---------------------------------------------------------------------------
# Рис. 2 — что воспроизводимо: суммарная утечка η, а не её разложение
# Источник чисел: RESEARCH_LOG.md, запись «2026-08-04 · [A] A3_refit_series»,
#   раздел «Числа», п.2 (таблица η(A0) · η₀(M5) · |ε_cross| · D_eff/D_phys).
# Имена датасетов даны в старой нотации, как в самой записи журнала A3;
# соответствие новым именам — data_pool/two_wgp_attenuator/README.md.
# ---------------------------------------------------------------------------
LEAK = [
    # (датасет, группа, η(A0) модельно-незав., η₀(M5), |ε_cross|)
    ("356att",          "аттенюатор", 0.01266, 0.00115, 0.00664),
    ("series1",         "аттенюатор", 0.01210, 0.06593, 0.04746),
    ("series2",         "аттенюатор", 0.01298, 0.04235, 0.04175),
    ("series3",         "аттенюатор", 0.01537, 0.00277, 0.00911),
    ("specac",          "одиночный",  0.03594, 0.00039, 0.00330),
    ("test_grid_40_20", "одиночный",  0.03290, 0.00971, 0.01023),
]


def fig_leakage():
    fig, ax = plt.subplots(figsize=(9.0, 5.0))
    idx = np.arange(len(LEAK))

    ax.plot(idx, [r[2] for r in LEAK], "o-", color=S1, markersize=10,
            markeredgecolor=SURFACE, markeredgewidth=2.0,
            label=r"$\eta$ — гармоники A0, модельно-независимо  ·  размах в группе $\times$1.3")
    ax.plot(idx, [r[3] for r in LEAK], "s--", color=S2, markersize=9,
            markeredgecolor=SURFACE, markeredgewidth=2.0,
            label=r"$\eta_0$ — параметр модели M5  ·  размах $\times$57")
    ax.plot(idx, [r[4] for r in LEAK], "^--", color=S3, markersize=9,
            markeredgecolor=SURFACE, markeredgewidth=2.0,
            label=r"$|\varepsilon_{cross}|$ — параметр модели M5  ·  размах $\times$7.1")

    ax.axvline(3.5, color=GRID, linewidth=1.5, linestyle=":")
    ax.text(1.5, 1.35e-4, "группа «аттенюатор» — один прибор ATT-11-16",
            ha="center", color=INK_2, fontsize=9)
    ax.text(4.5, 1.35e-4, "группа «одиночный WGP»", ha="center",
            color=INK_2, fontsize=9)

    ax.set_yscale("log")
    ax.set_ylim(1e-4, 0.5)
    ax.set_xlim(-0.4, len(LEAK) - 0.6)
    ax.set_xticks(idx)
    ax.set_xticklabels([r[0] for r in LEAK], rotation=15, ha="right", fontsize=9)
    _finish(ax,
            "Воспроизводима суммарная утечка, а не её разложение",
            "", "утечка (доля амплитуды, лог. шкала)")
    ax.legend(loc="upper left", frameon=False, labelcolor=INK, fontsize=9.5, ncol=1)
    fig.tight_layout()
    fig.savefig(OUT / "fig2_leakage_reproducibility.png", dpi=170)
    plt.close(fig)


# ---------------------------------------------------------------------------
# Рис. 3 — цена смешивания двух дней съёмки purewave
# Источник чисел: RESEARCH_LOG.md, запись «[A] 2026-08-04 — A7», разделы 2 и 3
#   (DS_rmse по диапазонам: полный 1.510 / день 1 0.490 / день 2 0.801 /
#    продакшн (0,90) 2.211 дБ).
# ---------------------------------------------------------------------------
RANGES = [
    ("день 1\n(−95…80°)", 0.490, S1),
    ("день 2\n(82…110°)", 0.801, S1),
    ("оба дня\n(−95…110°)", 1.510, S2),
    ("продакшн-обрезка\n(0…90°), смешивает дни", 2.211, S2),
]


def fig_purewave():
    fig, ax = plt.subplots(figsize=(9.0, 4.6))
    idx = np.arange(len(RANGES))
    vals = [r[1] for r in RANGES]
    cols = [r[2] for r in RANGES]

    bars = ax.bar(idx, vals, width=0.58, color=cols, edgecolor=SURFACE, linewidth=2.0)
    for b, v in zip(bars, vals):
        ax.text(b.get_x() + b.get_width() / 2, v + 0.06, f"{v:.3f} дБ",
                ha="center", color=INK, fontsize=11, fontweight="bold")

    ax.set_xticks(idx)
    ax.set_xticklabels([r[0] for r in RANGES], fontsize=9.5)
    ax.set_ylim(0, 3.1)
    _finish(ax,
            "Один образец, два дня: внутри дня согласовано, между днями — нет",
            "", "невязка в глубокой тени, DS_rmse (дБ)")
    ax.text(0.5, 2.72, "внутри одного дня (синий) — норма", ha="center",
            color=INK_2, fontsize=10)
    ax.text(2.5, 2.82, "как только дни смешиваются (оранжевый) —\nневязка растёт в 3–4.5 раза",
            ha="center", color=INK_2, fontsize=10)
    fig.tight_layout()
    fig.savefig(OUT / "fig3_purewave_two_days.png", dpi=170)
    plt.close(fig)


if __name__ == "__main__":
    fig_deff()
    fig_leakage()
    fig_purewave()
    print("готово:", ", ".join(sorted(p.name for p in OUT.glob("*.png"))))
