# -*- coding: utf-8 -*-
"""Единый стиль графиков зоны P — СВЕТЛАЯ тема.

Палитра и правила — `research/presentations/STYLE_MANIFEST.md §4`. Значения прогнаны
валидатором навыка `dataviz` на белой поверхности:
    validate_palette.js "#2a78d6,#eb6834,#1baf7a" --mode light --surface "#ffffff" --pairs all
    → светлотная полоса PASS, цветность PASS, дальтонизм (худшая пара ΔE 9.2) PASS,
      нормальное зрение (ΔE 24.0) PASS, контраст: слот 3 = 2.82:1 → действует «правило
      рельефа»: у слота 3 всегда есть прямая подпись или строка в таблице слайда.
"""

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

SURFACE = "#ffffff"   # поверхность графика (совпадает с фоном слайда)
INK = "#0f172a"       # основной текст
INK_2 = "#475569"     # вторичный текст, подписи точек
GRID = "#cbd5e1"      # сетка
TITLE = "#0b5aa8"     # заголовок графика (тот же тон, что заголовки слайдов)
PANEL = "#f1f5f9"     # заливка врезок

S1 = "#2a78d6"        # слот 1 — основная измеренная величина
S2 = "#eb6834"        # слот 2 — конкурирующая величина
S3 = "#1baf7a"        # слот 3 — третья серия (всегда с прямой подписью, см. выше)
REF = "#64748b"       # опорные линии, законы, эталоны — НЕ серия данных

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
    "grid.alpha": 0.9,
    "font.family": "DejaVu Sans",
    "font.size": 11,
    "axes.titlesize": 13,
    "lines.linewidth": 2.0,
})


def finish(ax, title, xlabel, ylabel):
    """Общая отделка осей: заголовок слева, приглушённая сетка позади данных."""
    ax.set_title(title, color=TITLE, pad=14, loc="left", fontweight="bold")
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.grid(True, linewidth=0.8)
    ax.set_axisbelow(True)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
