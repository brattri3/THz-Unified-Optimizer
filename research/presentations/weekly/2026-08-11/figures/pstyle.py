# -*- coding: utf-8 -*-
"""
Палитра зоны P — светлая тема.
Решение владельца 2026-08-04: все отчёты P строятся на БЕЛОМ фоне.
Правила — STYLE_MANIFEST.md §4.
"""
from pathlib import Path
import matplotlib.pyplot as plt

# --- Поверхность и текст ---
SURFACE = "#ffffff"
PANEL   = "#f1f5f9"
INK     = "#0f172a"
INK_2   = "#64748b"
GRID    = "#cbd5e1"
REF     = "#94a3b8"

# --- Категориальные слоты (НЕ менять порядок) ---
S1 = "#2a78d6"   # слот 1 — основная измеренная величина (синий)
S2 = "#eb6834"   # слот 2 — сравниваемая/конкурирующая (оранжевый)
S3 = "#1a9e6c"   # слот 3 — третья серия (бирюзовый)

plt.rcParams.update({
    "figure.facecolor": SURFACE,
    "axes.facecolor":   SURFACE,
    "axes.edgecolor":   INK_2,
    "xtick.color":      INK_2,
    "ytick.color":      INK_2,
    "text.color":       INK,
    "axes.labelcolor":  INK,
    "font.family":      "sans-serif",
    "font.size":        11,
    "axes.titlesize":   13,
    "axes.titleweight": "bold",
    "lines.linewidth":  2.0,
    "axes.grid":        True,
    "grid.color":       GRID,
    "grid.linewidth":   0.6,
    "legend.frameon":   False,
})


def finish(ax, title="", xlabel="", ylabel=""):
    """Стандартное завершение осей: заголовок, подписи, скрытые рёбра."""
    ax.set_title(title, color=INK, loc="left", pad=10)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
