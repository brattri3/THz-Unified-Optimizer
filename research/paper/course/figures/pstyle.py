# -*- coding: utf-8 -*-
"""Единый стиль графиков учебного курса (зона D).

Основа скопирована из зоны P (`research/presentations/weekly/2026-08-04/figures/pstyle.py`,
палитра прогнана валидатором dataviz: светлотная полоса / цветность / дальтонизм — PASS) и
адаптирована под ПЕЧАТНЫЙ документ:
  * выход в PDF (векторный), а не PNG — курс печатный, увеличение не должно пикселизовать;
  * шрифт мельче и линии тоньше, чем на слайдах: рисунок в тексте A4 меньше слайда;
  * шрифты вкладываются в PDF как TrueType (pdf.fonttype=42) — иначе часть просмотрщиков
    подставляет свои и кириллица в подписях ломается.

Палитра совпадает с цветами преамбулы LaTeX (`preamble.tex`), чтобы текст и иллюстрации
читались как один документ.
"""

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

SURFACE = "#ffffff"   # поверхность графика = фон страницы
INK = "#0f172a"       # основной текст
INK_2 = "#475569"     # вторичный текст, подписи точек
GRID = "#cbd5e1"      # сетка
TITLE = "#0b5aa8"     # заголовок графика (тон заголовков документа)
PANEL = "#f1f5f9"     # заливка врезок

S1 = "#2a78d6"        # слот 1 — основная величина
S2 = "#eb6834"        # слот 2 — конкурирующая величина
S3 = "#1baf7a"        # слот 3 — третья серия (всегда с прямой подписью)
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
    "font.size": 9,
    "axes.titlesize": 10.5,
    "legend.fontsize": 8.5,
    "lines.linewidth": 1.6,
    "pdf.fonttype": 42,
    "savefig.bbox": "tight",
    "savefig.pad_inches": 0.02,
})


def finish(ax, title=None, xlabel=None, ylabel=None, legend=False):
    """Общая отделка осей: заголовок слева, приглушённая сетка позади данных."""
    if title:
        ax.set_title(title, color=TITLE, pad=10, loc="left", fontweight="bold")
    if xlabel:
        ax.set_xlabel(xlabel)
    if ylabel:
        ax.set_ylabel(ylabel)
    ax.grid(True, linewidth=0.6)
    ax.set_axisbelow(True)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    if legend:
        leg = ax.legend(frameon=False, loc=legend if isinstance(legend, str) else "best")
        for text in leg.get_texts():
            text.set_color(INK)


def save(fig, name):
    """Сохранить рядом со скриптом под именем `name`.pdf и сообщить путь."""
    out = Path(__file__).resolve().parent / f"{name}.pdf"
    fig.savefig(out)
    plt.close(fig)
    print(f"[figures] {out.name}")
    return out
