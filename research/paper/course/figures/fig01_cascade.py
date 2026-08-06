# -*- coding: utf-8 -*-
"""Рис. 1.1 — схема измерения на языке матриц Джонса.

Показывает порядок элементов в тракте (поляризатор -> образец под углом theta ->
анализатор -> приёмник) и то, что вращается только образец: вся угловая зависимость
наблюдаемой рождается в нём. Отдельно подчёркнут порядок сомножителей — свет идёт
слева направо, матрицы в формуле стоят справа налево.

Схема, а не расчёт: никаких данных не используется.
"""

import numpy as np

import pstyle as ps
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, Rectangle

fig, ax = plt.subplots(figsize=(7.4, 2.9))
ax.set_xlim(-0.4, 10.4)
ax.set_ylim(-2.05, 1.75)
ax.axis("off")

BEAM_Y = 0.0

# --- луч ---
ax.plot([-0.2, 9.7], [BEAM_Y, BEAM_Y], color=ps.GRID, lw=1.4, zorder=1)
for x in (1.7, 4.0, 6.3, 8.4):
    ax.add_patch(FancyArrowPatch((x, BEAM_Y), (x + 0.35, BEAM_Y), arrowstyle="-|>",
                                 mutation_scale=11, color=ps.GRID, lw=1.4, zorder=1))


def element(x, angle_deg, color, label, sub, rotating=False):
    """Элемент тракта: наклонная штриховка = проволоки, подпись сверху и снизу."""
    w, h = 0.52, 1.05
    ax.add_patch(Rectangle((x - w / 2, BEAM_Y - h / 2), w, h, facecolor="white",
                           edgecolor=color, lw=1.6, zorder=3))
    a = np.deg2rad(angle_deg)
    for s in np.linspace(-0.42, 0.42, 7):
        # штрихи внутри рамки: направление = ориентация проволок
        dx, dy = 0.30 * np.cos(a), 0.30 * np.sin(a)
        px, py = x - s * np.sin(a), BEAM_Y + s * np.cos(a)
        ax.plot([px - dx, px + dx], [py - dy, py + dy], color=color, lw=1.1,
                alpha=0.75, zorder=4, solid_capstyle="butt")
    ax.text(x, BEAM_Y + 0.78, label, ha="center", va="bottom", fontsize=9,
            color=color, fontweight="bold")
    ax.text(x, BEAM_Y - 0.80, sub, ha="center", va="top", fontsize=7.6,
            color=ps.INK_2, linespacing=1.3)
    if rotating:
        # дуга поворота — единственный подвижный элемент схемы
        arc = np.linspace(np.deg2rad(35), np.deg2rad(145), 60)
        ax.plot(x + 0.72 * np.cos(arc), BEAM_Y + 0.72 * np.sin(arc),
                color=ps.S2, lw=1.3, zorder=5)
        ax.add_patch(FancyArrowPatch((x + 0.72 * np.cos(arc[-2]),
                                      BEAM_Y + 0.72 * np.sin(arc[-2])),
                                     (x + 0.72 * np.cos(arc[-1]),
                                      BEAM_Y + 0.72 * np.sin(arc[-1])),
                                     arrowstyle="-|>", mutation_scale=10,
                                     color=ps.S2, lw=1.3, zorder=5))
        ax.text(x, BEAM_Y + 1.34, r"$\theta$", ha="center", va="bottom",
                fontsize=11, color=ps.S2)


# --- источник ---
ax.text(-0.25, BEAM_Y + 0.30, "источник", ha="left", va="bottom", fontsize=8.5,
        color=ps.INK_2)
ax.text(-0.25, BEAM_Y - 0.30, "импульс", ha="left", va="top", fontsize=7.6,
        color=ps.INK_2)

element(2.2, 90, ps.INK_2, "$P_1$", "задаёт входную\nполяризацию")
element(4.6, 55, ps.S2, "образец", "объект исследования\n$\\mathrm{diag}(t_\\perp, t_\\parallel)$",
        rotating=True)
element(6.9, 90, ps.INK_2, "$P_2$", "анализатор")

# --- приёмник ---
ax.add_patch(Rectangle((9.0, BEAM_Y - 0.42), 0.62, 0.84, facecolor=ps.PANEL,
                       edgecolor=ps.S1, lw=1.6, zorder=3))
ax.text(9.31, BEAM_Y + 0.60, "приёмник", ha="center", va="bottom", fontsize=9,
        color=ps.S1, fontweight="bold")
ax.text(9.31, BEAM_Y - 0.62, "регистрирует\nполе $E(t)$", ha="center", va="top",
        fontsize=7.6, color=ps.INK_2, linespacing=1.3)

# --- порядок сомножителей ---
ax.annotate("", xy=(9.2, -1.72), xytext=(0.3, -1.72),
            arrowprops=dict(arrowstyle="-|>", color=ps.S3, lw=1.4))
ax.text(4.75, -1.62, "свет идёт в эту сторону; в формуле матрицы стоят "
                     "в обратном порядке:", ha="center", va="bottom",
        fontsize=8, color=ps.S3)
# Формула набирается обычным текстом, а не mathtext: последний не понимает ни
# \mathsf, ни кириллицу внутри \mathrm, а индексы здесь русские.
ax.text(4.75, -2.02, "E_прм  =  dᵀ · J(P₂) · J(образец, θ) · J(P₁) · E_вх",
        ha="center", va="bottom", fontsize=10, color=ps.INK)

ax.set_title("Вращается только образец — вся угловая зависимость рождается в нём",
             color=ps.TITLE, fontsize=10, loc="left", fontweight="bold", pad=6)

fig.tight_layout()
ps.save(fig, "fig01_cascade")
