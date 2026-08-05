# -*- coding: utf-8 -*-
"""Рис. 0.3 — геометрия решётки и конвенция углов проекта.

Показывает то место, где ошибаются чаще всего: проволочный поляризатор пропускает поле,
перпендикулярное проволокам, а не параллельное им (интуиция «щель в заборе» даёт
противоположный ответ). Плюс фиксирует конвенцию репозитория: угол элемента = ориентация
оси ПРОПУСКАНИЯ, проволоки лежат под углом +90 градусов к ней (research/two_wgp/model_2wgp.py).

Схема, а не расчёт: никаких данных не используется.
"""

import numpy as np

import pstyle as ps
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch

fig, (axA, axB) = plt.subplots(1, 2, figsize=(7.4, 3.6))

# =========================== панель A: два канала ===========================
axA.set_xlim(-1.45, 1.45)
axA.set_ylim(-1.55, 1.45)
axA.set_aspect("equal")
axA.axis("off")

# Проволоки укорочены сверху и снизу, чтобы подписи каналов легли на чистое поле,
# а не поверх штриховки: наложение текста на решётку делало рисунок нечитаемым.
for x in np.linspace(-1.0, 1.0, 11):
    axA.plot([x, x], [-0.35, 0.95], color=ps.INK_2, lw=3.0, solid_capstyle="butt")
axA.text(0, 1.10, "проволоки", ha="center", fontsize=9, color=ps.INK_2)

BOX = dict(facecolor="white", edgecolor="none", pad=1.2)

# поле ВДОЛЬ проволок -> тёмный канал
axA.add_patch(FancyArrowPatch((-0.62, 0.05), (-0.62, 0.72),
                              arrowstyle="-|>", mutation_scale=13,
                              color=ps.S2, lw=2.2, zorder=4))
axA.text(-0.62, 0.38, "поле $\\parallel$", ha="center", va="center", fontsize=9,
         color=ps.S2, zorder=5, bbox=BOX)
axA.text(-0.62, -0.55, "ток течёт\n$\\Rightarrow$ отражение\n(тёмный канал $t_\\parallel$)",
         ha="center", va="top", fontsize=8, color=ps.S2)

# поле ПОПЕРЁК проволок -> яркий канал
axA.add_patch(FancyArrowPatch((0.20, 0.38), (1.05, 0.38),
                              arrowstyle="-|>", mutation_scale=13,
                              color=ps.S3, lw=2.2, zorder=4))
axA.text(0.62, 0.38, "поле $\\perp$", ha="center", va="center", fontsize=9,
         color=ps.S3, zorder=5, bbox=BOX)
axA.text(0.62, -0.55, "току негде течь\n$\\Rightarrow$ проходит\n(яркий канал $t_\\perp$)",
         ha="center", va="top", fontsize=8, color=ps.S3)

axA.set_title("Пропускается то, что ПОПЕРЁК проволок",
              color=ps.TITLE, fontsize=9.5, loc="left", fontweight="bold", pad=6)

# ======================= панель B: конвенция углов =========================
axB.set_xlim(-1.55, 1.75)
axB.set_ylim(-1.55, 1.45)
axB.set_aspect("equal")
axB.axis("off")

theta = np.deg2rad(25.0)

# проволоки под углом theta + 90
n = np.array([np.cos(theta), np.sin(theta)])          # ось пропускания
w = np.array([-np.sin(theta), np.cos(theta)])         # направление проволок
for s in np.linspace(-0.95, 0.95, 11):
    p0 = n * s - w * 0.95
    p1 = n * s + w * 0.95
    axB.plot([p0[0], p1[0]], [p0[1], p1[1]], color=ps.GRID, lw=2.4,
             solid_capstyle="butt", zorder=1)

# оси
axB.add_patch(FancyArrowPatch((0, 0), tuple(1.12 * n), arrowstyle="-|>",
                              mutation_scale=13, color=ps.S1, lw=2.2, zorder=3))
axB.text(*(1.20 * n), r"ось пропускания $=\theta$", color=ps.S1, fontsize=9,
         ha="left", va="center")

axB.add_patch(FancyArrowPatch((0, 0), tuple(0.85 * w), arrowstyle="-|>",
                              mutation_scale=11, color=ps.INK_2, lw=1.6, zorder=3))
axB.text(*(0.95 * w), r"проволоки $=\theta+90^\circ$", color=ps.INK_2, fontsize=9,
         ha="center", va="bottom")

# опорная ось x и дуга угла
axB.plot([0, 1.15], [0, 0], color=ps.GRID, lw=1.0, zorder=2)
arc = np.linspace(0, theta, 60)
axB.plot(0.45 * np.cos(arc), 0.45 * np.sin(arc), color=ps.S1, lw=1.2, zorder=3)
axB.text(0.55 * np.cos(theta / 2), 0.55 * np.sin(theta / 2), r"$\theta$",
         color=ps.S1, fontsize=11, ha="left", va="center")

axB.set_title("Конвенция репозитория: угол = ось пропускания", color=ps.TITLE,
              fontsize=9.5, loc="left", fontweight="bold", pad=6)

fig.tight_layout(w_pad=3.0)
ps.save(fig, "fig00_geometry")
