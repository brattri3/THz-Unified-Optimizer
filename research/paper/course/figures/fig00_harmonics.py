# -*- coding: utf-8 -*-
"""Рис. 0.2 — почему измерение поля даёт пять чисел, а не три.

Левая панель: тождество cos^4(theta) = 3/8 + 1/2 cos2theta + 1/8 cos4theta —
идеальная угловая кривая раскладывается ровно на постоянную, 2-ю и 4-ю гармоники.
Правая панель: отношение амплитуд A4/A2. Для идеального cos^4 оно равно ровно 0.25;
измеренные значения по четырём наборам отклоняются менее чем на 2 %.

Чистая тригонометрия + УЖЕ ПОСЧИТАННЫЕ числа (ничего не пересчитывается):
A4/A2 = 0.2498 / 0.2458 / 0.2549 / 0.2549 — BRIEF_FOR_D_teaching_ladder.md §5.2,
артефакт research/results/two_wgp/a0_malus_harmonics_data.json.

ЖАНР (STYLE.md §3): на подписях осей образцы обозначены физически — через параметр
плотности D/P, а не внутренними именами наборов. Соответствие «подпись -> набор»
живёт в приложении «Происхождение числовых примеров».
"""

import numpy as np

import pstyle as ps
import matplotlib.pyplot as plt

# --- измеренные отношения A4/A2 (перенесены, не пересчитаны) ---
# att-11-16-356 -> D/P = 0.71; att-11-16-s1, -s2 -> D/P = 0.69 (повторные
# установки ОДНОГО прибора); test_grid_40_20 -> D/P = 0.46.
MEASURED = [
    ("плотная решётка\n$D/P = 0.71$", 0.2498),
    ("повторная установка\nтого же образца (1)", 0.2458),
    ("повторная установка\nтого же образца (2)", 0.2549),
    ("разреженный образец\n$D/P = 0.46$", 0.2549),
]
IDEAL = 0.25

th = np.linspace(-90, 90, 721)
r = np.deg2rad(th)

fig, (axL, axR) = plt.subplots(1, 2, figsize=(7.4, 3.2),
                               gridspec_kw={"width_ratios": [1.35, 1.0]})

# ---------------- левая панель: разложение ----------------
axL.plot(th, np.cos(r) ** 4, color=ps.INK, lw=2.0, label=r"$\cos^4\theta$ (сумма)")
axL.plot(th, np.full_like(th, 3 / 8), color=ps.REF, ls=":", label=r"постоянная $3/8$")
axL.plot(th, 0.5 * np.cos(2 * r), color=ps.S1, label=r"2-я гармоника $\frac{1}{2}\cos 2\theta$")
axL.plot(th, 0.125 * np.cos(4 * r), color=ps.S2, label=r"4-я гармоника $\frac{1}{8}\cos 4\theta$")

axL.set_xticks([-90, -45, 0, 45, 90])
axL.set_xlim(-90, 90)
# Запас сверху: иначе легенда из четырёх строк ложится прямо на пик cos^4.
axL.set_ylim(-0.62, 1.62)
axL.axhline(0, color=ps.GRID, lw=0.8)
ps.finish(axL,
          title="Идеальная угловая кривая = 3 слагаемых",
          xlabel=r"угол оси пропускания $\theta$, град",
          ylabel="вклад в $U(\\theta)$",
          legend="upper right")

# ---------------- правая панель: A4/A2 ----------------
names = [n for n, _ in MEASURED]
vals = [v for _, v in MEASURED]
y = np.arange(len(vals))[::-1]

axR.axvline(IDEAL, color=ps.REF, ls="--", lw=1.4)
axR.text(IDEAL - 0.0006, -0.5, "идеал 0.25", color=ps.REF, fontsize=8,
         ha="right", va="bottom")

axR.scatter(vals, y, s=48, color=ps.S1, zorder=3)
# подписи значений — сбоку от точки, а не над ней: над точкой они наезжают
# на опорную линию «идеал 0.25» и на соседние метки
for v, yy in zip(vals, y):
    side = -1 if v < IDEAL else 1
    axR.annotate(f"{v:.4f}", (v, yy), textcoords="offset points",
                 xytext=(9 * side, 0), ha="right" if side < 0 else "left",
                 va="center", fontsize=8, color=ps.INK_2)

axR.set_yticks(y)
axR.set_yticklabels(names, fontsize=7)
axR.set_xlim(0.2385, 0.2615)
axR.set_ylim(-0.6, len(vals) - 0.15)
ps.finish(axR,
          title=r"Измеренное $A_4/A_2$",
          xlabel=r"$A_4/A_2$")
axR.grid(axis="y", visible=False)

fig.tight_layout(w_pad=2.0)
ps.save(fig, "fig00_harmonics")
