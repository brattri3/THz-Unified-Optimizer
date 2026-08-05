# -*- coding: utf-8 -*-
"""Рис. 0.4 — параметр плотности D/P: по нему выстраиваются все различия образцов.

Верхняя панель: где на шкале D/P лежат семь наборов проекта и где проходит граница
применимости аналитики Бланко, объявленная самим автором (d/p < 0.5).
Нижняя панель: дефицит эффективного диаметра D_eff/D_phys против D/P вместе с
эмпирическим законом 1 - 0.85*(D/P).

Все числа ПЕРЕНЕСЕНЫ из BRIEF_FOR_D_teaching_ladder.md §11.3 (колонка «полный угловой
диапазон»; значения с обрезки (0,90) смещены и здесь НЕ используются). Ничего не
пересчитывается: скрипт только рисует.
"""

import numpy as np

import pstyle as ps
import matplotlib.pyplot as plt

# name, D/P, D_eff/D_phys (полный угловой диапазон)
SAMPLES = [
    ("purewave (день 1)", 0.416, 0.659),
    ("test_grid_40_20",   0.464, 0.586),
    ("specac",            0.562, 0.443),
    ("att-11-16-s1",      0.688, 0.528),
    ("att-11-16-s2",      0.688, 0.452),
    ("att-11-16-s3",      0.688, 0.371),
    ("att-11-16-356",     0.710, 0.373),
]
BLANCO_LIMIT = 0.5   # объявленная автором граница приемлемости d/p

fig, (axT, axB) = plt.subplots(2, 1, figsize=(6.4, 4.6),
                               gridspec_kw={"height_ratios": [0.42, 1.0]},
                               sharex=True)

dp = np.array([s[1] for s in SAMPLES])
ratio = np.array([s[2] for s in SAMPLES])

# ------------------- верх: шкала плотности -------------------
axT.axvspan(0.38, BLANCO_LIMIT, color=ps.S3, alpha=0.12)
axT.axvspan(BLANCO_LIMIT, 0.75, color=ps.S2, alpha=0.10)
axT.axvline(BLANCO_LIMIT, color=ps.S2, ls="--", lw=1.4)
axT.text(BLANCO_LIMIT - 0.008, 0.62, "разреженный режим", ha="right", fontsize=8, color=ps.S3)
axT.text(BLANCO_LIMIT + 0.008, 0.62, "плотный режим — вне границы Бланко $d/p<0.5$",
         ha="left", fontsize=8, color=ps.S2)

axT.scatter(dp, np.full_like(dp, 0.25), s=42, color=ps.S1, zorder=3)

# Три набора att-11-16-s1/s2/s3 имеют ОДИН И ТОТ ЖЕ D/P = 0.688 (это один прибор
# в разных условиях съёмки), поэтому точки совпадают — подписываем их одной меткой,
# иначе три подписи ложатся друг на друга.
labels = {}
for name, x, _ in SAMPLES:
    labels.setdefault(round(x, 3), []).append(name)

for x, names_at_x in labels.items():
    if len(names_at_x) > 1:
        stem = names_at_x[0].rsplit("-", 1)[0]
        text = stem + "-" + "/".join(n.rsplit("-", 1)[1] for n in names_at_x)
    else:
        text = names_at_x[0]
    axT.annotate(text, (x, 0.25), textcoords="offset points", xytext=(0, -12),
                 rotation=32, ha="right", va="top", fontsize=7, color=ps.INK_2)

axT.set_ylim(-0.9, 1.0)
axT.set_yticks([])
axT.set_title("Все семь наборов на шкале плотности $D/P$",
              color=ps.TITLE, fontsize=10.5, loc="left", fontweight="bold", pad=8)
for side in ("top", "right", "left"):
    axT.spines[side].set_visible(False)

# ------------------- низ: дефицит D_eff -------------------
grid = np.linspace(0.38, 0.75, 100)
axB.plot(grid, 1 - 0.85 * grid, color=ps.REF, ls="--",
         label=r"эмпирический закон $1-0.85\,(D/P)$")
axB.axvline(BLANCO_LIMIT, color=ps.S2, ls="--", lw=1.2, alpha=0.7)

axB.scatter(dp, ratio, s=52, color=ps.S1, zorder=3,
            label=r"подгонка $D_{\rm eff}/D_{\rm phys}$")
axB.axhline(1.0, color=ps.S3, lw=1.6)
axB.text(0.395, 1.02, "точный расчёт (решёточные суммы): $\\approx 1$",
         fontsize=8, color=ps.S3, va="bottom")

axB.set_ylim(0.30, 1.12)
axB.set_xlim(0.38, 0.75)
ps.finish(axB,
          title="Дефицит эффективного диаметра растёт с плотностью",
          xlabel=r"$D/P$ — доля периода, занятая металлом",
          ylabel=r"$D_{\rm eff}/D_{\rm phys}$",
          legend="lower left")

fig.tight_layout(h_pad=0.6)
ps.save(fig, "fig00_dp")
