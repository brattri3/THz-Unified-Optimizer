# -*- coding: utf-8 -*-
"""Рис. 2.2 -- карта применимости модели по D/P и положение образцов.

Геометрия чистая математика: закрашенные зоны -- это две границы, объявленные
авторами модели (D/P<0.1 -- заявленная точность, D/P<0.5 -- приемлемость при
дополнительном условии D/lambda<0.2, которое у нас выполнено с большим запасом
на всей рабочей полосе, см. текст параграфа). Точки образцов -- уже установленные
в курсе значения D/P (см. fig00_dp.py, Z_provenance.tex ключ 2.2), новых
вычислений над данными нет.
"""

import numpy as np

import pstyle as ps
import matplotlib.pyplot as plt

SAMPLES = [
    (0.42, "наиболее\nразреженный"),
    (0.46, "разреженный"),
    (0.56, "средней\nплотности"),
    (0.69, "плотный\n(повтор)"),
    (0.71, "плотный"),
]

fig, ax = plt.subplots(figsize=(7.4, 2.5))

ax.axvspan(0.0, 0.1, color=ps.S3, alpha=0.16, lw=0)
ax.axvspan(0.1, 0.5, color=ps.S2, alpha=0.14, lw=0)
ax.axvspan(0.5, 0.85, color="#c94b4b", alpha=0.12, lw=0)

ax.axvline(0.1, color=ps.S3, lw=1.2, ls="--")
ax.axvline(0.5, color=ps.S2, lw=1.2, ls="--")

ax.text(0.05, 1.32, "точно\n$D/P<0{,}1$", ha="center", va="bottom", fontsize=8, color=ps.S3)
ax.text(0.30, 1.32, "приемлемо\n$0{,}1<D/P<0{,}5$", ha="center", va="bottom", fontsize=8, color=ps.S2)
ax.text(0.675, 1.32, "вне заявленной\nобласти", ha="center", va="bottom", fontsize=8, color="#c94b4b")

for dp, lab in SAMPLES:
    ax.plot([dp], [0.5], marker="o", ms=7, color=ps.INK, zorder=5)
    ax.annotate(f"{lab}\n$D/P={dp:.2f}$".replace(".", "{,}"), (dp, 0.5),
                textcoords="offset points", xytext=(0, -34), fontsize=7.5,
                ha="center", color=ps.INK_2, linespacing=1.2)

ax.set_xlim(0.0, 0.85)
ax.set_ylim(-0.35, 1.1)
ax.set_yticks([])
ax.set_xlabel(r"$D/P$")
for side in ("top", "right", "left"):
    ax.spines[side].set_visible(False)
ax.spines["bottom"].set_color(ps.GRID)

fig.tight_layout()
ps.save(fig, "fig02_validity_map")
