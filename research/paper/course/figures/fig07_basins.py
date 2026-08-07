# -*- coding: utf-8 -*-
"""Рис. 7.2 -- множественные локальные минимумы и доля стартов, достигающих
глобального.

ВИЗУАЛИЗАЦИЯ УЖЕ ПОСЧИТАННЫХ ЧИСЕЛ (не новый расчёт): левая панель -- реальное
сечение целевой функции (AIC) вдоль диаметра D при остальных параметрах,
оптимизированных на каждом шаге сетки (профильное правдоподобие), из готового
артефакта research/results/identifiability/identifiability.json, поле
profile_Dscan (тот же расчёт используется в §9 для другой иллюстрации --
идентифицируемости). Правая панель -- доля успешных стартов 22.2%, число из
research/SYNTHESIS.md / BRIEF_FOR_D_teaching_ladder.md §8.1 (уже посчитано в
проекте, здесь только отрисовано).
"""

import json
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[4]
sys.path.insert(0, str(REPO))

import numpy as np

import pstyle as ps
import matplotlib.pyplot as plt

ART = REPO / "research" / "results" / "identifiability" / "identifiability.json"
d = json.loads(ART.read_text(encoding="utf-8"))
scan = d["profile_Dscan"]
D = np.array([p["D_um"] for p in scan])
aic = np.array([p["aic"] for p in scan])
d_best = d["profile_min_D"]

fig, (axL, axR) = plt.subplots(1, 2, figsize=(7.4, 3.2), gridspec_kw={"width_ratios": [1.6, 1]})

axL.plot(D, aic, "-o", color=ps.S1, ms=4, lw=1.4)
i_best = int(np.argmin(aic))
axL.plot([D[i_best]], [aic[i_best]], "*", color=ps.S3, ms=15, zorder=5,
         label=f"глобальный минимум, $D={D[i_best]:.1f}$ мкм")
axL.axvline(4.694, color=ps.REF, lw=1.0, ls=":", label="локальный минимум,\nближайший к тёплому старту")
ps.finish(axL, title="Сечение целевой функции вдоль $D$",
          xlabel=r"$D$, мкм", ylabel="AIC", legend="lower right")

frac = 22.2
axR.bar([0, 1], [frac, 100 - frac], color=[ps.S1, ps.GRID], width=0.6)
axR.set_xticks([0, 1])
axR.set_xticklabels(["в глобальную\nобласть", "в другие\nобласти"], fontsize=8.5)
for xi, v in zip([0, 1], [frac, 100 - frac]):
    axR.text(xi, v + 1.5, f"{v:.1f}%", ha="center", fontsize=9, color=ps.INK)
axR.set_ylim(0, 100)
ps.finish(axR, title="Доля стартов сетки")

fig.tight_layout(w_pad=2.6)
ps.save(fig, "fig07_basins")
