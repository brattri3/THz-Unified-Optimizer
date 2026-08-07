# -*- coding: utf-8 -*-
"""Рис. 9.1 -- профиль AIC(D) и компенсирующее движение eta0(D).

ВИЗУАЛИЗАЦИЯ УЖЕ ПОСЧИТАННЫХ ЧИСЕЛ (не новый расчёт): оба ряда читаются из
готового артефакта профильного правдоподобия
research/results/identifiability/identifiability.json, поле profile_Dscan
(на каждом шаге сетки по D остальные параметры, включая eta0, уже
переоптимизированы -- это тот же расчёт, что использован в §7 для другой
иллюстрации, здесь -- другая проекция тех же чисел: D против eta0).
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
eta0 = np.array([p["eta0"] for p in scan])
free_D = d["free_fit"]["D_um"]
free_D_err = d["free_fit"]["D_um_stderr"]

fig, (axA, axE) = plt.subplots(1, 2, figsize=(7.4, 3.2), sharex=True)

i_best = int(np.argmin(aic))
axA.plot(D, aic, "-o", color=ps.S1, ms=4, lw=1.4)
axA.plot([D[i_best]], [aic[i_best]], "*", color=ps.S3, ms=15, zorder=5,
         label=f"минимум профиля, $D={D[i_best]:.1f}$ мкм")
axA.axvspan(free_D - free_D_err, free_D + free_D_err, color=ps.S2, alpha=0.15,
            label=f"свободная оценка\n${free_D:.2f}\\pm{free_D_err:.2f}$ мкм")
ps.finish(axA, title=r"Профиль $\mathrm{AIC}(D)$", xlabel=r"$D$, мкм", ylabel="AIC",
          legend="lower right")

axE.plot(D, eta0, "-o", color=ps.S2, ms=4, lw=1.4)
ps.finish(axE, title=r"Компенсирующее $\eta_0(D)$", xlabel=r"$D$, мкм", ylabel=r"$\eta_0$")

fig.tight_layout(w_pad=2.6)
ps.save(fig, "fig09_profile")
