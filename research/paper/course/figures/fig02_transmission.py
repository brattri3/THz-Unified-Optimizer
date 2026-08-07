# -*- coding: utf-8 -*-
"""Рис. 2.1 — модельные t_perp(nu), t_par(nu) по формулам Бланко (без Друде).

МОДЕЛЬНЫЙ РАСЧЁТ (не данные!): вызывает штатную реализацию модели
(unified_optimizer/model_blanco.py, compute_t_perp/compute_t_par, use_drude=False --
чистая электростатическая часть, без поправки на конечную проводимость, которая
вводится отдельным слоем позже в курсе). Геометрия -- P=25 мкм, две плотности D/P
из уже установленного в курсе диапазона образцов (0.42 разреженный .. 0.71 плотный,
см. figures/fig00_dp.py и Z_provenance.tex, ключ 2.2). Частотная сетка -- рабочая
полоса спектрометра, 0.2-3.0 ТГц.

Данные не используются, подгонок нет -- только формула, уже реализованная в проекте.
"""

import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[4]
sys.path.insert(0, str(REPO))

import numpy as np

import pstyle as ps
import matplotlib.pyplot as plt

from unified_optimizer import model_blanco as mb

P_UM = 25.0
DP = [(0.42, "разреженная, $D/P=0{,}42$", ps.S1),
      (0.71, "плотная, $D/P=0{,}71$", ps.S2)]

freq_thz = np.linspace(0.2, 3.0, 400)
p_over_lambda = P_UM * freq_thz / 300.0  # lambda_um = 300 / freq_thz

fig, (axT, axL) = plt.subplots(1, 2, figsize=(7.4, 3.1))

for dp, lab, col in DP:
    t_perp = np.array([mb.compute_t_perp(pl, dp, use_drude=False) for pl in p_over_lambda])
    t_par = np.array([mb.compute_t_par(pl, dp, use_drude=False) for pl in p_over_lambda])
    axT.plot(freq_thz, np.abs(t_perp), color=col, lw=1.8, label=lab)
    axL.plot(freq_thz, 20 * np.log10(np.abs(t_par)), color=col, lw=1.8, label=lab)

axT.set_ylim(0.94, 1.005)
ps.finish(axT,
          title=r"Канал пропускания $|t_\perp|$",
          xlabel="частота, ТГц", ylabel=r"$|t_\perp|$",
          legend="lower left")

ps.finish(axL,
          title=r"Канал подавления $|t_\parallel|$, дБ",
          xlabel="частота, ТГц", ylabel=r"$20\lg|t_\parallel|$, дБ",
          legend="upper right")

fig.tight_layout(w_pad=2.2)
ps.save(fig, "fig02_transmission")
