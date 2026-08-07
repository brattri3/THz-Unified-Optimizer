# -*- coding: utf-8 -*-
"""Рис. 3.1 -- вклад членов потерь (Друде и степенное рассеяние).

МОДЕЛЬНЫЙ РАСЧЁТ (не данные): левая панель вызывает штатную реализацию модели
(unified_optimizer/model_blanco.py, compute_t_perp с use_drude=True/False) --
геометрия D/P=0.56 (образец средней плотности, см. Z_provenance.tex ключ 2.2),
частотная сетка -- рабочая полоса спектрометра. Правая панель -- чистая
математика: формула степенных потерь exp(-0.5*(L/4.343)*nu^gamma) для трёх
показателей gamma (два из них -- числа примера 3.1/3.2 параграфа, ключ 3.1/3.2),
без всякой подгонки под данные.

Данные не используются, подгонок нет.
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
DP = 0.56  # образец средней плотности (ключ 2.2)

freq_thz = np.linspace(0.2, 1.5, 300)
p_over_lambda = P_UM * freq_thz / 300.0  # lambda_um = 300 / freq_thz

fig, (axD, axG) = plt.subplots(1, 2, figsize=(7.4, 3.1))

t_perp_drude = np.array([
    mb.compute_t_perp(pl, DP, freq_thz=f, use_drude=True)
    for pl, f in zip(p_over_lambda, freq_thz)
])
t_perp_ideal = np.array([
    mb.compute_t_perp(pl, DP, use_drude=False) for pl in p_over_lambda
])

axD.plot(freq_thz, np.abs(t_perp_ideal), color=ps.REF, lw=1.4, ls="--",
         label="идеальный проводник (§2)")
axD.plot(freq_thz, np.abs(t_perp_drude), color=ps.S1, lw=1.8,
         label="с поправкой Друде")
ps.finish(axD,
          title=r"Канал пропускания $|t_\perp|$",
          xlabel="частота, ТГц", ylabel=r"$|t_\perp|$",
          legend="lower left")

# --- правая панель: чистая математика формулы (★) степенных потерь ---
L_DB = 3.0  # амплитудный масштаб потерь на верхней границе полосы -- иллюстративный
gammas = [(2.0, "рэлеевское, $\\gamma=2$ (разреженный образец, ключ 3.2)", ps.S3),
          (0.91, "$\\gamma=0{,}91$ (плотный, с просачиванием, ключ 3.1)", ps.S1),
          (4.24, "$\\gamma=4{,}24$ (плотный, без просачивания, ключ 3.1)", ps.S2)]

freq_g = np.linspace(0.2, 1.5, 300)
for gamma, lab, col in gammas:
    # калибруем L так, чтобы все три кривые совпадали на верхней границе полосы --
    # иллюстрация ФОРМЫ хода по частоте при одинаковых "потерях на верхнем крае",
    # а не подгонка под какие-либо измерения
    L = L_DB / (1.5 ** gamma)
    loss_amp_db = L * freq_g ** gamma
    axG.plot(freq_g, loss_amp_db, color=col, lw=1.8, label=lab)

ps.finish(axG,
          title=r"Степенные потери, $L\cdot\nu^\gamma$",
          xlabel="частота, ТГц", ylabel="потери, дБ (амплитуда)",
          legend="upper left")

fig.tight_layout(w_pad=2.2)
ps.save(fig, "fig03_drude")
