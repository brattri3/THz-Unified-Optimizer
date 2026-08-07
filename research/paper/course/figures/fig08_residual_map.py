# -*- coding: utf-8 -*-
"""Рис. 8.2 -- карта остатка амплитуды (угол x частота), плотная решётка.

ВИЗУАЛИЗАЦИЯ УЖЕ ПОСЧИТАННЫХ ЧИСЕЛ (не новый расчёт): остаток (измерение
минус предсказание, дБ) на сетке угол x частота -- готовый артефакт
research/results/residuals/356att_grids.npz (поле amp_res_db). Это ИМЕННО
та карта, которая исторически впервые локализовала дефицит тени, приведший
к введению явного просачивания в модель (§1, §6) -- модель на этом срезе
ещё БЕЗ явного просачивания (config 0-90, тот же расчёт, что дал числа
ключа 8.1/8.2 в тексте).
"""

import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[4]
sys.path.insert(0, str(REPO))

import numpy as np

import pstyle as ps
import matplotlib.pyplot as plt

RES = REPO / "research" / "results" / "residuals"
d = np.load(RES / "356att_grids.npz")
angles = d["angles"]
freqs = d["freqs_THz"]
res = d["amp_res_db"]  # shape (n_angles, n_freqs)

fig, ax = plt.subplots(figsize=(6.6, 3.4))

vmax = np.nanpercentile(np.abs(res), 97)
im = ax.pcolormesh(freqs, angles, res, cmap="RdBu_r", vmin=-vmax, vmax=vmax,
                    shading="nearest")
cb = fig.colorbar(im, ax=ax, pad=0.02)
cb.set_label("остаток амплитуды, дБ")

ps.finish(ax, title="Остаток без явного просачивания (плотная решётка)",
          xlabel="частота, ТГц", ylabel=r"угол $\theta$, град.")

fig.tight_layout()
ps.save(fig, "fig08_residual_map")
