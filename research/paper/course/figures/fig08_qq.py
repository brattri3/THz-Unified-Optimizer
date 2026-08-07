# -*- coding: utf-8 -*-
"""Рис. 8.1 -- график квантиль-квантиль остатка амплитуды, два образца.

ВИЗУАЛИЗАЦИЯ УЖЕ ПОСЧИТАННЫХ ЧИСЕЛ (не новый расчёт): остатки амплитуды
(измерение минус предсказание модели, дБ) читаются из готовых артефактов
research/results/residuals/{356att,test_grid_40_20}_grids.npz (поле
amp_res_db, сетка угол x частота). Здесь только строится Q-Q график по уже
вычисленным значениям остатка -- стандартная статистическая процедура над
готовыми числами, не подгонка и не новый физический расчёт.
"""

import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[4]
sys.path.insert(0, str(REPO))

import numpy as np
from scipy import stats

import pstyle as ps
import matplotlib.pyplot as plt

RES = REPO / "research" / "results" / "residuals"

fig, axes = plt.subplots(1, 2, figsize=(7.4, 3.3))

datasets = [
    ("356att", "плотная решётка, $D/P=0{,}71$", ps.S1, axes[0]),
    ("test_grid_40_20", "разреженный образец, $D/P=0{,}46$", ps.S3, axes[1]),
]

for name, label, col, ax in datasets:
    d = np.load(RES / f"{name}_grids.npz")
    res = d["amp_res_db"].ravel()
    res = res[np.isfinite(res)]
    res = (res - res.mean()) / res.std(ddof=1)

    (osm, osr), (slope, intercept, r) = stats.probplot(res, dist="norm")
    ax.plot(osm, osr, "o", color=col, ms=3.5, alpha=0.75)
    lims = [osm.min(), osm.max()]
    ax.plot(lims, lims, color=ps.REF, lw=1.2, ls="--", label="идеальная нормальность")
    ps.finish(ax, title=label, xlabel="теоретический квантиль",
              ylabel="квантиль остатка (норм.)",
              legend="upper left" if ax is axes[0] else False)

fig.tight_layout(w_pad=2.4)
ps.save(fig, "fig08_qq")
