# -*- coding: utf-8 -*-
"""Рис. 4.2 -- согласие двух независимых оценок офсета (2-я и 4-я гармоника)
при относительном взвешивании sigma_i ~ U_i.

ВИЗУАЛИЗАЦИЯ УЖЕ ПОСЧИТАННЫХ ЧИСЕЛ: theta0 из 2-й и из 4-й гармоники для всех
четырёх образцов читаются из готового артефакта
research/results/two_wgp/a0_malus_harmonics_data.json (поле parseval_order4).
Подписи образцов -- физическое описание по D/P (Z_provenance.tex), без
внутренних имён наборов.
"""

import json
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[4]
sys.path.insert(0, str(REPO))

import numpy as np

import pstyle as ps
import matplotlib.pyplot as plt

ART = REPO / "research" / "results" / "two_wgp" / "a0_malus_harmonics_data.json"
d = json.loads(ART.read_text(encoding="utf-8"))["datasets"]

order = ["356att", "series1", "series2", "test_grid_40_20"]
labels = ["плотная,\n$D/P=0{,}71$", "повторная (а),\n$D/P=0{,}69$",
          "повторная (б),\n$D/P=0{,}69$", "разреженная,\n$D/P=0{,}46$"]

h2 = [d[k]["parseval_order4"]["theta0_deg_from_h2"] for k in order]
h4 = [d[k]["parseval_order4"]["theta0_deg_from_h4_mod45"] for k in order]
diff = [d[k]["parseval_order4"]["theta0_h4_minus_h2_deg"] for k in order]

x = np.arange(len(order))
fig, (axA, axB) = plt.subplots(1, 2, figsize=(7.4, 3.0))

w = 0.32
axA.bar(x - w / 2, h2, width=w, color=ps.S1, label=r"из 2-й гармоники")
axA.bar(x + w / 2, h4, width=w, color=ps.S2, label=r"из 4-й гармоники")
axA.axhline(0, color=ps.GRID, lw=0.8)
axA.set_xticks(x)
axA.set_xticklabels(labels, fontsize=7.5)
ps.finish(axA, title=r"Оценка офсета $\theta_0$", ylabel="градусы",
          legend="lower right")

axB.bar(x, np.abs(diff), color=ps.S3, width=0.5)
axB.set_xticks(x)
axB.set_xticklabels(labels, fontsize=7.5)
for xi, v in zip(x, diff):
    axB.text(xi, abs(v) + 0.002, f"{v:+.3f}°", ha="center", va="bottom", fontsize=7.5, color=ps.INK_2)
ps.finish(axB, title=r"Расхождение $|\theta_0^{(4)}-\theta_0^{(2)}|$", ylabel="градусы")

fig.tight_layout(w_pad=2.4)
ps.save(fig, "fig04_weights")
