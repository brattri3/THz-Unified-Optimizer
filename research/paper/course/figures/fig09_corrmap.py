# -*- coding: utf-8 -*-
"""Рис. 9.2 -- корреляция свободной оценки D со всеми остальными параметрами.

ВИЗУАЛИЗАЦИЯ УЖЕ ПОСЧИТАННЫХ ЧИСЕЛ (не новый расчёт): значения читаются из
готового артефакта research/results/identifiability/identifiability.json,
поле free_fit.correl_with_D (ковариационная матрица подгонки lmfit в точке
свободного минимума). Подписи параметров -- физическое описание (без
внутренних имён полей), с расшифровкой в Z_provenance.tex, ключ 9.4.
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
correl = d["free_fit"]["correl_with_D"]

label_map = {
    "D_um~loss_factor": "коэфф. потерь",
    "D_um~gamma": "показатель рассеяния $\\gamma$",
    "D_um~angle_offset": "офсет угла",
    "D_um~tau_ps": "фазовая задержка",
    "D_um~tau_par_ps": "задержка канала\nподавления",
    "D_um~eta0": "амплитуда\nпросачивания $\\eta_0$",
    "D_um~eta_exp": "частотный ход\nпросачивания",
    "D_um~delta0": "фаза\nпросачивания $\\delta_0$",
    "D_um~tau_leak": "задержка\nпросачивания $\\tau_{\\rm leak}$",
}

keys = list(correl.keys())
labels = [label_map.get(k, k) for k in keys]
vals = [correl[k] for k in keys]

order = np.argsort(np.abs(vals))
labels = [labels[i] for i in order]
vals = [vals[i] for i in order]

fig, ax = plt.subplots(figsize=(6.2, 3.8))
colors = [ps.S2 if v < 0 else ps.S1 for v in vals]
y = np.arange(len(vals))
ax.barh(y, vals, color=colors)
ax.set_yticks(y)
ax.set_yticklabels(labels, fontsize=8.5)
ax.axvline(0, color=ps.INK_2, lw=0.8)
ax.set_xlim(-1.05, 1.05)
ps.finish(ax, title=r"Корреляция с $D$ (свободная подгонка)", xlabel="коэффициент корреляции")

fig.tight_layout()
ps.save(fig, "fig09_corrmap")
