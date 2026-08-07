# -*- coding: utf-8 -*-
"""Рис. 6.1 -- лестница dAIC для плотной решётки (D/P=0.71), один и два элемента тракта.

ВИЗУАЛИЗАЦИЯ УЖЕ ПОСЧИТАННЫХ ЧИСЕЛ (не новый расчёт): все значения AIC читаются
из готового артефакта лестницы моделей
research/results/two_wgp/a2_356att_ladder.json (поле dAIC_vs_M5_1wgp для каждого
варианта; M5_1wgp -- опорная модель курса, dAIC=0 по определению этого поля).
Подписи вариантов переведены на физическое описание ступени (без внутренних
кодовых имён), согласно STYLE.md; соответствие "подпись -> имя варианта в
артефакте" -- в комментариях ниже и в Z_provenance.tex, ключ 6.4.
"""

import json
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[4]
sys.path.insert(0, str(REPO))

import numpy as np

import pstyle as ps
import matplotlib.pyplot as plt

ART = REPO / "research" / "results" / "two_wgp" / "a2_356att_ladder.json"
variants = {v["label"]: v["dAIC_vs_M5_1wgp"] for v in json.loads(ART.read_text(encoding="utf-8"))["variants"]}

# один элемент тракта (идеальный анализатор)
single = [
    ("без просачивания,\n$D$ свободен", variants["M3ref_1wgp"]),
    ("опорная модель\n(+ просачивание)", variants["M5_1wgp"]),
    ("+ перекрёстный\nчлен", variants["M5_eps"]),
]

# два реальных элемента тракта
two = [
    ("выровн. ось,\nбез просачивания", variants["W2_shared_noleak"]),
    ("выровн. ось,\nс просачиванием", variants["W2_shared"]),
    ("разъюстированная\nось", variants["W2_mis"]),
    ("+ наклон\nизлучателя", variants["W2_mis_psi"]),
    ("+ перекрёстный\nчлен", variants["W2_mis_psi_eps_seedEps"]),
]

fig, (axL, axR) = plt.subplots(1, 2, figsize=(7.4, 3.6), sharex=True)

for ax, data, title, col in [
    (axL, single, "Один элемент тракта", ps.S1),
    (axR, two, "Два элемента тракта", ps.S2),
]:
    labels = [d[0] for d in data]
    vals = [d[1] for d in data]
    y = np.arange(len(data))[::-1]
    ax.hlines(y, 0, vals, color=col, lw=2.2)
    ax.plot(vals, y, "o", color=col, ms=6, zorder=5)
    ax.axvline(0, color=ps.INK_2, lw=1.0, ls="-")
    ax.set_yticks(y)
    ax.set_yticklabels(labels, fontsize=8)
    ax.set_xscale("symlog", linthresh=10)
    ps.finish(ax, title=title, xlabel=r"$\Delta$AIC относительно опорной модели")

axL.annotate("опорная\nмодель = 0", xy=(0, 1), xytext=(60, 1.35), fontsize=7.5,
             color=ps.INK_2, ha="left")

fig.tight_layout(w_pad=3.0)
ps.save(fig, "fig06_ladder")
