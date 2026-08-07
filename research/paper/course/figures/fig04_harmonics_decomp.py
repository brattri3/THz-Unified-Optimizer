# -*- coding: utf-8 -*-
"""Рис. 4.1 -- разложение U(theta) на 2-ю и 4-ю гармоники: полное и урезанное.

ВИЗУАЛИЗАЦИЯ УЖЕ ПОСЧИТАННЫХ ЧИСЕЛ (не новый расчёт): коэффициенты c0,c2,c4
(через I_perp, I_par, Re_cross) и офсеты theta0 читаются из готового артефакта
гармонического разложения (research/results/two_wgp/a0_malus_harmonics_data.json,
образец "356att" -- в курсе "плотная решётка, D/P=0.71", см. Z_provenance.tex).
Ниже -- только алгебраическая пересборка кривой U(theta) из этих чисел по
формулам параграфа (c2=(I_perp-I_par)/2, c4=((I_perp+I_par)/2-Re_cross)/4,
c0=(I_perp+I_par)/2-c4) и её сравнение с версией без 4-й гармоники (order=2,
поле "parseval_order2_malus" того же артефакта). Никакой подгонки здесь не
происходит -- коэффициенты уже получены и записаны в артефакт.
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
d = json.loads(ART.read_text(encoding="utf-8"))
s = d["datasets"]["356att"]  # плотная решётка, D/P=0.71 (Z_provenance.tex)

p4 = s["parseval_order4"]
p2 = s["parseval_order2_malus"]
angles_meas = np.array(s["angles_deg"])

I_perp, I_par, Re_cross = p4["I_perp"], p4["I_par"], p4["Re_cross"]
c2 = (I_perp - I_par) / 2.0
c4 = ((I_perp + I_par) / 2.0 - Re_cross) / 4.0
c0 = (I_perp + I_par) / 2.0 - c4
th0_h2 = np.deg2rad(p4["theta0_deg_from_h2"])
th0_h4 = np.deg2rad(p4["theta0_deg_from_h4_mod45"])

theta = np.linspace(-90, 90, 721)
th = np.deg2rad(theta)

U4 = c0 + c2 * np.cos(2 * (th - th0_h2)) + c4 * np.cos(4 * (th - th0_h4))

th0_2 = np.deg2rad(p2["theta0_deg"])
U2 = p2["U_min"] + (p2["U_max"] - p2["U_min"]) * np.cos(th - th0_2) ** 2

# в dB относительно максимума полной (5-параметрической) кривой -- та же логика,
# что и в fig01_malus_floor: тень видна только в логарифмическом масштабе
U4_db = 10 * np.log10(np.clip(U4, 1e-9, None) / U4.max())
U2_db = 10 * np.log10(np.clip(U2, 1e-9, None) / U4.max())

fig, (axL, axD) = plt.subplots(1, 2, figsize=(7.4, 3.1))

axL.plot(theta, U4 / U4.max(), color=ps.S1, lw=1.8, label="полное разложение (5 коэфф.)")
axL.plot(theta, U2 / U4.max(), color=ps.S2, lw=1.4, ls="--", label="без 4-й гармоники")
axL.plot(angles_meas, np.interp(angles_meas, theta, U4 / U4.max()),
         "o", color=ps.INK, ms=3, zorder=5, label="измеренные углы")
ps.finish(axL, title="Линейный масштаб", xlabel=r"$\theta$, град.",
          ylabel=r"$U/U_{\max}$", legend="upper right")

axD.plot(theta, U4_db, color=ps.S1, lw=1.8, label="полное разложение")
axD.plot(theta, U2_db, color=ps.S2, lw=1.4, ls="--", label="без 4-й гармоники")
axD.set_xlim(60, 90)
axD.set_ylim(U4_db[np.searchsorted(theta, 60):].min() - 2, U4_db[np.searchsorted(theta, 60):].max() + 1)
ps.finish(axD, title="Децибелы, зона тени", xlabel=r"$\theta$, град.",
          ylabel=r"$10\lg(U/U_{\max})$, дБ", legend="lower left")

fig.tight_layout(w_pad=2.2)
ps.save(fig, "fig04_harmonics_decomp")
