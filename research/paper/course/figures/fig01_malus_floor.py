# -*- coding: utf-8 -*-
"""Рис. 1.2 — пол просачивания: почему его видно только в логарифмическом масштабе.

МОДЕЛЬНЫЙ РАСЧЁТ (не данные!) по формуле U(theta) = |cos^2 t + eta e^{j d} sin^2 t|^2
с двумя реалистичными значениями просачивания — крайними точками измеренного диапазона
eta = 0.012...0.036 (плотная и разреженная решётки). Фаза просачивания взята нулевой:
она смещает кривую вблизи 45 градусов, но не меняет уровень пола, ради которого
строится рисунок.

Слева линейный масштаб: идеальная и реальная кривые неразличимы.
Справа децибелы: вместо нуля при 90 градусах — полка на уровне eta^2.

Данные не используются, подгонок нет — только формула из текста параграфа.
"""

import numpy as np

import pstyle as ps
import matplotlib.pyplot as plt

ETA = [(0.012, "плотная решётка, $\\eta = 0{,}012$", ps.S1),
       (0.036, "разреженная, $\\eta = 0{,}036$", ps.S2)]

th = np.linspace(0, 90, 1801)
r = np.deg2rad(th)
ideal = np.cos(r) ** 4


def u_leaky(eta, delta=0.0):
    amp = np.cos(r) ** 2 + eta * np.exp(1j * delta) * np.sin(r) ** 2
    return np.abs(amp) ** 2


fig, (axL, axR) = plt.subplots(1, 2, figsize=(7.4, 3.1))

# ---------------- линейный масштаб ----------------
axL.plot(th, ideal, color=ps.REF, lw=2.2, ls="--", label="идеальный: $\\cos^4\\theta$")
for eta, lab, col in ETA:
    axL.plot(th, u_leaky(eta), color=col, lw=1.6, label=lab)

axL.set_xticks([0, 30, 45, 60, 90])
axL.set_xlim(0, 90)
axL.set_ylim(-0.03, 1.05)
ps.finish(axL,
          title="Линейный масштаб: разницы не видно",
          xlabel=r"угол оси пропускания $\theta$, град",
          ylabel=r"$U(\theta)$, отн. ед.",
          legend="upper right")

# ---------------- децибелы ----------------
FLOOR_DB = -70.0
with np.errstate(divide="ignore"):
    ideal_db = np.maximum(10 * np.log10(ideal), FLOOR_DB)
axR.plot(th, ideal_db, color=ps.REF, lw=2.2, ls="--")
for eta, lab, col in ETA:
    axR.plot(th, 10 * np.log10(u_leaky(eta)), color=col, lw=1.6)
    db = 20 * np.log10(eta)
    axR.axhline(db, color=col, lw=0.9, ls=":", alpha=0.8)
    axR.annotate(f"пол $\\eta^2$: −{abs(db):.0f} дБ", (2, db),
                 textcoords="offset points", xytext=(0, 4), fontsize=8,
                 color=col, va="bottom")

axR.set_xticks([0, 30, 45, 60, 90])
axR.set_xlim(0, 90)
axR.set_ylim(-52, 4)
ps.finish(axR,
          title="Децибелы: пол просачивания очевиден",
          xlabel=r"угол оси пропускания $\theta$, град",
          ylabel=r"$10\lg U(\theta)$, дБ")
# Подпись к идеальной кривой — левее места, где она круто уходит вниз,
# иначе текст ложится прямо на неё.
axR.annotate("идеальная кривая\nуходит в $-\\infty$", (76, -47), fontsize=8,
             color=ps.REF, ha="right", va="bottom", linespacing=1.25)

fig.tight_layout(w_pad=2.2)
ps.save(fig, "fig01_malus_floor")
