# -*- coding: utf-8 -*-
"""Рис. 0.1 — соотношение масштабов: длина волны против периода решётки.

Что показывает: во всём рабочем диапазоне проекта (0.2–1.5 ТГц) длина волны на
порядок-два больше периода решётки (15–40 мкм). Отсюда подволновой режим: высшие
дифракционные порядки эванесцентны, решётка законно заменяется анизотропным листом.

Чистая кинематика (lambda = c / nu) — никаких данных и никаких подгонок.
Числа диапазонов: BRIEF_FOR_D_teaching_ladder.md §2 п.5.
"""

import numpy as np

import pstyle as ps
import matplotlib.pyplot as plt

C = 299_792_458.0  # м/с

nu = np.linspace(0.15, 1.6, 400)          # ТГц
lam = C / (nu * 1e12) * 1e6               # мкм

fig, ax = plt.subplots(figsize=(5.4, 3.2))

# Длина волны
ax.plot(nu, lam, color=ps.S1, label=r"длина волны $\lambda = c/\nu$")

# Полоса периодов решёток проекта
ax.axhspan(15, 40, color=ps.S3, alpha=0.18)
ax.text(1.53, 26, "период решёток\nпроекта $P = 15\\ldots40$ мкм",
        ha="right", va="center", fontsize=8, color=ps.S3)

# Граница, за которой появились бы распространяющиеся порядки: P = lambda
ax.plot(nu, lam, color=ps.S1)

# Рабочий диапазон по частоте
ax.axvspan(0.2, 1.5, color=ps.S1, alpha=0.06)
ax.text(0.85, 1750, "рабочий диапазон 0.2–1.5 ТГц",
        ha="center", fontsize=8, color=ps.INK_2)

ax.set_yscale("log")
ax.set_ylim(10, 2600)
ax.set_xlim(0.15, 1.6)

# Аннотация отношения масштабов на краях диапазона
for nu0, txt in ((0.2, r"$\lambda \approx 1500$ мкм"), (1.5, r"$\lambda \approx 200$ мкм")):
    lam0 = C / (nu0 * 1e12) * 1e6
    ax.plot([nu0], [lam0], "o", color=ps.S1, ms=5)
    ax.annotate(txt, (nu0, lam0), textcoords="offset points",
                xytext=(8, 8), fontsize=8, color=ps.S1)

ps.finish(ax,
          title="Волна на порядок-два длиннее периода решётки",
          xlabel=r"частота $\nu$, ТГц",
          ylabel=r"масштаб длины, мкм (лог)",
          legend=False)

ps.save(fig, "fig00_scales")
