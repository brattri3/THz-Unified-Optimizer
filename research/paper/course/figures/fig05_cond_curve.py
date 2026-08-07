# -*- coding: utf-8 -*-
"""Рис. 5.1 -- cond(X^T X) как функция углового диапазона.

ЧИСТАЯ ГЕОМЕТРИЯ МАТРИЦЫ ПЛАНА (не данные, не подгонка): для базиса
{1, cos2*theta, sin2*theta, cos4*theta, sin4*theta} (§4) число обусловленности
X^T X зависит только от того, НА КАКИХ УГЛАХ стоят строки матрицы. Это разрешённое
самостоятельное вычисление курса (PLAN.md §5: "иллюстрации чистой математики...
геометрия матрицы плана" -- явное исключение из правила "новых расчётов не
делаем"). Две выделенные точки (cond=2.14 при 19 угол -90..90 и cond=236.3 при
10 углов 0..90, шаг 10 град.) воспроизводят уже известные числа
(research/results/two_wgp/a0_malus_linearization.json) как проверку правильности
расчёта -- сама кривая между ними является НОВЫМ элементом, но того же самого
типа чистой геометрии, без обращения к каким-либо измерениям.
"""

import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[4]
sys.path.insert(0, str(REPO))

import numpy as np

import pstyle as ps
import matplotlib.pyplot as plt

STEP = 10.0  # шаг между углами, как в реальном плане измерения


def design_matrix(angles_deg):
    th = np.deg2rad(np.asarray(angles_deg))
    return np.column_stack([
        np.ones_like(th), np.cos(2 * th), np.sin(2 * th),
        np.cos(4 * th), np.sin(4 * th),
    ])


def cond_for_angles(angles_deg):
    X = design_matrix(angles_deg)
    XtX = X.T @ X
    ev = np.linalg.eigvalsh(XtX)
    return ev[-1] / ev[0]


half_widths = np.arange(20, 91, 2.0)

cond_sym = []
cond_onesided = []
for half in half_widths:
    n_sym = int(round(2 * half / STEP)) + 1
    angles_sym = np.linspace(-half, half, n_sym)
    cond_sym.append(cond_for_angles(angles_sym))

    n_one = int(round(half / STEP)) + 1
    angles_one = np.linspace(0, half, n_one)
    cond_onesided.append(cond_for_angles(angles_one))

fig, (axS, axO) = plt.subplots(1, 2, figsize=(7.4, 3.1), sharey=True)

axS.plot(half_widths, cond_sym, color=ps.S1, lw=1.8)
axS.set_yscale("log")
axS.axhline(2.14, color=ps.REF, lw=1.0, ls=":")
axS.plot([90], [2.1380711874576988], "o", color=ps.INK, ms=4, zorder=5)
axS.annotate("19 углов,\n$\\pm90^\\circ$: 2,14", xy=(90, 2.14), xytext=(48, 6),
             fontsize=8, color=ps.INK_2)
ps.finish(axS, title=r"Симметричный диапазон $\pm\Theta$",
          xlabel=r"$\Theta$, град.", ylabel=r"$\mathrm{cond}(X^TX)$")

axO.plot(half_widths, cond_onesided, color=ps.S2, lw=1.8)
axO.set_yscale("log")
axO.axhline(236.3, color=ps.REF, lw=1.0, ls=":")
axO.plot([90], [236.29477120206943], "o", color=ps.INK, ms=4, zorder=5)
axO.annotate("10 углов,\n$[0,90^\\circ]$: 236,3", xy=(90, 236.3), xytext=(38, 420),
             fontsize=8, color=ps.INK_2)
ps.finish(axO, title=r"Диапазон $[0,\Theta]$",
          xlabel=r"$\Theta$, град.")

fig.tight_layout(w_pad=2.2)
ps.save(fig, "fig05_cond_curve")
