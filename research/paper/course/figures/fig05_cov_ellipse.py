# -*- coding: utf-8 -*-
"""Рис. 5.2 -- эллипсы ковариации (a2,b2) для двух планов измерения.

ЧИСТАЯ ГЕОМЕТРИЯ МАТРИЦЫ ПЛАНА (не данные): Cov = (X^T X)^-1 при единичном и
одинаковом для всех точек шуме (sigma=1) -- условный масштаб, важна ФОРМА, а
не абсолютный размер эллипса. Два плана -- те же, что в fig05_cond_curve.py и
в таблице параграфа (19 углов -90..90 против 10 углов 0..90, шаг 10 град.).
"""

import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[4]
sys.path.insert(0, str(REPO))

import numpy as np

import pstyle as ps
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse


def design_matrix(angles_deg):
    th = np.deg2rad(np.asarray(angles_deg))
    return np.column_stack([
        np.ones_like(th), np.cos(2 * th), np.sin(2 * th),
        np.cos(4 * th), np.sin(4 * th),
    ])


def cov_a2b2(angles_deg):
    X = design_matrix(angles_deg)
    cov = np.linalg.inv(X.T @ X)
    idx = [1, 2]  # столбцы a2, b2
    return cov[np.ix_(idx, idx)]


angles_full = np.arange(-90, 91, 10.0)
angles_half = np.arange(0, 91, 10.0)

fig, ax = plt.subplots(figsize=(4.6, 4.3))

for angles, lab, col in [
    (angles_full, "полный диапазон, 19 углов", ps.S1),
    (angles_half, "урезанный, 10 углов", ps.S2),
]:
    cov = cov_a2b2(angles)
    vals, vecs = np.linalg.eigh(cov)
    order = np.argsort(vals)[::-1]
    vals, vecs = vals[order], vecs[:, order]
    angle_deg = np.degrees(np.arctan2(vecs[1, 0], vecs[0, 0]))
    width, height = 2 * np.sqrt(vals) * 2  # масштаб 2-сигма, условные единицы
    e = Ellipse((0, 0), width=width, height=height, angle=angle_deg,
                facecolor=col, edgecolor=col, alpha=0.18, lw=1.8)
    ax.add_patch(e)
    e2 = Ellipse((0, 0), width=width, height=height, angle=angle_deg,
                 facecolor="none", edgecolor=col, lw=1.8, label=lab)
    ax.add_patch(e2)

ax.set_xlim(-4, 4)
ax.set_ylim(-4, 4)
ax.set_aspect("equal")
ps.finish(ax, title=r"Эллипс ковариации $(a_2,b_2)$, усл. масштаб $2\sigma$",
          xlabel=r"$a_2$", ylabel=r"$b_2$", legend="upper right")

fig.tight_layout()
ps.save(fig, "fig05_cov_ellipse")
