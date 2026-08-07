# -*- coding: utf-8 -*-
"""Рис. 7.1 -- иллюстрация алгоритма Левенберга-Марквардта на учебной задаче.

ЧИСТАЯ МАТЕМАТИКА / УЧЕБНАЯ ИЛЛЮСТРАЦИЯ (не данные и не параметры проекта):
цель -- показать МЕХАНИКУ алгоритма (переключение между градиентным спуском
и методом Гаусса-Ньютона, роль стартовой точки), а не какой-либо конкретный
параметр опорной модели курса. Игрушечная задача -- подгонка двухпараметрической
экспоненты y=a*exp(b*x) к синтетическим (без шума) точкам той же формы; это
классический пример из учебников по нелинейному МНК, использован здесь
исключительно как демонстрационный полигон для алгоритма.
"""

import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[4]
sys.path.insert(0, str(REPO))

import numpy as np

import pstyle as ps
import matplotlib.pyplot as plt

# --- игрушечная задача: y = a*exp(b*x), "истинные" a=2.0, b=0.30 ---
A_TRUE, B_TRUE = 2.0, 0.30
x_data = np.linspace(0, 6, 12)
y_data = A_TRUE * np.exp(B_TRUE * x_data)


def chi2(a, b):
    model = a * np.exp(b * x_data)
    return np.sum((y_data - model) ** 2)


def residuals_jac(a, b):
    model = a * np.exp(b * x_data)
    r = y_data - model
    J = np.column_stack([-np.exp(b * x_data), -a * x_data * np.exp(b * x_data)])
    return r, J


def levenberg_marquardt(a0, b0, n_iter=25, lam0=1e-2):
    p = np.array([a0, b0], dtype=float)
    lam = lam0
    path = [p.copy()]
    for _ in range(n_iter):
        r, J = residuals_jac(*p)
        g = J.T @ r
        H = J.T @ J
        step = np.linalg.solve(H + lam * np.diag(np.diag(H)) + 1e-12 * np.eye(2), -g)
        p_new = p + step
        if chi2(*p_new) < chi2(*p):
            p = p_new
            lam = max(lam / 3, 1e-6)
        else:
            lam = lam * 4
        path.append(p.copy())
        if np.linalg.norm(step) < 1e-8:
            break
    return np.array(path)

# сетка для линий уровня
aa = np.linspace(0.2, 4.5, 220)
bb = np.linspace(-0.05, 0.65, 220)
AA, BB = np.meshgrid(aa, bb)
Z = np.array([[chi2(a, b) for a in aa] for b in bb])
Z = np.log10(Z + 1e-6)

fig, (axG, axB) = plt.subplots(1, 2, figsize=(7.4, 3.4))

for ax, (a0, b0), title in [
    (axG, (0.6, 0.55), "удачный старт"),
    (axB, (4.2, 0.02), "неудачный старт (плато)"),
]:
    ax.contour(AA, BB, Z, levels=18, colors=ps.GRID, linewidths=0.6)
    path = levenberg_marquardt(a0, b0)
    ax.plot(path[:, 0], path[:, 1], "-o", color=ps.S1, ms=3, lw=1.4,
            label="траектория ЛМ")
    ax.plot([a0], [b0], "s", color=ps.S2, ms=6, zorder=5, label="старт")
    ax.plot([A_TRUE], [B_TRUE], "*", color=ps.S3, ms=13, zorder=5, label="минимум")
    ps.finish(ax, title=title, xlabel="параметр $a$", ylabel="параметр $b$",
              legend="upper right" if ax is axG else "lower left")

fig.tight_layout(w_pad=2.6)
ps.save(fig, "fig07_lm_trajectory")
