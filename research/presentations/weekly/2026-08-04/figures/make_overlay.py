#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Наложение модели M5 и измерений для `purewave` — в линейном и децибельном виде.

⚠ ОСОБЫЙ СТАТУС. По умолчанию зона P расчётов не ведёт (`STATE.md` п.1). Этот скрипт —
исключение по **прямому разрешению владельца от 2026-08-04**: «разрешаю писать и исполнять
любые скрипты для построения графиков на основе данных из зон других чатов».

Что скрипт делает и чего НЕ делает:
  · НЕ фитит. Параметры модели берутся ГОТОВЫМИ из артефакта зоны A
    `research/results/two_wgp/a7_purewave.json`, ветка `fit.day1_-95_80.best.params`
    (мультистарт по сетке из 21 старта, победитель — холодный старт).
  · Загружает измерения тем же билдером зоны A (`a3_eps_stability.build_masked`), теми же
    масками (водяная + динамический диапазон), что и сам фит, — чтобы картинка показывала
    ровно то, что видел оптимизатор.
  · Сводит спектр в одно число на угол: энергия Парсеваля `U(θ) = Σ_ν |T(θ,ν)|²` по
    частотам, валидным ОДНОВРЕМЕННО для всех углов набора (69 из 88 точек полосы
    0.21–1.49 ТГц) — иначе точки нельзя сравнивать между собой.
  · Нормирует обе кривые на `U(0°)` эксперимента, поэтому вертикальная шкала — «доля от
    максимума пропускания», а не абсолютные единицы.

Запуск из корня репозитория:
    .venv\\Scripts\\python.exe research/presentations/weekly/2026-08-04/figures/make_overlay.py
"""

import json
import sys
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[4]   # figures → 2026-08-04 → weekly → presentations → research → корень
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(REPO))
sys.path.insert(0, str(REPO / "research" / "two_wgp"))
sys.path.insert(0, str(REPO / "research" / "experiments"))

import lmfit                                    # noqa: E402
import a3_eps_stability as a3                   # noqa: E402  зона A, только чтение
from pstyle import SURFACE, INK, INK_2, S1, S2, S3, finish  # noqa: E402

ART = REPO / "research" / "results" / "two_wgp" / "a7_purewave.json"
RANGE = (-95.0, 80.0)          # день 1 — дни не сшиваются, мешать их нельзя (A7)
BRANCH = "day1_-95_80"


def load_params():
    """Готовые параметры лучшего бассейна M5_eps из артефакта A7 (без фита)."""
    d = json.loads(ART.read_text(encoding="utf-8"))
    best = d["fit"][BRANCH]["best"]
    P = lmfit.Parameters()
    for name, rec in best["params"].items():
        P.add(name, value=float(rec["v"]), vary=bool(rec["vary"]))
    return P, best


def compute():
    P, best = load_params()
    # День 1: диапазон, на котором оптимизировалась модель
    ang1, freqs, exp1, valid1, _meta = a3.build_masked("purewave", angles_limit=RANGE)
    common = valid1.all(axis=0)                      # частоты, годные на ВСЕХ углах дня 1

    U_exp1 = (np.abs(exp1) ** 2 * common).sum(axis=1)
    dense = np.linspace(-95.0, 115.0, 600)
    U_the = (np.abs(a3.model_grid(P, dense, freqs)) ** 2 * common).sum(axis=1)
    U_the_at_pts = (np.abs(a3.model_grid(P, ang1, freqs)) ** 2 * common).sum(axis=1)

    norm = U_exp1[np.argmin(np.abs(ang1))]           # нормировка на U(0°) эксперимента дня 1

    # День 2: все точки 82…110°, показываем на том же графике
    ang2, freqs2, exp2, valid2, _meta2 = a3.build_masked("purewave", angles_limit=(82.0, 110.0))
    common2 = valid2.all(axis=0)
    U_exp2 = (np.abs(exp2) ** 2 * common2).sum(axis=1)

    return {
        "angles1": ang1, "U_exp1": U_exp1 / norm,
        "angles2": ang2, "U_exp2": U_exp2 / norm,
        "dense": dense, "U_the": U_the / norm,
        "U_the_at_pts": U_the_at_pts / norm,
        "n_common": int(common.sum()), "n_freqs": int(len(freqs)),
        "best": best,
    }



def fig_linear(d):
    fig, ax = plt.subplots(figsize=(9.0, 4.8))
    ax.plot(d["dense"], d["U_the"], "-", color=S2, linewidth=2.2,
            label="модель M5 (параметры из фита дня 1)")
    ax.plot(d["angles1"], d["U_exp1"], "o", markersize=9, color=S1,
            markeredgecolor=SURFACE, markeredgewidth=1.5, linestyle="none",
            label="день 1 (− 95°…+80°)")
    ax.plot(d["angles2"], d["U_exp2"], "^", markersize=9, color=S3,
            markeredgecolor=SURFACE, markeredgewidth=1.5, linestyle="none",
            label="день 2 (+82°…+110°) — не сшиваются")
    ax.set_xlim(-100, 115)
    ax.set_ylim(-0.05, 1.14)
    finish(ax, "Пропускание vs. угол — линейная шкала",
           "угол поворота поляризатора θ, град.",
           "U(θ) / U(0°)")
    ax.legend(loc="upper left", frameon=False, fontsize=10)
    fig.tight_layout()
    fig.savefig(HERE / "fig4_purewave_overlay_linear.png", dpi=170)
    plt.close(fig)


def fig_db(d):
    fig, ax = plt.subplots(figsize=(9.0, 4.8))
    ax.plot(d["dense"], 10 * np.log10(d["U_the"]), "-", color=S2, linewidth=2.2,
            label="модель M5 (параметры из фита дня 1)")
    ax.plot(d["angles1"], 10 * np.log10(d["U_exp1"]), "o", markersize=9, color=S1,
            markeredgecolor=SURFACE, markeredgewidth=1.5, linestyle="none",
            label="день 1 (−95°…+80°)")
    ax.plot(d["angles2"], 10 * np.log10(d["U_exp2"]), "^", markersize=9, color=S3,
            markeredgecolor=SURFACE, markeredgewidth=1.5, linestyle="none",
            label="день 2 (+82°…+110°) — не сшиваются")
    ax.set_xlim(-100, 115)
    ax.set_ylim(-40, 4)
    finish(ax, "То же в децибелах — видно, где модель расходится с данными",
           "угол поворота поляризатора θ, град.",
           "10·lg[ U(θ)/U(0°) ], дБ")
    ax.legend(loc="lower center", frameon=False, fontsize=10)

    dev = 10 * np.log10(d["U_exp1"] / d["U_the_at_pts"])
    i = int(np.argmax(np.abs(dev)))
    ax.annotate(f"{dev[i]:+.1f} дБ", xy=(d["angles1"][i], 10 * np.log10(d["U_exp1"][i])),
                xytext=(d["angles1"][i] + 13, 10 * np.log10(d["U_exp1"][i]) - 7),
                color=INK, fontsize=10.5,
                arrowprops=dict(arrowstyle="->", color=INK_2, linewidth=1.2))
    fig.tight_layout()
    fig.savefig(HERE / "fig5_purewave_overlay_db.png", dpi=170)
    plt.close(fig)


if __name__ == "__main__":
    d = compute()
    fig_linear(d)
    fig_db(d)

    dev = 10 * np.log10(d["U_exp1"] / d["U_the_at_pts"])
    print(f"частот в сравнении: {d['n_common']} из {d['n_freqs']} в полосе")
    print(f"AIC фита (из артефакта): {d['best']['aic']:.1f}; "
          f"D_eff = {d['best']['D_um']:.3f} мкм; "
          f"DS_rmse = {d['best']['deepshadow_rmse_db']:.3f} дБ")
    print(f"точек день 1: {len(d['angles1'])}, день 2: {len(d['angles2'])}")
    print("готово: fig4, fig5")
