#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""A16 — единообразный пересчёт геометрии всех образцов и перенастройка закона H5.

ЗАДАЧА ВЛАДЕЛЬЦА: «применить метод оценки диаметра по фото и пересчитать все
отношения правильно, откорректировать H5».

ЧТО ВЫЯСНИЛОСЬ ПО ХОДУ — три вещи, каждая меняет метод
------------------------------------------------------
1. МАСШТАБ ПЛАВАЕТ МЕЖДУ КАДРАМИ. Отношение периодов в пикселях high/mid должно
   быть равно отношению калибровочных масштабов 5.098, а фактически: 5.02
   (test_grid_33_11), 5.29 (test_grid_40_20), 5.39 (purewave), **4.39** (specac).
   Общий масштаб `high` применять нельзя.
   РЕШЕНИЕ: каждый кадр большого увеличения калибруется по СОБСТВЕННОМУ периоду
   образца (`мкм/px = P_мкм / P_px`), период берётся со среднего увеличения, где
   гребёнка разрешена и периодов много. Зум-дрейф сокращается тождественно.

2. РЕГРЕССИЯ «ШИРИНА ~ РАЗМЫТИЕ» ПЕРЕОЦЕНИВАЕТ ПОПРАВКУ. Сильная корреляция
   (r = 0.7…0.9) создаётся ХВОСТОМ сильно расфокусированных проволок. Если
   ограничиться слабо размытыми (b < 2.5 мкм), наклон разваливается:
   +0.47 / +1.88 / −0.34 / +0.16 у четырёх образцов, то есть в рабочей области
   ширина от размытия почти не зависит. Значит линейная модель годится для
   ДИАГНОСТИКИ (эффект есть), но не для экстраполяции.
   ⚠ СЛЕДСТВИЕ: число `D = 10.91` из A14b ПЕРЕОЦЕНЕНО поправкой и отзывается.
   РЕШЕНИЕ: оценивать `D` как медиану ширины по САМЫМ РЕЗКИМ проволокам
   (нижняя четверть по размытию). Эта оценка устойчива к выбору порога
   (проверено на 1.6/1.8/2.0/2.2/2.5 мкм), в отличие от экстраполяции.

3. ЧЕТЫРЕ АТТЕНЮАТОРНЫЕ ТОЧКИ ИЗМЕРИТЬ НЕЧЕМ. У `att-11-16-*` есть только кадры
   среднего увеличения, проволоки почти сомкнуты (`D/P ~ 0.69`), FFT цепляется за
   кластеры яркости (62/103/146 px вместо ожидаемых ~25). Они остаются с
   ПАСПОРТНОЙ геометрией и в перенастройку закона не входят.

ЧЕСТНАЯ ОГОВОРКА, встроенная в вывод: оценка по самым резким проволокам меряет
`D_истинный + аппаратное уширение` (~1 мкм, предел резкости прибора по гребёнке).
Это уширение общее для всех образцов, поэтому в СРАВНЕНИИ между образцами оно
почти сокращается, но абсолютные `D` систематически завышены на ~1 мкм.

Запуск ИЗ КОРНЯ репозитория:
    .venv/Scripts/python.exe research/two_wgp/a16_h5_recompute.py
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt                                      # noqa: E402

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "research" / "two_wgp"))
sys.path.insert(0, str(REPO / "research" / "experiments"))
import a16_geometry_all as a16                                       # noqa: E402

OUT = REPO / "research" / "results" / "two_wgp"
SHARP_FRACTION = 0.25          # доля самых резких проволок для оценки D
BLUR_FLOOR_UM = a16.BLUR_FLOOR_UM

# D_eff из мультистартовых фитов (полный угловой диапазон).
# Источники: a3_eps_stability.json, a7_purewave.json, a14_refit_measured_geometry.json.
D_EFF = {
    "purewave": 6.733,             # день 1, ПЕРЕФИТ при P = 24.72 (было 6.982 при P = 25.5)
    "specac": 6.202,
    "test_grid_40_20": 10.548,
    "test_grid_33_11": 7.233,      # перефит на P = 31.75
}
# Аттенюаторные точки — паспортная геометрия, диаметр не измерялся.
ATT = [("att-11-16-356", 15.5, 11.0, 4.107),
       ("att-11-16-s3", 16.0, 11.0, 4.079),
       ("att-11-16-s2", 16.0, 11.0, 4.967),
       ("att-11-16-s1", 16.0, 11.0, 5.811)]
PUBLISHED = {"purewave": (25.5, 10.6), "specac": (24.9, 14.0),
             "test_grid_40_20": (38.8, 18.0), "test_grid_33_11": (31.75, 10.91)}


def sharpest_D(name):
    """Медиана ширины по самым резким проволокам + разброс по порогам."""
    cfg = a16.SAMPLES[name]
    dd = a16.MIC / cfg["dir"]
    # период со среднего увеличения
    Ps = []
    for f in cfg["mid"]:
        prof = a16.profile(dd / f)
        per = a16.period_px(prof)
        peaks, _ = __import__("micrograph_period").peak_spacings(prof, per)
        P_lat, _r, _n, _m = __import__("micrograph_period").lattice_fit(peaks, per)
        Ps.append(P_lat * a16.UM_PER_PX_MID)
    P_um = float(np.mean(Ps))
    P_err = float(np.std(Ps, ddof=1)) if len(Ps) > 1 else 0.02 * P_um

    rows = []
    for f in cfg["high"]:
        prof = a16.profile(dd / f)
        per_px = a16.period_px(prof)
        r, _ = a16.wires_of(dd / f, P_um / per_px)
        rows += r
    b = np.array([x["blur_um"] for x in rows])
    w = np.array([x["w10_um"] for x in rows])

    k = max(3, int(SHARP_FRACTION * len(w)))
    idx = np.argsort(b)[:k]
    D = float(np.median(w[idx]))
    # разброс оценки по альтернативным порогам — честная систематика метода
    alts = []
    for cut in (1.6, 1.8, 2.0, 2.2, 2.5):
        m = b < cut
        if m.sum() >= 3:
            alts.append(float(np.median(w[m])))
    D_err = float(max(np.std(alts, ddof=1) if len(alts) > 1 else 0.3,
                      np.std(w[idx], ddof=1) / np.sqrt(k)))
    return {"P_um": P_um, "P_err_um": P_err, "D_um": D, "D_err_um": D_err,
            "n_wires": len(w), "n_sharp": k,
            "blur_sharp_um": [float(b[idx].min()), float(b[idx].max())],
            "D_alternatives": alts, "median_all_wires": float(np.median(w))}


def fit_h5(dp, ratio, through_one=True):
    """Закон H5: ratio = 1 - k*(D/P) (через единицу) или свободная прямая."""
    dp = np.asarray(dp, float)
    y = np.asarray(ratio, float)
    if through_one:
        k = float(np.sum(dp * (1 - y)) / np.sum(dp * dp))
        resid = y - (1 - k * dp)
        dof = max(len(y) - 1, 1)
        k_err = float(np.sqrt(np.sum(resid ** 2) / dof / np.sum(dp * dp)))
        return {"form": "1 - k*(D/P)", "k": k, "k_err": k_err,
                "rms": float(np.sqrt(np.mean(resid ** 2)))}
    A = np.vstack([dp, np.ones_like(dp)]).T
    c, *_ = np.linalg.lstsq(A, y, rcond=None)
    resid = y - A @ c
    return {"form": "a + b*(D/P)", "b": float(c[0]), "a": float(c[1]),
            "rms": float(np.sqrt(np.mean(resid ** 2)))}


def main():
    print(f"аппаратное уширение (предел резкости по гребёнке): {BLUR_FLOOR_UM:.2f} мкм\n")
    geo, rows = {}, []
    print(f"{'образец':<18}{'P нов.':>9}{'P стар.':>9}{'D нов.':>9}{'D стар.':>9}"
          f"{'D/P нов.':>10}{'D/P стар.':>11}")
    for name in a16.SAMPLES:
        g = sharpest_D(name)
        geo[name] = g
        Pp, Dp = PUBLISHED[name]
        dp_new = g["D_um"] / g["P_um"]
        ratio = D_EFF[name] / g["D_um"]
        rows.append({"sample": name, "D_over_P": dp_new,
                     "D_eff_over_D_phys": ratio, "measured": True})
        print(f"{name:<18}{g['P_um']:>9.2f}{Pp:>9.2f}{g['D_um']:>9.2f}{Dp:>9.2f}"
              f"{dp_new:>10.3f}{Dp / Pp:>11.3f}")
    print()
    for nm, P, D, deff in ATT:
        rows.append({"sample": nm, "D_over_P": D / P,
                     "D_eff_over_D_phys": deff / D, "measured": False})

    meas = [r for r in rows if r["measured"]]
    dp_m = [r["D_over_P"] for r in meas]
    print(f"⚠ диапазон D/P у измеренных образцов сжался: "
          f"{min(dp_m):.3f}…{max(dp_m):.3f} (было "
          f"{min(PUBLISHED[n][1] / PUBLISHED[n][0] for n in PUBLISHED):.3f}…"
          f"{max(PUBLISHED[n][1] / PUBLISHED[n][0] for n in PUBLISHED):.3f})")

    fits = {
        "measured_only": fit_h5([r["D_over_P"] for r in meas],
                                [r["D_eff_over_D_phys"] for r in meas]),
        "all_incl_passport": fit_h5([r["D_over_P"] for r in rows],
                                    [r["D_eff_over_D_phys"] for r in rows]),
        "old_law_k": 0.85,
    }
    print(f"\n--- перенастройка закона H5 (форма 1 - k*(D/P)) ---")
    print(f"  только измеренные (4 точки):  k = {fits['measured_only']['k']:.3f} "
          f"± {fits['measured_only']['k_err']:.3f}   rms {fits['measured_only']['rms']:.3f}")
    print(f"  все 8 (4 с паспортной геом.): k = {fits['all_incl_passport']['k']:.3f} "
          f"± {fits['all_incl_passport']['k_err']:.3f}   rms {fits['all_incl_passport']['rms']:.3f}")
    print(f"  прежний закон:                k = 0.850")

    # ---------- рисунок ----------
    fig, ax = plt.subplots(figsize=(8.2, 5.6))
    xs = np.linspace(0.15, 0.75, 100)
    ax.plot(xs, 1 - 0.85 * xs, "--", lw=1.5, color="#909090",
            label="прежний закон  1 − 0.85·D/P")
    k = fits["measured_only"]["k"]
    ax.plot(xs, 1 - k * xs, "-", lw=1.8, color="#c0504d",
            label=f"перенастроенный  1 − {k:.2f}·D/P (4 измеренные точки)")
    for r in rows:
        col, mk = ("#c0504d", "s") if r["measured"] else ("#909090", "o")
        ax.plot(r["D_over_P"], r["D_eff_over_D_phys"], mk, ms=9 if r["measured"] else 7,
                color=col, mfc=col if r["measured"] else "none", mew=1.5)
        ax.annotate(r["sample"].replace("att-11-16-", "att-").replace("test_grid_", "grid "),
                    (r["D_over_P"], r["D_eff_over_D_phys"]), fontsize=7.5,
                    xytext=(4, 4), textcoords="offset points", color=col)
    for n in PUBLISHED:
        P, D = PUBLISHED[n]
        ax.plot(D / P, D_EFF[n] / D, "x", ms=8, color="#404040", mew=1.4)
    ax.plot([], [], "x", ms=8, color="#404040", mew=1.4, label="было (прежняя геометрия)")
    ax.plot([], [], "o", ms=7, mfc="none", color="#909090", mew=1.5,
            label="att-11-16: паспорт, диаметр НЕ измерен")
    ax.set_xlabel("D/P (по микрофото; att — паспорт)")
    ax.set_ylabel(r"$D_{eff}/D_{phys}$")
    ax.set_title("Пересчёт геометрии по единой методике и перенастройка закона H5\n"
                 "⚠ диапазон D/P сжался, наклон определён слабо", fontsize=11)
    ax.grid(alpha=0.3)
    ax.legend(fontsize=8, loc="upper right")
    fig.tight_layout()
    png = OUT / "a16_h5_recompute.png"
    fig.savefig(png, dpi=155)

    payload = {"method": "медиана ширины на 10 % по самым резким 25 % проволок; "
                         "каждый кадр самокалиброван по периоду образца",
               "instrument_broadening_um": BLUR_FLOOR_UM,
               "geometry": geo, "points": rows, "h5_fits": fits,
               "published_geometry": PUBLISHED, "D_eff_used": D_EFF,
               "attenuator_note": a16.NOT_MEASURABLE}
    js = OUT / "a16_h5_recompute.json"
    js.write_text(json.dumps(payload, ensure_ascii=False, indent=2, default=float),
                  encoding="utf-8")
    print(f"\nартефакты: {js}\n           {png}")


if __name__ == "__main__":
    main()
