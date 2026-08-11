# -*- coding: utf-8 -*-
"""A19 — бюджет погрешностей точек (D/P, D_eff/D_phys) и взвешенная подгонка H5.

ИСТОЧНИКИ ПОГРЕШНОСТЕЙ (всё из уже посчитанных артефактов, ничего не выдумано):

  диаметр D
    стат.  std/sqrt(n) по проволокам образца            a18_<набор>.json: pooled.D
    сист.  конвенция границы: размах скана порога       a18_<набор>.json: threshold_effect.dD_um
           0.2..0.8 контраста как равномерное -> /sqrt(12)
  период P
    стат.  стандартная ошибка среднего по кадрам        a18_<набор>.json: P_density.sem
    сист.  расхождение двух способов счёта              |P_density - P_best| (§10.1 отчёта)
  эффективный диаметр D_eff
    фит    стандартная ошибка параметра                 a19_refit.json: params.D_um.err
    от P   период ПРИКОЛОТ, его неопределённость течёт
           в D_eff с чувствительностью S = dlnD_eff/dlnP,
           измеренной по паре прогонов A16 -> A19        a19_refit.json + a16_h5_recompute.json

КЛЮЧЕВОЕ ПРО КОНВЕНЦИЮ ГРАНИЦЫ. Она ОДНА И ТА ЖЕ для всех образцов, то есть двигает
все точки согласованно, а не независимо. Класть её в погрешность каждой точки при
проверке разброса вокруг прямой неправильно — она даёт сдвиг наклона k, а не рассеяние.
Поэтому она считается отдельно: подгонка повторяется при конвенции владельца, и разница
в k объявляется систематикой на k.

ОШИБКИ ПО ОСЯМ КОРРЕЛИРОВАНЫ: x = D/P и y = D_eff/D содержат один и тот же D в
противоположных местах, а P входит и в x, и (через пин) в D_eff. Поэтому строятся
эллипсы, а не кресты, и подгонка использует эффективную дисперсию остатка.

Запуск из корня репозитория:
  .venv\\Scripts\\python.exe research\\results\\a18_microscopy\\figs_report\\make_h5_errors.py
"""
import json
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patheffects as pe
from matplotlib.patches import Ellipse

HERE = Path(__file__).resolve().parent
RES = HERE.parent
REPO = HERE.parents[3]
A16 = REPO / "research" / "results" / "two_wgp" / "a16_h5_recompute.json"

SURFACE, INK, INK2, MUTED, GRID = "#fcfcfb", "#0b0b0b", "#52514e", "#8a8a85", "#e3e3df"
S1, S2, S3, RED = "#2a78d6", "#eb6834", "#1baf7a", "#c1272d"

EXCLUDED = {"purewave"}          # провисание проволоки в зоне пучка (владелец, 2026-08-11)
OWNER_EDGE_SHIFT_UM = 0.6
ORDER = ["test_grid_33_11", "specac", "purewave", "test_grid_40_20"]
NAMES = {"purewave": "purewave", "specac": "specac",
         "test_grid_33_11": "test_grid 33/11", "test_grid_40_20": "test_grid 40/20"}


def budget():
    a16 = json.loads(A16.read_text(encoding="utf-8"))
    refit = {r["sample"]: r for r in
             json.loads((RES / "a19_refit.json").read_text(encoding="utf-8"))["rows"]}
    runs = json.loads((RES / "a19_refit.json").read_text(encoding="utf-8"))["runs"]

    rows = []
    for ds in ORDER:
        a18 = json.loads((RES / f"a18_{ds}.json").read_text(encoding="utf-8"))
        D = a18["pooled"]["D"]["median"]
        sD_stat = a18["pooled"]["D"]["std"] / np.sqrt(a18["pooled"]["D"]["n"])
        sD_sys = a18["threshold_effect"]["dD_um"] / np.sqrt(12.0)

        P = a18["P_density"]["mean"]
        sP_stat = a18["P_density"]["sem"]
        sP_sys = abs(a18["P_density"]["mean"] - a18["P_best_um"])
        sP = float(np.hypot(sP_stat, sP_sys))

        De = refit[ds]["D_eff_um"]
        sDe_fit = runs[ds]["best"]["params"]["D_um"]["err"]
        # чувствительность приколотого периода: измерена по паре прогонов A16 -> A19
        P_old = a16["geometry"][ds]["P_um"]
        De_old = a16["D_eff_used"][ds]
        S = float(np.log(De / De_old) / np.log(P / P_old))
        sDe_fromP = abs(S) * De * sP / P
        sDe = float(np.hypot(sDe_fit, sDe_fromP))

        x, y = D / P, De / D
        # производные; сист. конвенции сюда НЕ входит — она общая, считается отдельно
        dx_dD, dx_dP = 1.0 / P, -D / P ** 2
        dy_dD, dy_dP = -De / D ** 2, (De / D) * S / P
        sx = float(np.sqrt((dx_dD * sD_stat) ** 2 + (dx_dP * sP) ** 2))
        sy = float(np.sqrt((dy_dD * sD_stat) ** 2 + (sDe_fit / D) ** 2
                           + (dy_dP * sP) ** 2))
        cov = float(dx_dD * dy_dD * sD_stat ** 2 + dx_dP * dy_dP * sP ** 2)
        rho = cov / (sx * sy)

        rows.append(dict(sample=ds, D=D, sD_stat=sD_stat, sD_sys=sD_sys, n_wires=a18["pooled"]["D"]["n"],
                         P=P, sP_stat=sP_stat, sP_sys=sP_sys, sP=sP,
                         De=De, sDe_fit=sDe_fit, S=S, sDe_fromP=sDe_fromP, sDe=sDe,
                         x=x, y=y, sx=sx, sy=sy, cov=cov, rho=rho,
                         excluded=ds in EXCLUDED))
    return rows


def wfit(rows, k0=0.93, iters=8):
    """Взвешенный МНК y = 1 - k*x с эффективной дисперсией остатка.

    Var(r) = sy^2 + k^2 sx^2 + 2k cov  — учитывает и погрешность по x, и корреляцию.
    """
    x = np.array([r["x"] for r in rows])
    y = np.array([r["y"] for r in rows])
    sx = np.array([r["sx"] for r in rows])
    sy = np.array([r["sy"] for r in rows])
    cov = np.array([r["cov"] for r in rows])
    k = k0
    for _ in range(iters):
        var = sy ** 2 + k ** 2 * sx ** 2 + 2 * k * cov
        w = 1.0 / var
        k = float(np.sum(w * x * (1 - y)) / np.sum(w * x * x))
    var = sy ** 2 + k ** 2 * sx ** 2 + 2 * k * cov
    w = 1.0 / var
    r = y - (1 - k * x)
    chi2 = float(np.sum(w * r ** 2))
    dof = len(x) - 1
    return dict(k=k, k_err=float(1 / np.sqrt(np.sum(w * x * x))),
                chi2=chi2, dof=dof, chi2_red=chi2 / dof if dof else float("nan"),
                resid=r.tolist(), sigma_eff=np.sqrt(var).tolist())


def owner_edge_fit(rows):
    """Тот же расчёт при конвенции края владельца: D -> D + 0.6 мкм (общий сдвиг)."""
    shifted = []
    for r in rows:
        D = r["D"] + OWNER_EDGE_SHIFT_UM
        s = dict(r)
        s.update(x=D / r["P"], y=r["De"] / D)
        shifted.append(s)
    return wfit(shifted)


def main():
    rows = budget()
    keep = [r for r in rows if not r["excluded"]]

    fit = wfit(keep)
    fit_all = wfit(rows)
    fit_edge = owner_edge_fit(keep)
    k_sys_edge = abs(fit_edge["k"] - fit["k"])

    print("БЮДЖЕТ ПОГРЕШНОСТЕЙ (1 сигма)\n")
    print(f"{'образец':<17}{'D':>7}{'ст.D':>7}{'сис.D':>7}{'P':>8}{'ст.P':>7}{'сис.P':>7}"
          f"{'D_eff':>8}{'фит':>7}{'S':>6}{'от P':>7}")
    for r in rows:
        print(f"{r['sample']:<17}{r['D']:>7.2f}{r['sD_stat']:>7.3f}{r['sD_sys']:>7.3f}"
              f"{r['P']:>8.2f}{r['sP_stat']:>7.3f}{r['sP_sys']:>7.3f}"
              f"{r['De']:>8.3f}{r['sDe_fit']:>7.3f}{r['S']:>6.2f}{r['sDe_fromP']:>7.3f}")

    print(f"\n{'образец':<17}{'x=D/P':>9}{'sx':>8}{'y':>9}{'sy':>8}{'rho':>8}{'ост.':>9}{'ост/sig':>9}")
    for r, res, se in zip(keep, fit["resid"], fit["sigma_eff"]):
        print(f"{r['sample']:<17}{r['x']:>9.4f}{r['sx']:>8.4f}{r['y']:>9.4f}"
              f"{r['sy']:>8.4f}{r['rho']:>8.2f}{res:>+9.4f}{res / se:>+9.2f}")
    for r in rows:
        if r["excluded"]:
            pred = 1 - fit["k"] * r["x"]
            se = np.sqrt(r["sy"] ** 2 + fit["k"] ** 2 * r["sx"] ** 2 + 2 * fit["k"] * r["cov"])
            print(f"{r['sample']:<17}{r['x']:>9.4f}{r['sx']:>8.4f}{r['y']:>9.4f}"
                  f"{r['sy']:>8.4f}{r['rho']:>8.2f}{r['y'] - pred:>+9.4f}"
                  f"{(r['y'] - pred) / se:>+9.2f}  ИСКЛЮЧЁН")

    print(f"\nвзвешенная подгонка y = 1 - k*(D/P), без purewave:")
    print(f"  k = {fit['k']:.4f} ± {fit['k_err']:.4f} (стат) ± {k_sys_edge:.4f} (конвенция края)")
    print(f"  chi2 = {fit['chi2']:.3f} при dof = {fit['dof']}  ->  chi2/dof = {fit['chi2_red']:.3f}")
    print(f"  та же подгонка при конвенции владельца: k = {fit_edge['k']:.4f}")
    print(f"  с purewave: k = {fit_all['k']:.4f} ± {fit_all['k_err']:.4f}, "
          f"chi2/dof = {fit_all['chi2_red']:.2f}")

    # ------------------------------- рисунок -------------------------------
    fig, ax = plt.subplots(figsize=(7.2, 4.8), facecolor=SURFACE)
    ax.set_facecolor(SURFACE)
    grid = np.linspace(0.28, 0.76, 200)
    k = fit["k"]
    ax.plot(grid, 1 - k * grid, color=S1, ls="--", lw=1.8,
            label=f"закон H5: $1-{k:.3f}\\,(D/P)$")
    band = fit["k_err"] + k_sys_edge
    ax.fill_between(grid, 1 - (k - band) * grid, 1 - (k + band) * grid,
                    color=S1, alpha=0.10, lw=0,
                    label="коридор $k$: стат. + конвенция края")
    ax.axvline(0.5, color=S2, ls="--", lw=1.2, alpha=0.75)
    ax.text(0.505, 0.775, "граница применимости\nмодели Бланко $D/P<0.5$",
            fontsize=8, color=S2, va="top")

    for r in rows:
        c = RED if r["excluded"] else S1
        C = np.array([[r["sx"] ** 2, r["cov"]], [r["cov"], r["sy"] ** 2]])
        val, vec = np.linalg.eigh(C)
        ang = np.degrees(np.arctan2(vec[1, np.argmax(val)], vec[0, np.argmax(val)]))
        w, h = 2 * np.sqrt(np.maximum(val[::-1], 0))
        ax.add_patch(Ellipse((r["x"], r["y"]), w, h, angle=ang, facecolor=c,
                             alpha=0.22, edgecolor=c, lw=1.2, zorder=4))
        ax.scatter([r["x"]], [r["y"]], s=90 if r["excluded"] else 46,
                   marker="X" if r["excluded"] else "o", color=c, zorder=6)

    ax.scatter([], [], s=46, color=S1, label="эллипс $1\\sigma$: стат. D, P и $D_{\\rm eff}$")
    ax.scatter([], [], s=90, marker="X", color=RED,
               label="исключён: провисание проволоки в зоне пучка")

    halo = [pe.withStroke(linewidth=3.0, foreground=SURFACE)]
    OFF = {"test_grid_33_11": ((0, 14), "center"), "specac": ((-12, -6), "right"),
           "purewave": ((13, 4), "left"), "test_grid_40_20": ((10, -14), "left")}
    for r in rows:
        off, ha = OFF[r["sample"]]
        ax.annotate(NAMES[r["sample"]], (r["x"], r["y"]), textcoords="offset points",
                    xytext=off, ha=ha, fontsize=8,
                    color=RED if r["excluded"] else INK2, path_effects=halo)

    ax.set_xlim(0.28, 0.62)
    ax.set_ylim(0.42, 0.80)
    ax.set_xlabel("$D/P$ — доля периода, занятая металлом", fontsize=10, color=INK)
    ax.set_ylabel("$D_{\\rm eff}/D_{\\rm phys}$", fontsize=10, color=INK)
    ax.set_title("Закон H5 с полным бюджетом погрешностей",
                 fontsize=11.5, color=INK, loc="left", fontweight="bold", pad=10)
    ax.grid(True, color=GRID, lw=0.8)
    ax.set_axisbelow(True)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    for s in ("left", "bottom"):
        ax.spines[s].set_color(GRID)
    ax.tick_params(colors=INK2, labelsize=9)
    ax.legend(loc="lower left", fontsize=8, framealpha=0.95, facecolor=SURFACE,
              edgecolor=GRID)
    fig.tight_layout()
    out = HERE / "r12_h5_errors.png"
    fig.savefig(out, dpi=200, facecolor=SURFACE)
    plt.close(fig)

    (RES / "a19_error_budget.json").write_text(json.dumps({
        "note": "конвенция границы — ОБЩАЯ систематика, вынесена в погрешность k, "
                "а не в разброс точек; ошибки по x и y коррелированы через D и через "
                "приколотый период",
        "owner_edge_shift_um": OWNER_EDGE_SHIFT_UM,
        "rows": rows,
        "fit_excl_purewave": fit, "fit_all": fit_all, "fit_owner_edge": fit_edge,
        "k_sys_edge": k_sys_edge,
    }, ensure_ascii=False, indent=2, default=float), encoding="utf-8")
    print(f"\n{out.name}, a19_error_budget.json")


if __name__ == "__main__":
    main()
