# -*- coding: utf-8 -*-
"""A20 — медианный коэффициент экстинкции против плотности намотки D/P.

Экстинкция — из `a20_extinction.json` (посчитана прямо из угловых сканов).
Плотность: у четырёх образцов — измеренная по микрофото (с погрешностью из
`a19_error_budget.json`), у att-11-16 — паспортная, погрешности нет.

Вертикальный отрезок у точки — НЕ погрешность, а размах экстинкции по полосе:
от значения на верхнем крае полосы до максимума. Экстинкция сильно зависит от
частоты, и одной медианой это скрывать нельзя.

Запуск из корня репозитория:
  .venv\\Scripts\\python.exe research\\results\\a18_microscopy\\figs_report\\make_er_density.py
"""
import json
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patheffects as pe

HERE = Path(__file__).resolve().parent
RES = HERE.parent

SURFACE, INK, INK2, MUTED, GRID = "#fcfcfb", "#0b0b0b", "#52514e", "#8a8a85", "#e3e3df"
S1, S2, S3, RED = "#2a78d6", "#eb6834", "#1baf7a", "#c1272d"

EXCLUDED = {"purewave"}
LABEL = {"purewave": "purewave", "specac": "specac",
         "test_grid_33_11": "test_grid 33/11", "test_grid_40_20": "test_grid 40/20"}


def main():
    er = {r["dataset"]: r for r in
          json.loads((RES / "a20_extinction.json").read_text(encoding="utf-8"))["rows"]}
    bud = {r["sample"]: r for r in
           json.loads((RES / "a19_error_budget.json").read_text(encoding="utf-8"))["rows"]}
    att_dp = {"att-11-16-356": 11.0 / 15.5, "att-11-16-s1": 11.0 / 16.0,
              "att-11-16-s2": 11.0 / 16.0, "att-11-16-s3": 11.0 / 16.0}

    pts = []
    for ds, r in er.items():
        if ds in bud:
            x, sx, kind = bud[ds]["x"], bud[ds]["sx"], "meas"
        else:
            x, sx, kind = att_dp[ds], 0.0, "att"
        pts.append(dict(ds=ds, x=x, sx=sx, kind=kind,
                        y=r["ER_db_median_band"], y_hi=r["ER_db_max"],
                        y_lo=r["ER_db_at_band_edges"][1],
                        excluded=ds in EXCLUDED))

    fit_pts = [p for p in pts if not p["excluded"]]
    X = np.array([p["x"] for p in fit_pts])
    Y = np.array([p["y"] for p in fit_pts])
    b, a = np.polyfit(X, Y, 1)
    resid = Y - (a + b * X)
    n = len(X)
    s = float(np.sqrt(np.sum(resid ** 2) / (n - 2)))
    b_err = float(s / np.sqrt(np.sum((X - X.mean()) ** 2)))
    rr = float(1 - np.sum(resid ** 2) / np.sum((Y - Y.mean()) ** 2))

    fig, ax = plt.subplots(figsize=(7.2, 4.8), facecolor=SURFACE)
    ax.set_facecolor(SURFACE)
    g = np.linspace(0.29, 0.75, 100)
    ax.plot(g, a + b * g, color=S1, ls="--", lw=1.8,
            label=f"линейный тренд: {b:.0f} ± {b_err:.0f} дБ на единицу $D/P$"
                  f"  ($R^2={rr:.2f}$)")

    for p in pts:
        c = RED if p["excluded"] else (S3 if p["kind"] == "att" else S1)
        ax.plot([p["x"], p["x"]], [p["y_lo"], p["y_hi"]], color=c, lw=1.2,
                alpha=0.45, zorder=3, solid_capstyle="round")
        if p["sx"] > 0:
            ax.plot([p["x"] - p["sx"], p["x"] + p["sx"]], [p["y"], p["y"]],
                    color=c, lw=1.2, alpha=0.75, zorder=3)
        ax.scatter([p["x"]], [p["y"]],
                   s=110 if p["excluded"] else (58 if p["kind"] == "att" else 66),
                   marker="X" if p["excluded"] else ("s" if p["kind"] == "att" else "o"),
                   color=c, zorder=6)

    ax.scatter([], [], s=66, color=S1, label="измерено по микрофото")
    ax.scatter([], [], s=58, marker="s", color=S3, label="att-11-16: $D/P$ по паспорту")
    ax.scatter([], [], s=110, marker="X", color=RED, label="purewave: исключён из тренда")
    ax.plot([], [], color=MUTED, lw=1.2, alpha=0.6,
            label="вертикаль — размах по полосе 0.3–1.5 ТГц")

    halo = [pe.withStroke(linewidth=3.0, foreground=SURFACE)]
    OFF = {"purewave": ((13, 2), "left"), "specac": ((-12, 2), "right"),
           "test_grid_33_11": ((12, -2), "left"), "test_grid_40_20": ((0, -16), "center")}
    for p in pts:
        if p["ds"] in OFF:
            off, ha = OFF[p["ds"]]
            ax.annotate(LABEL[p["ds"]], (p["x"], p["y"]), textcoords="offset points",
                        xytext=off, ha=ha, fontsize=8, path_effects=halo,
                        color=RED if p["excluded"] else INK2)
    ax.annotate("att-11-16\n(4 установки)", (0.688, 39.0), textcoords="offset points",
                xytext=(-14, 6), ha="right", fontsize=8, color=INK2,
                linespacing=1.2, path_effects=halo)

    ax.set_xlim(0.29, 0.76)
    ax.set_ylim(22, 50)
    ax.set_xlabel("$D/P$ — плотность намотки (доля периода, занятая металлом)",
                  fontsize=10, color=INK)
    ax.set_ylabel("коэффициент экстинкции, дБ", fontsize=10, color=INK)
    ax.set_title("Чем плотнее намотка, тем глубже гашение",
                 fontsize=11.5, color=INK, loc="left", fontweight="bold", pad=10)
    ax.grid(True, color=GRID, lw=0.8)
    ax.set_axisbelow(True)
    for sd in ("top", "right"):
        ax.spines[sd].set_visible(False)
    for sd in ("left", "bottom"):
        ax.spines[sd].set_color(GRID)
    ax.tick_params(colors=INK2, labelsize=9)
    ax.legend(loc="upper left", fontsize=8, framealpha=0.95, facecolor=SURFACE,
              edgecolor=GRID)
    fig.tight_layout()
    out = HERE / "r14_er_vs_density.png"
    fig.savefig(out, dpi=200, facecolor=SURFACE)
    plt.close(fig)

    print(f"{'образец':<17}{'D/P':>8}{'медиана ER':>12}{'размах':>16}")
    for p in sorted(pts, key=lambda q: q["x"]):
        mark = "  ИСКЛ" if p["excluded"] else ("  паспорт" if p["kind"] == "att" else "")
        print(f"{p['ds']:<17}{p['x']:>8.3f}{p['y']:>12.1f}"
              f"{p['y_lo']:>9.1f}..{p['y_hi']:<6.1f}{mark}")
    print(f"\nтренд: ER[дБ] = {a:.1f} + {b:.1f}*(D/P),  наклон {b:.1f} ± {b_err:.1f} дБ, "
          f"R^2 = {rr:.3f}, n = {n}")
    print(f"purewave: измерено {er['purewave']['ER_db_median_band']:.1f} дБ, "
          f"тренд при D/P={bud['purewave']['x']:.3f} даёт "
          f"{a + b * bud['purewave']['x']:.1f} дБ")
    print(f"\n{out.name}")


if __name__ == "__main__":
    main()
