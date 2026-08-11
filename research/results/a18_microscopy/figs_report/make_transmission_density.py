# -*- coding: utf-8 -*-
"""A21 -- пропускание в ярком положении против плотности намотки, две панели.

Верх: интегральное пропускание во временной области (энергия импульса сигнала к
энергии опорного) -- величина без выбора полосы усреднения.
Низ: спад пропускания по частоте, дБ/ТГц -- «насколько сильно падает с частотой».

Запуск из корня репозитория:
  .venv\\Scripts\\python.exe research\\results\\a18_microscopy\\figs_report\\make_transmission_density.py
"""
import json
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patheffects as pe
from matplotlib.lines import Line2D

HERE = Path(__file__).resolve().parent          # figs_report/
RES  = HERE.parent                              # a18_microscopy/

# --- Палитра (светлый фон) ---
SURFACE = "#fcfcfb"
INK     = "#0b0b0b"
INK2    = "#52514e"
MUTED   = "#8a8a85"
GRID    = "#e3e3df"
S1      = "#2a78d6"   # измеренный образец (синий)
S3      = "#1baf7a"   # аттенюаторы (паспортная геометрия, зелёный)
RED     = "#c1272d"   # исключённые

EXCLUDED = {"purewave"}

LABEL = {
    "purewave":        "purewave",
    "specac":          "specac",
    "test_grid_33_11": "33/11",
    "test_grid_40_20": "40/20",
}

ATT_DP = {
    "att-11-16-356": 11.0 / 15.5,
    "att-11-16-s1":  11.0 / 16.0,
    "att-11-16-s2":  11.0 / 16.0,
    "att-11-16-s3":  11.0 / 16.0,
}


def lin(x, y):
    b, a = np.polyfit(x, y, 1)
    r    = y - (a + b * x)
    n    = len(x)
    s    = np.sqrt(np.sum(r ** 2) / max(n - 2, 1))
    ss_x = np.sum((x - x.mean()) ** 2)
    b_err = float(s / np.sqrt(ss_x)) if ss_x else np.nan
    ss_y  = np.sum((y - y.mean()) ** 2)
    r2    = float(1 - np.sum(r ** 2) / ss_y) if ss_y else np.nan
    return float(a), float(b), b_err, r2


def ax_style(ax):
    ax.set_facecolor(SURFACE)
    ax.grid(True, color=GRID, lw=0.8, zorder=0)
    ax.set_axisbelow(True)
    for sd in ("top", "right"):
        ax.spines[sd].set_visible(False)
    for sd in ("left", "bottom"):
        ax.spines[sd].set_color(GRID)
    ax.tick_params(colors=INK2, labelsize=9)


def main():
    tr  = {r["dataset"]: r for r in
           json.loads((RES / "a21_transmission.json").read_text(encoding="utf-8"))["rows"]}
    bud = {r["sample"]:  r for r in
           json.loads((RES / "a19_error_budget.json").read_text(encoding="utf-8"))["rows"]}

    pts = []
    for ds, r in tr.items():
        if ds in bud:
            x, sx, kind = bud[ds]["x"], bud[ds]["sx"], "meas"
        elif ds in ATT_DP:
            x, sx, kind = ATT_DP[ds], 0.0, "att"
        else:
            continue

        pts.append(dict(
            ds=ds, x=x, sx=sx, kind=kind, excluded=(ds in EXCLUDED),
            T=r["T_time_pct"],
            Tb=r["T_band_median_pct"],
            sl=r["slope_db_per_thz"],
            sl_err=r["slope_err_db_per_thz"],
        ))

    keep = [p for p in pts if not p["excluded"]]
    Xk   = np.array([p["x"]  for p in keep])
    aT, bT, bT_err, rT = lin(Xk, np.array([p["T"]  for p in keep]))
    aS, bS, bS_err, rS = lin(Xk, np.array([p["sl"] for p in keep]))

    fig, (axT, axS) = plt.subplots(
        2, 1, figsize=(7.6, 7.2), sharex=True,
        facecolor=SURFACE, gridspec_kw={"hspace": 0.12},
    )
    for ax in (axT, axS):
        ax_style(ax)

    g = np.linspace(0.27, 0.77, 120)
    halo = [pe.withStroke(linewidth=3.2, foreground=SURFACE)]

    # =================== Верхняя панель: уровень ===================
    axT.axhline(100, color=MUTED, ls=":", lw=1.2, zorder=1)
    axT.text(0.283, 100.4, "100 %", fontsize=8, color=MUTED)
    axT.plot(g, aT + bT * g, color=S1, ls="--", lw=1.8, zorder=2,
             label=f"тренд: ${bT:+.1f}\\pm{bT_err:.1f}$ %/ед. $D/P$,  $R^2={rT:.2f}$")

    for p in pts:
        c  = RED if p["excluded"] else (S3 if p["kind"] == "att" else S1)
        mk = "X" if p["excluded"] else ("s" if p["kind"] == "att" else "o")
        if p["sx"] > 0:
            axT.plot([p["x"] - p["sx"], p["x"] + p["sx"]], [p["T"]] * 2, color=c, lw=1.3)
        axT.scatter([p["x"]], [p["T"]], s=105 if p["excluded"] else 65,
                    marker=mk, color=c, zorder=5, linewidths=0)

    axT.set_ylabel("интегральное пропускание, %", fontsize=10, color=INK)
    axT.set_title(
        "Уровень пропускания не зависит от $D/P$,\n"
        "а его спад по частоте — зависит",
        fontsize=11.5, color=INK, loc="left", fontweight="bold", pad=10,
    )
    axT.set_ylim(87, 113)

    # =================== Нижняя панель: спад по частоте ===========
    axS.axhline(0, color=MUTED, ls=":", lw=1.2, zorder=1)
    axS.text(0.283, 0.06, "0 — нет зависимости от частоты", fontsize=8, color=MUTED)
    axS.plot(g, aS + bS * g, color=S1, ls="--", lw=1.8, zorder=2,
             label=f"тренд: ${bS:+.3f}\\pm{bS_err:.3f}$ (дБ/ТГц)/(ед. $D/P$),  $R^2={rS:.2f}$")

    for p in pts:
        c  = RED if p["excluded"] else (S3 if p["kind"] == "att" else S1)
        mk = "X" if p["excluded"] else ("s" if p["kind"] == "att" else "o")
        # вертикальная погрешность наклона
        axS.plot([p["x"]] * 2, [p["sl"] - p["sl_err"], p["sl"] + p["sl_err"]],
                 color=c, lw=1.5, zorder=3)
        # горизонтальная погрешность D/P
        if p["sx"] > 0:
            axS.plot([p["x"] - p["sx"], p["x"] + p["sx"]], [p["sl"]] * 2, color=c, lw=1.2)
        axS.scatter([p["x"]], [p["sl"]], s=105 if p["excluded"] else 65,
                    marker=mk, color=c, zorder=5, linewidths=0)

    axS.set_ylabel("спад пропускания, дБ/ТГц\n(< 0 — падает с ростом частоты)", fontsize=10, color=INK)
    axS.set_xlabel("$D/P$ — плотность намотки (проволока / шаг)", fontsize=10, color=INK)
    axS.set_ylim(-0.95, 0.72)

    # =================== Подписи образцов =========================
    OFF_T = {
        "purewave":        (( 11,   4), "left"),
        "specac":          ((-11,   4), "right"),
        "test_grid_33_11": ((-11,   4), "right"),
        "test_grid_40_20": (( 11,   4), "left"),
    }
    OFF_S = {
        "purewave":        (( 11,  -8), "left"),
        "specac":          ((-11,   5), "right"),
        "test_grid_33_11": ((-11,   5), "right"),
        "test_grid_40_20": (( 11,   5), "left"),
    }
    for p in pts:
        if p["ds"] not in LABEL:
            continue
        c = RED if p["excluded"] else INK2
        for ax, field, ofd in ((axT, "T", OFF_T), (axS, "sl", OFF_S)):
            if p["ds"] in ofd:
                off, ha = ofd[p["ds"]]
                ax.annotate(LABEL[p["ds"]], (p["x"], p[field]),
                            textcoords="offset points", xytext=off,
                            ha=ha, fontsize=8.5, color=c, path_effects=halo)

    # Групповая подпись аттенюаторов
    att_pts = [p for p in pts if p["kind"] == "att"]
    if att_pts:
        ax_x = np.mean([p["x"] for p in att_pts])
        for ax, field in ((axT, "T"), (axS, "sl")):
            yv = np.mean([p[field] for p in att_pts])
            ax.annotate("att-11-16\n(4 уст.)", (ax_x, yv),
                        textcoords="offset points", xytext=(-14, -20),
                        ha="right", va="top", fontsize=8, color=INK2,
                        linespacing=1.25, path_effects=halo)

    # =================== Общая легенда маркеров ===================
    legend_handles = [
        Line2D([0], [0], marker="o", color="w", markerfacecolor=S1, markersize=8,
               label="образец (D/P из микрофото)"),
        Line2D([0], [0], marker="s", color="w", markerfacecolor=S3, markersize=8,
               label="ATT-11-16 (D/P паспорт)"),
        Line2D([0], [0], marker="X", color="w", markerfacecolor=RED, markersize=9,
               label="исключён из тренда (purewave)"),
    ]
    # Добавляем линию тренда к легенде верхней панели
    trend_handle, _ = axT.get_legend_handles_labels()
    axT.legend(handles=legend_handles + trend_handle,
               loc="lower right", fontsize=8.5, framealpha=0.96,
               facecolor=SURFACE, edgecolor=GRID)
    axS.legend(loc="upper right", fontsize=8.5, framealpha=0.96,
               facecolor=SURFACE, edgecolor=GRID)

    axS.set_xlim(0.27, 0.77)
    fig.tight_layout()
    out = HERE / "r15_transmission_vs_density.png"
    fig.savefig(out, dpi=200, facecolor=SURFACE)
    plt.close(fig)

    # =================== Таблица в консоль ========================
    print(f"\n{'образец':<20} {'D/P':>6} {'T врем,%':>9} {'T полос,%':>10} {'спад дБ/ТГц':>14}")
    print("-" * 66)
    for p in sorted(pts, key=lambda q: q["x"]):
        mark = "  [ИСКЛ]" if p["excluded"] else ("  [паспорт]" if p["kind"] == "att" else "")
        print(f"{p['ds']:<20} {p['x']:>6.3f} {p['T']:>9.1f} {p['Tb']:>10.1f}"
              f"  {p['sl']:>+6.3f} +/- {p['sl_err']:.3f}{mark}")

    print(f"\nТренд уровень : {bT:+.2f} +/- {bT_err:.2f} %/(ед. D/P),           R2 = {rT:.3f}")
    print(f"Тренд спад    : {bS:+.3f} +/- {bS_err:.3f} (дБ/ТГц)/(ед. D/P),  R2 = {rS:.3f}")
    print(f"\nГотово: {out}")


if __name__ == "__main__":
    main()
