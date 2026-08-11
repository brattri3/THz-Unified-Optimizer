# -*- coding: utf-8 -*-
"""A21b -- то же самое, что r15, но при ТЁМНОМ положении (~50°).

Контрольный эксперимент: если тренд «спад vs D/P» есть только при ярком
положении (0°) и исчезает при 50°, значит он физически реален, а не артефакт.

Запуск из корня репозитория:
  .venv\\Scripts\\python.exe research\\results\\a18_microscopy\\figs_report\\make_transmission_dark.py
"""
import re
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patheffects as pe
from matplotlib.lines import Line2D
import numpy as np

HERE  = Path(__file__).resolve().parent   # figs_report/
RES   = HERE.parent                       # a18_microscopy/
REPO  = RES.parents[2]
DATA  = REPO / "data_pool"

SURFACE = "#fcfcfb"
INK     = "#0b0b0b"
INK2    = "#52514e"
MUTED   = "#8a8a85"
GRID    = "#e3e3df"
S1      = "#2a78d6"
S3      = "#1baf7a"
RED     = "#c1272d"

import sys
sys.path.insert(0, str(REPO / "research" / "two_wgp"))
sys.path.insert(0, str(REPO / "research" / "experiments"))
sys.path.insert(0, str(REPO))
import a3_eps_stability as a3

BAND      = (0.3, 1.5)
TARGET    = 50.0     # целевой угол «тёмного» положения
TOL       = 10.0     # допуск ±10°

EXCLUDED = {"purewave"}   # purewave — дефектный образец, всё равно отмечаем

LABEL = {
    "purewave":        "purewave\n(55°)",
    "specac":          "specac\n(50°)",
    "test_grid_33_11": "33/11\n(50°)",
    "test_grid_40_20": "40/20\n(50°)",
}

# D/P из a19_error_budget (те же, что в r15)
import json
bud = {r["sample"]: r for r in
       json.loads((RES / "a19_error_budget.json").read_text(encoding="utf-8"))["rows"]}


def ax_style(ax):
    ax.set_facecolor(SURFACE)
    ax.grid(True, color=GRID, lw=0.8, zorder=0)
    ax.set_axisbelow(True)
    for sd in ("top", "right"):
        ax.spines[sd].set_visible(False)
    for sd in ("left", "bottom"):
        ax.spines[sd].set_color(GRID)
    ax.tick_params(colors=INK2, labelsize=9)


def get_point(ds, target=TARGET, tol=TOL):
    """Возвращает T_time и slope для угла ближайшего к target."""
    pool = DATA / ds
    best_diff, best = 999.0, None
    for sig_f in sorted(pool.glob("*_sig.txt")):
        m = re.search(r"_(-?\d+)deg_rep(\d+)", sig_f.name)
        if not m:
            continue
        angle = float(m.group(1))
        if abs(angle - target) > tol:
            continue
        bg_f = pool / sig_f.name.replace("_sig.txt", "_bg.txt")
        if not bg_f.exists():
            continue
        diff = abs(angle - target)
        if diff < best_diff:
            best_diff = diff
            best = (angle, sig_f, bg_f)
    if best is None:
        return None

    angle, sf, bf = best
    sig = np.loadtxt(sf)
    bg  = np.loadtxt(bf)
    T_time = (np.trapezoid(sig[:, 1]**2, sig[:, 0]) /
              np.trapezoid(bg[:, 1]**2, bg[:, 0]))

    # slope через DataManager (спектральный расчёт корректен при любом угле)
    dm = a3.DataManager(a3.config.DATA_DIR)
    dm.scan_directory()
    data_dm = dm.get_data_for_dataset(ds)
    closest = (min(data_dm.keys(), key=lambda a: abs(a - target))
               if data_dm else None)
    sl, sl_err = None, None
    if closest is not None and abs(closest - target) <= tol:
        ts, Es, tb, Eb = data_dm[closest]
        f, _, _, tr = a3.get_transmission_spectra(ts, Es, tb, Eb)
        f = np.asarray(f)
        T_sp = np.abs(tr)**2
        mask = (f >= BAND[0]) & (f <= BAND[1]) & np.isfinite(T_sp) & (T_sp > 0)
        fb, Tb2 = f[mask], T_sp[mask]
        if len(fb) > 2:
            db   = 10.0 * np.log10(Tb2)
            coef = np.polyfit(fb, db, 1)
            sl   = float(coef[0])
            resid  = db - np.polyval(coef, fb)
            s_res  = np.sqrt(np.sum(resid**2) / max(len(fb) - 2, 1))
            sl_err = float(s_res / np.sqrt(np.sum((fb - fb.mean())**2)))

    return dict(ds=ds, angle=angle, T=100.0 * T_time, slope=sl, slope_err=sl_err)


def main():
    SETS = ["purewave", "specac", "test_grid_33_11", "test_grid_40_20"]

    pts = []
    for ds in SETS:
        r = get_point(ds)
        if r is None:
            print(f"  {ds}: нет данных при ~{TARGET}°")
            continue
        x   = bud[ds]["x"]  if ds in bud else None
        sx  = bud[ds]["sx"] if ds in bud else 0.0
        if x is None:
            continue
        pts.append(dict(
            ds       = ds,
            angle    = r["angle"],
            x        = x,
            sx       = sx,
            excluded = (ds in EXCLUDED),
            T        = r["T"],
            sl       = r["slope"]    or 0.0,
            sl_err   = r["slope_err"] or 0.0,
        ))

    # ---------------------------------------------------------------- Figure -
    fig, (axT, axS) = plt.subplots(
        2, 1, figsize=(7.6, 7.2), sharex=True,
        facecolor=SURFACE, gridspec_kw={"hspace": 0.12},
    )
    for ax in (axT, axS):
        ax_style(ax)

    halo = [pe.withStroke(linewidth=3.2, foreground=SURFACE)]

    # Добавляем аннотацию «тёмное положение»
    for ax in (axT, axS):
        ax.annotate(
            "тёмное положение (~50°)\nконтрольный эксперимент",
            xy=(0.98, 0.97), xycoords="axes fraction",
            ha="right", va="top", fontsize=8.5, color=MUTED,
            style="italic",
        )

    # =================== Верхняя панель: T ===================
    axT.axhline(100, color=MUTED, ls=":", lw=1.2, zorder=1)
    axT.text(0.283, 100.6, "100 %", fontsize=8, color=MUTED, va="bottom")

    for p in pts:
        c  = RED if p["excluded"] else S1
        mk = "X" if p["excluded"] else "o"
        sz = 110 if p["excluded"] else 80
        if p["sx"] > 0:
            axT.plot([p["x"] - p["sx"], p["x"] + p["sx"]], [p["T"]] * 2,
                     color=c, lw=1.3, zorder=3)
        axT.scatter([p["x"]], [p["T"]], s=sz, marker=mk,
                    color=c, zorder=5, linewidths=0)

    axT.set_ylabel("интегральное пропускание, %", fontsize=10, color=INK)
    axT.set_title(
        "Контроль: пропускание и спад при тёмном положении (~50°)\n"
        "сравнение с ярким положением (0°) — тренд должен исчезнуть",
        fontsize=11.0, color=INK, loc="left", fontweight="bold", pad=10,
    )
    t_vals = [p["T"] for p in pts]
    margin = max(5, (max(t_vals) - min(t_vals)) * 0.6)
    axT.set_ylim(min(t_vals) - margin, max(t_vals) + margin)

    # =================== Нижняя панель: спад =================
    axS.axhline(0, color=MUTED, ls=":", lw=1.2, zorder=1)
    axS.text(0.283, 0.12, "0 — нет зависимости от ν", fontsize=8, color=MUTED, va="bottom")

    for p in pts:
        c  = RED if p["excluded"] else S1
        mk = "X" if p["excluded"] else "o"
        sz = 110 if p["excluded"] else 80
        axS.plot([p["x"]] * 2, [p["sl"] - p["sl_err"], p["sl"] + p["sl_err"]],
                 color=c, lw=1.5, zorder=3)
        if p["sx"] > 0:
            axS.plot([p["x"] - p["sx"], p["x"] + p["sx"]], [p["sl"]] * 2,
                     color=c, lw=1.2, zorder=3)
        axS.scatter([p["x"]], [p["sl"]], s=sz, marker=mk,
                    color=c, zorder=5, linewidths=0)

    axS.set_ylabel("спад пропускания, дБ/ТГц\n(< 0 — падает с ростом ν)", fontsize=10, color=INK)
    axS.set_xlabel("$D/P$ — плотность намотки (проволока / шаг)", fontsize=10, color=INK)
    sl_vals = [p["sl"] for p in pts]
    sl_errs = [p["sl_err"] for p in pts]
    y_lo = min(v - e for v, e in zip(sl_vals, sl_errs)) - 0.15
    y_hi = max(v + e for v, e in zip(sl_vals, sl_errs)) + 0.15
    axS.set_ylim(y_lo, y_hi)

    # =================== Подписи =====
    OFF_T = {
        "purewave":        (( 11,  5), "left"),
        "specac":          ((-11,  5), "right"),
        "test_grid_33_11": ((-11,  5), "right"),
        "test_grid_40_20": (( 11,  5), "left"),
    }
    OFF_S = {
        "purewave":        (( 11, -12), "left"),
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
                            ha=ha, fontsize=8, color=c, path_effects=halo,
                            linespacing=1.3)

    # =================== Легенда ====
    legend_handles = [
        Line2D([0], [0], marker="o", color="w", markerfacecolor=S1, markersize=8,
               label="образец (D/P из микрофото)"),
        Line2D([0], [0], marker="X", color="w", markerfacecolor=RED, markersize=9,
               label="исключён (purewave — дефектный)"),
    ]
    for ax in (axT, axS):
        ax.legend(handles=legend_handles,
                  loc="lower right" if ax is axT else "upper right",
                  fontsize=8.5, framealpha=0.96, facecolor=SURFACE, edgecolor=GRID)

    axS.set_xlim(0.27, 0.65)
    fig.tight_layout()

    out = HERE / "r16_transmission_dark_50deg.png"
    fig.savefig(out, dpi=200, facecolor=SURFACE)
    plt.close(fig)

    # =================== Таблица ====
    print(f"\n{'образец':<20} {'факт.угол':>10} {'D/P':>6} {'T,%':>7} {'спад дБ/ТГц':>14}")
    print("-" * 62)
    for p in sorted(pts, key=lambda q: q["x"]):
        mark = "  [ИСКЛ]" if p["excluded"] else ""
        print(f"{p['ds']:<20} {p['angle']:>10.1f} {p['x']:>6.3f} {p['T']:>7.1f}"
              f"  {p['sl']:>+6.3f} ± {p['sl_err']:.3f}{mark}")

    print(f"\nГотово: {out}")


if __name__ == "__main__":
    main()
