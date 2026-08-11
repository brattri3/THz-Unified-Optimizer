# -*- coding: utf-8 -*-
"""A19 — закон H5 в полном квадрате 0..1 с равным масштабом осей.

Смысл рисунка: при отношении осей 1:1 видно, что прямая идёт почти под 45°, то есть
угловой коэффициент близок к единице. Показаны ВСЕ точки, включая att-11-16 (16/11),
для которого микрофото нет и геометрия паспортная.

Подписи: <образец> <P>/<D> (D/P по паспорту ≠ D/P по измерению).

Бюджет погрешностей и подгонка берутся из `make_h5_errors.py` — здесь только рисование.

Запуск из корня репозитория:
  .venv\\Scripts\\python.exe research\\results\\a18_microscopy\\figs_report\\make_h5_unit.py
"""
import json
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patheffects as pe
from matplotlib.patches import Ellipse

from make_h5_errors import (budget, wfit, owner_edge_fit, EXCLUDED,
                            SURFACE, INK, INK2, MUTED, GRID, S1, S2, S3, RED)

HERE = Path(__file__).resolve().parent
RES = HERE.parent
REPO = HERE.parents[3]
A16 = REPO / "research" / "results" / "two_wgp" / "a16_h5_recompute.json"

# Паспорт: период/диаметр, мкм. specac и purewave — 25/10 (владелец, 2026-08-11).
PASSPORT = {
    "purewave":        (25.0, 10.0),
    "specac":          (25.0, 10.0),
    "test_grid_33_11": (33.0, 11.0),
    "test_grid_40_20": (40.0, 20.0),
}
NAMES = {"purewave": "purewave", "specac": "specac",
         "test_grid_33_11": "test_grid", "test_grid_40_20": "test_grid"}


def main():
    rows = budget()
    keep = [r for r in rows if not r["excluded"]]
    fit = wfit(keep)
    k, k_stat = fit["k"], fit["k_err"]
    k_sys = abs(owner_edge_fit(keep)["k"] - k)
    k_tot = float(np.hypot(k_stat, k_sys))

    ang = np.degrees(np.arctan(k))
    ang_lo = np.degrees(np.arctan(k - k_tot))
    ang_hi = np.degrees(np.arctan(k + k_tot))

    att = [p for p in json.loads(A16.read_text(encoding="utf-8"))["points"]
           if not p["measured"]]

    fig, ax = plt.subplots(figsize=(7.0, 7.0), facecolor=SURFACE)
    ax.set_facecolor(SURFACE)

    g = np.linspace(0, 1, 200)
    ax.plot(g, 1 - g, color=MUTED, ls=":", lw=1.4, label="$45^\\circ$, то есть $k=1$")
    ax.plot(g, 1 - k * g, color=S1, ls="--", lw=1.9,
            label=f"закон H5: $1-{k:.3f}\\,(D/P)$")
    ax.fill_between(g, 1 - (k - k_tot) * g, 1 - (k + k_tot) * g,
                    color=S1, alpha=0.12, lw=0,
                    label=f"$k = {k:.3f} \\pm {k_stat:.3f}_{{\\rm стат}} "
                          f"\\pm {k_sys:.3f}_{{\\rm конв}}$")
    ax.axvline(0.5, color=S2, ls="--", lw=1.1, alpha=0.7)
    ax.text(0.508, 0.965, "граница применимости Бланко $D/P<0.5$",
            fontsize=8, color=S2, va="top", rotation=90)

    for r in rows:
        c = RED if r["excluded"] else S1
        C = np.array([[r["sx"] ** 2, r["cov"]], [r["cov"], r["sy"] ** 2]])
        val, vec = np.linalg.eigh(C)
        a = np.degrees(np.arctan2(vec[1, np.argmax(val)], vec[0, np.argmax(val)]))
        w, h = 2 * np.sqrt(np.maximum(val[::-1], 0))
        ax.add_patch(Ellipse((r["x"], r["y"]), w, h, angle=a, facecolor=c,
                             alpha=0.25, edgecolor=c, lw=1.1, zorder=4))
        ax.scatter([r["x"]], [r["y"]], s=95 if r["excluded"] else 52,
                   marker="X" if r["excluded"] else "o", color=c, zorder=6)

    ax.scatter([p["D_over_P"] for p in att], [p["D_eff_over_D_phys"] for p in att],
               s=58, marker="s", color=S3, zorder=5,
               label="att-11-16: паспорт, микрофото нет")
    ax.scatter([], [], s=52, color=S1, label="измерено по микрофото, эллипс $1\\sigma$")
    ax.scatter([], [], s=95, marker="X", color=RED, label="исключён из подгонки")

    halo = [pe.withStroke(linewidth=3.0, foreground=SURFACE)]
    OFF = {"test_grid_33_11": ((-10, 10), "right"), "specac": ((-13, -2), "right"),
           "purewave": ((-13, -14), "right"), "test_grid_40_20": ((0, -18), "center")}
    for r in rows:
        P, D = PASSPORT[r["sample"]]
        txt = (f"{NAMES[r['sample']]} {P:.0f}/{D:.0f}\n"
               f"({D / P:.3f} $\\neq$ {r['x']:.3f})")
        off, ha = OFF[r["sample"]]
        ax.annotate(txt, (r["x"], r["y"]), textcoords="offset points", xytext=off,
                    ha=ha, fontsize=8, linespacing=1.25,
                    color=RED if r["excluded"] else INK2, path_effects=halo)
    ax.annotate("att 16/11\n(0.688 — паспорт,\nизмерения нет)", (0.710, 0.373),
                textcoords="offset points", xytext=(16, -12), ha="left", va="top",
                fontsize=8, color=INK2, linespacing=1.25, path_effects=halo)

    ax.annotate(f"наклон {ang:.1f}° (±{max(ang_hi - ang, ang - ang_lo):.1f}°)",
                (0.80, 1 - k * 0.80), textcoords="offset points", xytext=(-6, -16),
                ha="right", fontsize=9, color=S1, path_effects=halo)

    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.set_aspect("equal", adjustable="box")
    ax.set_xticks(np.arange(0, 1.01, 0.1))
    ax.set_yticks(np.arange(0, 1.01, 0.1))
    ax.set_xlabel("$D/P$ — доля периода, занятая металлом", fontsize=10, color=INK)
    ax.set_ylabel("$D_{\\rm eff}/D_{\\rm phys}$", fontsize=10, color=INK)
    ax.set_title("Закон H5 в равном масштабе осей: наклон почти 45°",
                 fontsize=11.5, color=INK, loc="left", fontweight="bold", pad=10)
    ax.grid(True, color=GRID, lw=0.7)
    ax.set_axisbelow(True)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    for s in ("left", "bottom"):
        ax.spines[s].set_color(GRID)
    ax.tick_params(colors=INK2, labelsize=9)
    ax.legend(loc="upper right", fontsize=8, framealpha=0.95, facecolor=SURFACE,
              edgecolor=GRID)
    fig.tight_layout()
    out = HERE / "r13_h5_unit_square.png"
    fig.savefig(out, dpi=200, facecolor=SURFACE)
    plt.close(fig)

    print(f"k = {k:.4f} ± {k_stat:.4f} (стат) ± {k_sys:.4f} (конвенция края)")
    print(f"полная погрешность k: ±{k_tot:.4f}")
    print(f"наклон: {ang:.2f}°  [{ang_lo:.2f}° .. {ang_hi:.2f}°]   до 45° "
          f"не хватает {45 - ang:.2f}° ({(45 - ang) / (ang_hi - ang):.1f} сигмы)")
    print(f"\n{'образец':<17}{'паспорт D/P':>13}{'измерено':>10}{'разница':>10}")
    for r in rows:
        P, D = PASSPORT[r["sample"]]
        print(f"{r['sample']:<17}{D / P:>13.3f}{r['x']:>10.3f}{r['x'] - D / P:>+10.3f}")
    print(f"\n{out.name}")


if __name__ == "__main__":
    main()
