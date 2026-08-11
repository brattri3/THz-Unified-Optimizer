# -*- coding: utf-8 -*-
"""A19: пересчёт графика D_eff/D_phys против D/P на геометрии из микрофото A18.

Что меняется по сравнению с A16: только ГЕОМЕТРИЯ (P и D из новых микрофото).
Значения D_eff НЕ пересчитываются — берутся из a16_h5_recompute.json (поле D_eff_used),
то есть из подгонок, где период приколот к прежнему значению. Для purewave и
test_grid_40_20 период сдвинулся на 4 %, поэтому их D_eff предварителен: точки
рисуются открытыми маркерами.

Точки att-11-16-* не трогаются: новых микрофото для прибора нет, геометрия паспортная.

Запуск из корня репозитория:
  .venv\\Scripts\\python.exe research\\results\\a18_microscopy\\figs_report\\make_h5_dp.py
"""
import json
import sys
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patheffects as pe

HERE = Path(__file__).resolve().parent
RES = HERE.parent                      # research/results/a18_microscopy
REPO = HERE.parents[3]                 # корень репозитория
A16 = REPO / "research" / "results" / "two_wgp" / "a16_h5_recompute.json"

SURFACE, INK, INK2, MUTED, GRID = "#fcfcfb", "#0b0b0b", "#52514e", "#8a8a85", "#e3e3df"
S1, S2, S3, S4 = "#2a78d6", "#eb6834", "#1baf7a", "#eda100"
RED = "#c1272d"   # только для исключённой точки; форма маркера тоже отличается

# Образцы, исключённые из подгонки закона H5, и причина (владелец, 2026-08-11).
EXCLUDED = {
    "purewave": "восстановленный поляризатор: в средней части проволока провисла "
                "из-за отклеивания от оправы. Микрофото сняты в сохранившейся "
                "заводской зоне, а пучок проходит через провисшую середину — "
                "измеренная геометрия не описывает то, что видит THz-луч",
}

# Образцы, у которых период сдвинулся настолько, что D_eff требует перезапуска подгонки
NEEDS_REFIT = {"purewave", "test_grid_40_20"}
# Сдвиг диаметра при конвенции края владельца (f ~ 0.25) относительно полувысоты, мкм
OWNER_EDGE_SHIFT_UM = 0.6

NAMES_RU = {
    "purewave": "purewave",
    "specac": "specac",
    "test_grid_33_11": "test_grid 33/11",
    "test_grid_40_20": "test_grid 40/20",
}


def load_a18():
    """Геометрия из артефактов A18: период по плотности, диаметр — медиана по образцу."""
    out = {}
    for name in NAMES_RU:
        p = RES / f"a18_{name}.json"
        d = json.loads(p.read_text(encoding="utf-8"))
        out[name] = {
            "P": d["P_density"]["mean"],
            "P_se": d["P_density"]["sem"],
            "D": d["pooled"]["D"]["median"],
            "D_sd": d["pooled"]["D"].get("std", float("nan")),
        }
    return out


def fit_k(x, y):
    """МНК для однопараметрической формы y = 1 - k*x (та же форма, что в A16)."""
    x, y = np.asarray(x, float), np.asarray(y, float)
    k = float(np.sum(x * (1.0 - y)) / np.sum(x * x))
    r = y - (1.0 - k * x)
    rms = float(np.sqrt(np.mean(r ** 2)))
    n = len(x)
    s = float(np.sqrt(np.sum(r ** 2) / max(n - 1, 1)))
    k_err = float(s / np.sqrt(np.sum(x * x)))
    return k, k_err, rms


def main():
    # --final: только итоговое расположение точек, без стрелок «было -> стало»
    HISTORY = "--final" not in sys.argv
    a16 = json.loads(A16.read_text(encoding="utf-8"))
    deff = a16["D_eff_used"]
    geo16 = a16["geometry"]
    att = [p for p in a16["points"] if not p["measured"]]
    a18 = load_a18()

    rows = []
    for name, g in a18.items():
        D_eff = deff[name]
        rows.append({
            "sample": name,
            "P_um": g["P"], "P_se_um": g["P_se"], "D_um": g["D"], "D_sd_um": g["D_sd"],
            "D_eff_um": D_eff,
            "D_over_P": g["D"] / g["P"],
            "D_eff_over_D_phys": D_eff / g["D"],
            "old_D_over_P": geo16[name]["D_um"] / geo16[name]["P_um"],
            "old_ratio": D_eff / geo16[name]["D_um"],
            "deff_provisional": name in NEEDS_REFIT,
            # чувствительность к конвенции края (вопрос A-1)
            "D_over_P_owner_edge": (g["D"] + OWNER_EDGE_SHIFT_UM) / g["P"],
            "ratio_owner_edge": D_eff / (g["D"] + OWNER_EDGE_SHIFT_UM),
        })

    # Если перефит на новой геометрии уже выполнен (A19), D_eff берём оттуда:
    # период приколот, поэтому арифметический пересчёт неправомерен.
    refit_path = RES / "a19_refit.json"
    if refit_path.exists():
        refit = {r["sample"]: r
                 for r in json.loads(refit_path.read_text(encoding="utf-8"))["rows"]}
        for r in rows:
            f = refit[r["sample"]]
            r["D_eff_um"] = f["D_eff_um"]
            r["D_eff_over_D_phys"] = f["ratio_new"]
            r["ratio_owner_edge"] = f["D_eff_um"] / (r["D_um"] + OWNER_EDGE_SHIFT_UM)
            r["deff_provisional"] = False
            r["deff_refitted"] = True

    keep = [r for r in rows if r["sample"] not in EXCLUDED]
    xm = [r["D_over_P"] for r in keep]
    ym = [r["D_eff_over_D_phys"] for r in keep]
    xa = xm + [p["D_over_P"] for p in att]
    ya = ym + [p["D_eff_over_D_phys"] for p in att]
    xall = [r["D_over_P"] for r in rows]
    yall = [r["D_eff_over_D_phys"] for r in rows]
    xo = [r["old_D_over_P"] for r in rows]
    yo = [r["old_ratio"] for r in rows]
    xe = [r["D_over_P_owner_edge"] for r in keep]
    ye = [r["ratio_owner_edge"] for r in keep]

    fits = {
        "A19_measured_only": dict(zip(("k", "k_err", "rms"), fit_k(xm, ym))),
        "A19_measured_incl_purewave": dict(zip(("k", "k_err", "rms"), fit_k(xall, yall))),
        "A19_all_incl_passport": dict(zip(("k", "k_err", "rms"), fit_k(xa, ya))),
        "A19_measured_owner_edge": dict(zip(("k", "k_err", "rms"), fit_k(xe, ye))),
        "A16_measured_only": a16["h5_fits"]["measured_only"],
        "A16_all_incl_passport": a16["h5_fits"]["all_incl_passport"],
    }

    # ------------------------------- рисунок -------------------------------
    fig, ax = plt.subplots(figsize=(7.2, 4.8), facecolor=SURFACE)
    ax.set_facecolor(SURFACE)

    grid = np.linspace(0.28, 0.76, 200)
    k16 = a16["h5_fits"]["measured_only"]["k"]
    k19 = fits["A19_measured_only"]["k"]
    if HISTORY:
        ax.plot(grid, 1 - k16 * grid, color=MUTED, ls=":", lw=1.6,
                label=f"A16, прежняя геометрия: $1-{k16:.3f}\\,(D/P)$")
    ax.plot(grid, 1 - k19 * grid, color=S1, ls="--", lw=1.8,
            label=("A19, микрофото: " if HISTORY else "закон H5: ")
                  + f"$1-{k19:.3f}\\,(D/P)$")
    ax.axvline(0.5, color=S2, ls="--", lw=1.2, alpha=0.75)
    ax.text(0.505, 0.845 if HISTORY else 0.775,
            "граница применимости\nмодели Бланко $D/P<0.5$",
            fontsize=8, color=S2, va="top")

    if HISTORY:
        # смещение точек: было -> стало
        for r in rows:
            ax.annotate("", xy=(r["D_over_P"], r["D_eff_over_D_phys"]),
                        xytext=(r["old_D_over_P"], r["old_ratio"]),
                        arrowprops=dict(arrowstyle="->", color=MUTED, lw=1.0,
                                        alpha=0.85, shrinkA=3, shrinkB=5))
        ax.scatter(xo, yo, s=34, facecolors="none", edgecolors=MUTED, lw=1.2,
                   zorder=3, label="прежняя геометрия (A16)")

    excl = [r for r in rows if r["sample"] in EXCLUDED]
    solid = [r for r in keep if not r["deff_provisional"]]
    prov = [r for r in keep if r["deff_provisional"]]
    ax.scatter([r["D_over_P"] for r in solid], [r["D_eff_over_D_phys"] for r in solid],
               s=68, color=S1, zorder=5,
               label=("геометрия из микрофото, подгонка перезапущена" if not HISTORY
                      else "микрофото A18 + перефит"))
    if prov:
        ax.scatter([r["D_over_P"] for r in prov], [r["D_eff_over_D_phys"] for r in prov],
                   s=68, facecolors=SURFACE, edgecolors=S1, lw=2.0, zorder=5,
                   label="то же, $D_{\\rm eff}$ ждёт перезапуска подгонки" if not HISTORY
                         else "микрофото A18, $D_{\\rm eff}$ ждёт перезапуска подгонки")
    # Исключённая точка: красный + другая форма маркера (не только цветом)
    ax.scatter([r["D_over_P"] for r in excl], [r["D_eff_over_D_phys"] for r in excl],
               s=110, marker="X", color=RED, zorder=6,
               label="исключён: провисание проволоки в зоне пучка")
    ax.scatter([p["D_over_P"] for p in att], [p["D_eff_over_D_phys"] for p in att],
               s=52, marker="s", color=S3, zorder=4,
               label="att-11-16 (паспорт, микрофото нет)")

    # смещения подписей подобраны вручную: точки и стрелки идут плотно
    LABEL_OFF = {
        "purewave":        ((12, 7), "left"),
        "specac":          ((-10, -13), "right"),
        # без стрелок место сверху свободно; со стрелками там проходит линия A16
        "test_grid_33_11": ((-9, 9), "right") if HISTORY else ((0, 12), "center"),
        "test_grid_40_20": ((10, -13), "left"),
    }
    halo = [pe.withStroke(linewidth=3.0, foreground=SURFACE)]
    for r in rows:
        off, ha = LABEL_OFF[r["sample"]]
        ax.annotate(NAMES_RU[r["sample"]], (r["D_over_P"], r["D_eff_over_D_phys"]),
                    textcoords="offset points", xytext=off, ha=ha, fontsize=8,
                    color=RED if r["sample"] in EXCLUDED else INK2,
                    path_effects=halo)
    ax.annotate("att-11-16\n(4 установки)", (0.695, 0.418), textcoords="offset points",
                xytext=(14, 0), ha="left", va="center", fontsize=8, color=INK2,
                linespacing=1.15, path_effects=halo)

    ax.set_xlim(0.28, 0.78)
    ax.set_ylim(0.30, 0.90 if HISTORY else 0.80)
    ax.set_xlabel("$D/P$ — доля периода, занятая металлом", fontsize=10, color=INK)
    ax.set_ylabel("$D_{\\rm eff}/D_{\\rm phys}$", fontsize=10, color=INK)
    ax.set_title("Дефицит эффективного диаметра на уточнённой геометрии" if HISTORY
                 else "Чем плотнее решётка, тем меньше эффективный диаметр",
                 fontsize=11.5, color=INK, loc="left", fontweight="bold", pad=10)
    ax.grid(True, color=GRID, lw=0.8)
    ax.set_axisbelow(True)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    for side in ("left", "bottom"):
        ax.spines[side].set_color(GRID)
    ax.tick_params(colors=INK2, labelsize=9)
    ax.legend(loc="lower left", fontsize=8, framealpha=0.95, facecolor=SURFACE,
              edgecolor=GRID)

    fig.tight_layout()
    out_png = HERE / ("r11_h5_dp.png" if HISTORY else "r11b_h5_dp_final.png")
    fig.savefig(out_png, dpi=200, facecolor=SURFACE)
    plt.close(fig)

    payload = {
        "source": "геометрия A18 (микрофото 2026-08-10) + D_eff из a16_h5_recompute.json",
        "caveat": "D_eff НЕ пересчитывался; для purewave и test_grid_40_20 период сдвинулся "
                  "на 4 %, их D_eff требует перезапуска подгонки (вопрос A-2)",
        "owner_edge_shift_um": OWNER_EDGE_SHIFT_UM,
        "excluded_from_fit": EXCLUDED,
        "points": rows,
        "att_points": att,
        "h5_fits": fits,
    }
    out_json = RES / "a19_h5_dp.json"
    out_json.write_text(json.dumps(payload, ensure_ascii=False, indent=2), encoding="utf-8")

    print(f"{'образец':<18}{'D/P было':>10}{'стало':>9}{'ratio было':>12}{'стало':>9}")
    for r in rows:
        flag = " *" if r["deff_provisional"] else ""
        print(f"{r['sample']:<18}{r['old_D_over_P']:>10.3f}{r['D_over_P']:>9.3f}"
              f"{r['old_ratio']:>12.3f}{r['D_eff_over_D_phys']:>9.3f}{flag}")
    print("\nзакон H5  y = 1 - k*(D/P):")
    for tag, f in fits.items():
        print(f"  {tag:<26} k = {f['k']:.4f} ± {f.get('k_err', float('nan')):.4f}"
              f"   rms = {f['rms']:.4f}")
    print(f"\n{out_png.name}, {out_json.name}")


if __name__ == "__main__":
    main()
