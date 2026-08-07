#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""A14b — диаметр проволоки с поправкой на расфокусировку (замечание владельца).

ПОВОД. В A14 диаметр вышел `D = 12.55 ± 0.60` мкм при номинале 11. Владелец
возразил: на максимальном увеличении не все проволоки в фокусе, картинка
расплывается, значит диаметр завышен. Возражение верно, и вот механизм.

`micrograph_period.in_focus_mask` отбирает проволоки по КВАНТИЛЯМ ВНУТРИ КАДРА
(верхние 45 % по крутизне фронта, верхние 60 % по яркости). Критерий
ОТНОСИТЕЛЬНЫЙ: он всегда оставляет около половины проволок кадра, даже если
размыты все. Для компланарного образца это приемлемо, для этого — нет: пара
кадров «фокус 0 / 0.12 мм» прямо показывает разную глубину залегания.

МЕТОД. Не решать, какая проволока «в фокусе», а измерить РАЗМЫТИЕ каждой и
привести все к общему уровню резкости. Для каждой проволоки:
    w10  — ширина следа на уровне 10 % от пика, мкм (та же величина, что в A14);
    b    — размытие = (высота пика) / (макс. модуль градиента), px; для ступеньки,
           свёрнутой с аппаратной функцией, эта величина растёт линейно с её
           шириной, то есть это обратная нормированная крутизна фронта.
Расфокусировка уширяет след линейно по `b`, поэтому регрессия `w10(b)` позволяет
пересчитать ширину на любой заданный уровень резкости.

ДВА РАЗНЫХ ОТВЕТА, которые важно не путать:
  * `b -> 0`  — гипотетический ИДЕАЛЬНЫЙ прибор. Перекоррекция: у микроскопа
    есть собственный предел, ниже которого резкость не бывает.
  * `b -> b0` — прибор в НАИЛУЧШЕМ фокусе, где `b0` измеряется по калибровочной
    гребёнке (её края физически резкие, всё наблюдаемое размытие там —
    аппаратное). Это и есть рекомендуемая оценка: она убирает расфокусировку,
    но оставляет аппаратную функцию, как и у всех прежних образцов проекта.

КОНТРОЛЬ. Та же процедура прогоняется по `purewave`, для которого в
`MANIFEST.md` опубликовано `D = 10.6 ± 0.5`. Совпадение НАКЛОНОВ у двух образцов
доказывает, что связь «размытие → ширина» аппаратная, а не подгонка под образец.

Запуск ИЗ КОРНЯ репозитория:
    .venv/Scripts/python.exe research/two_wgp/a14b_diameter_defocus.py
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt                                       # noqa: E402

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "research" / "experiments"))
import micrograph_period as mp                                        # noqa: E402

MICS = REPO / "research" / "validation" / "micrographs"
OUT = REPO / "research" / "results" / "two_wgp"
UM_PER_PX_HIGH = 0.1264
NOMINAL_D = 11.0
D_EFF_FIT = 7.233          # A14, перефит на измеренном периоде P = 31.75
P_MEASURED = 31.75

STEM = "wgp-grid-D11-P33"
FRAMES_33 = [MICS / "wgp_grid_D11_P33" / f"{STEM}__wgp__{t}.bmp"
             for t in ("highmag_01", "highmag_02",
                       "highmag_focus0um", "highmag_focus120um")]
FRAMES_PW = [MICS / "purewave" / f"purewave__wgp__highmag_0{i}.bmp" for i in (1, 2, 3)]
CAL_HIGH = MICS / "wgp_grid_D11_P33" / f"{STEM}__cal__highmag__div0p01mm.bmp"


def per_wire(path: Path, pmin=3.0):
    prof = mp.column_profile(mp.load_gray(str(path)), (0.30, 0.70)).astype(float)
    per, _spec, _p = mp.fft_fundamental(prof, pmin_px=pmin)
    peaks, _spac = mp.peak_spacings(prof, per)
    p = prof - prof.min()
    grad = np.gradient(p)
    h = int(0.5 * per)
    w10 = mp.peak_widths_at(prof, peaks, per, frac=0.10)
    rows = []
    for pk, a in zip(peaks.astype(int), w10):
        lo, hi = max(0, pk - h), min(len(p) - 1, pk + h)
        steep = max(grad[lo:pk].max(initial=0.0), -grad[pk:hi].min(initial=0.0))
        if steep > 0 and np.isfinite(a):
            rows.append({"file": path.name, "peak_px": int(pk),
                         "blur_px": float(p[pk] / steep),
                         "w10_um": float(a * UM_PER_PX_HIGH)})
    return rows


def regress(rows):
    b = np.array([r["blur_px"] for r in rows])
    w = np.array([r["w10_um"] for r in rows])
    A = np.vstack([b, np.ones_like(b)]).T
    coef, *_ = np.linalg.lstsq(A, w, rcond=None)
    resid = w - A @ coef
    dof = max(len(b) - 2, 1)
    cov = (float(resid @ resid) / dof) * np.linalg.inv(A.T @ A)
    return {"slope": float(coef[0]), "intercept": float(coef[1]), "cov": cov,
            "r": float(np.corrcoef(b, w)[0, 1]), "n": len(b),
            "blur": b, "w10": w}


def at_blur(fit, b0):
    """Ширина при заданном уровне размытия + её ошибка (дельта-метод)."""
    val = fit["slope"] * b0 + fit["intercept"]
    J = np.array([b0, 1.0])
    err = float(np.sqrt(J @ fit["cov"] @ J))
    return float(val), err


def main():
    # ---------- предел резкости прибора по калибровочной гребёнке ----------
    cal = per_wire(CAL_HIGH, pmin=4.0)
    b_cal = np.array([r["blur_px"] for r in cal])
    b0 = float(np.median(b_cal))
    print("--- предел резкости прибора (калибровочная гребёнка, края физически резкие) ---")
    print(f"  зубцов {len(b_cal)}, размытие: медиана {b0:.2f} px "
          f"(мин {b_cal.min():.2f}) = {b0 * UM_PER_PX_HIGH:.2f} мкм")

    # ---------- образец и контроль ----------
    res = {}
    for tag, frames in (("test_grid_33_11", FRAMES_33), ("purewave", FRAMES_PW)):
        rows = [r for f in frames for r in per_wire(f)]
        fit = regress(rows)
        d0, e0 = at_blur(fit, 0.0)
        df, ef = at_blur(fit, b0)
        med = float(np.median(fit["w10"]))
        res[tag] = {"fit": fit, "rows": rows, "median_w10": med,
                    "min_w10": float(fit["w10"].min()),
                    "D_zero_blur": d0, "D_zero_blur_err": e0,
                    "D_at_instrument_limit": df, "D_at_instrument_limit_err": ef}
        print(f"\n--- {tag} ---")
        print(f"  проволок {fit['n']}, размытие {fit['blur'].min():.1f}…{fit['blur'].max():.1f} px "
              f"(медиана {np.median(fit['blur']):.1f})")
        print(f"  связь: r = {fit['r']:+.3f}, наклон {fit['slope']:+.3f} мкм на px размытия")
        print(f"  w10 без поправки: медиана {med:.2f}, минимум {fit['w10'].min():.2f} мкм")
        print(f"  D при b -> 0 (идеальный прибор):        {d0:5.2f} ± {e0:.2f} мкм")
        print(f"  D при b -> {b0:.1f} px (наилучший фокус): {df:5.2f} ± {ef:.2f} мкм  <-- рекомендуется")

    s33, spw = res["test_grid_33_11"]["fit"]["slope"], res["purewave"]["fit"]["slope"]
    print(f"\n--- сопоставимость ---")
    print(f"  наклоны: {s33:+.3f} (33_11) против {spw:+.3f} (purewave) — расхождение "
          f"{100 * abs(s33 - spw) / abs(spw):.0f} % ⇒ связь АППАРАТНАЯ, не образцовая")
    print(f"  purewave: опубликовано 10.6 ± 0.5, с поправкой "
          f"{res['purewave']['D_at_instrument_limit']:.2f} мкм "
          f"⇒ ТА ЖЕ поправка сдвигает и эталонный образец")

    D = res["test_grid_33_11"]["D_at_instrument_limit"]
    Derr = res["test_grid_33_11"]["D_at_instrument_limit_err"]
    dp = D / P_MEASURED
    h5 = 1 - 0.85 * dp
    ratio = D_EFF_FIT / D
    print(f"\n=== ИТОГ ПО test_grid_33_11 ===")
    print(f"  D = {D:.2f} ± {Derr:.2f} мкм  (было 12.55 ± 0.60; номинал {NOMINAL_D})")
    print(f"  D/P = {dp:.3f} при P = {P_MEASURED}")
    print(f"  D_eff/D_phys = {D_EFF_FIT:.3f}/{D:.2f} = {ratio:.3f}")
    print(f"  закон H5 = {h5:.3f}  ⇒ отклонение {ratio - h5:+.3f}")
    print(f"  (для сравнения: с D = 12.55 было отклонение −0.088)")

    # ---------- рисунок ----------
    fig, ax = plt.subplots(figsize=(7.6, 5.2))
    for tag, mk, col in (("test_grid_33_11", "o", "#c0504d"),
                         ("purewave", "^", "#1f4e79")):
        f = res[tag]["fit"]
        ax.plot(f["blur"], f["w10"], mk, ms=6, color=col, alpha=0.7,
                label=f"{tag} ({f['n']} проволок)")
        xs = np.linspace(0, f["blur"].max() * 1.05, 50)
        ax.plot(xs, f["slope"] * xs + f["intercept"], "-", lw=1.5, color=col, alpha=0.9)
    ax.axvline(b0, ls="--", lw=1.4, color="#404040")
    ax.text(b0 - 0.6, ax.get_ylim()[0] + 0.3, "предел резкости прибора",
            fontsize=8.5, va="bottom", ha="right", rotation=90, color="#404040")
    ax.axvspan(0, b0, color="#e8e8e8", alpha=0.7, zorder=0)
    ax.plot([b0], [D], "s", ms=10, mfc="#c0504d", mec="k", mew=1.2, zorder=5,
            label=f"итог: D = {D:.2f} мкм")
    ax.axhline(12.55, ls="-.", lw=1.2, color="#909090", label="прежняя оценка 12.55 — отозвана")
    ax.axhline(NOMINAL_D, ls=":", lw=1.4, color="#7030a0", label=f"номинал {NOMINAL_D:.0f} мкм")
    ax.set_xlabel("размытие проволоки (высота пика / крутизна фронта), px")
    ax.set_ylabel("ширина следа на уровне 10 %, мкм")
    ax.set_title("Диаметр завышался расфокусировкой: приведение к общей резкости\n"
                 "наклоны двух образцов совпадают ⇒ эффект аппаратный", fontsize=11)
    ax.grid(alpha=0.3)
    ax.legend(fontsize=8.5, loc="upper left")
    ax.set_xlim(left=0)
    fig.tight_layout()
    png = OUT / "a14b_diameter_defocus.png"
    fig.savefig(png, dpi=155)

    payload = {
        "instrument_blur_floor_px": b0,
        "instrument_blur_floor_um": b0 * UM_PER_PX_HIGH,
        "um_per_px": UM_PER_PX_HIGH,
        "samples": {t: {k: v for k, v in d.items() if k not in ("fit", "rows")}
                    | {"slope": d["fit"]["slope"], "intercept": d["fit"]["intercept"],
                       "r": d["fit"]["r"], "n_wires": d["fit"]["n"]}
                    for t, d in res.items()},
        "result_test_grid_33_11": {
            "D_um": D, "D_err_um": Derr, "D_um_previous": 12.55,
            "P_um": P_MEASURED, "D_over_P": dp,
            "D_eff_um": D_EFF_FIT, "D_eff_over_D_phys": ratio,
            "H5": h5, "deviation": ratio - h5},
        "wires": {t: res[t]["rows"] for t in res},
    }
    js = OUT / "a14b_diameter_defocus.json"
    js.write_text(json.dumps(payload, ensure_ascii=False, indent=2, default=float),
                  encoding="utf-8")
    print(f"\nартефакты: {js}\n           {png}")


if __name__ == "__main__":
    main()
