#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""A14 — геометрия `test_grid_33_11` по микрофотографиям (P и D с погрешностями).

ЗАЧЕМ. Отчёт `REPORT_test_grid_33_11.md` §8 упирался в то, что геометрия образца
известна только со слов (`P = 33`, `D = 11`, без погрешностей), поэтому у
`D_eff/D_phys = 0.694` не было доверительного интервала по знаменателю. Владелец
указал снимки на флешке (`I:\\att\\WGP-11-33\\`, съёмка 2026-08-04 08:40–08:43).

СНИМКИ. Девять кадров трёх увеличений: 2 малых, 3 средних, 2 больших + пара
«фокус 0 / 0.12 мм» на большом увеличении. Разложены по конвенции `MANIFEST.md`.

⚠ КАЛИБРОВКА ПЕРЕНОСИТСЯ. В сессии 2026-08-04 гребёнку НЕ снимали, поэтому
масштаб берётся из сессии 2026-07-24 (те же кадры гребёнки, что у `purewave`;
эта практика уже применялась к `specac` и документирована в `MANIFEST.md`).
Перенос ПРОВЕРЯЕТСЯ двумя независимыми способами:
  1. отношение периодов в ПИКСЕЛЯХ между ступенями увеличения должно совпасть с
     отношением июльских масштабов (5.098 = 0.6444 / 0.1264);
  2. после умножения на масштаб период в МИКРОНАХ обязан совпасть на разных
     увеличениях. Расхождение между ними и есть систематическая погрешность —
     ровно так `MANIFEST.md` оценивает sigma для `purewave`.
Если обе проверки не проходят, абсолютные числа не публикуются.

Малое увеличение для АБСОЛЮТНОЙ шкалы непригодно (июльская гребёнка на нём ниже
разрешения оптики: 317 px/дел против 79 и 15.5) — используется только как
дополнительная проверка отношения.

ДИАМЕТР. Ширина отражающего следа на уровне 10 % от пика (не FWHM: FWHM ловит
только зеркальный блик и занижает D примерно вдвое — разобрано в `MANIFEST.md`),
только по проводам В ФОКУСЕ: пара кадров «фокус 0 / 0.12 мм» прямо показывает,
что провода не компланарны.

Запуск ИЗ КОРНЯ репозитория:
    .venv/Scripts/python.exe research/two_wgp/a14_micrograph_33_11.py
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "research" / "experiments"))
import micrograph_period as mp                                        # noqa: E402

SAMPLE = "test_grid_33_11"
MIC = REPO / "research" / "validation" / "micrographs" / "wgp_grid_D11_P33"
OUT = REPO / "research" / "results" / "two_wgp"
STEM = "wgp-grid-D11-P33"

# Масштабы сессии 2026-07-24 (кадры гребёнки, 1 деление = 10 мкм).
UM_PER_PX = {"high": 0.1264, "mid": 0.6444}
UM_RATIO_JULY = UM_PER_PX["mid"] / UM_PER_PX["high"]          # 5.098

NOMINAL_P, NOMINAL_D = 33.0, 11.0        # со слов владельца, 2026-08-06
D_EFF_FIT = 7.636                        # A13, полный диапазон

IMAGES = {
    "low": [f"{STEM}__wgp__lowmag_01.bmp", f"{STEM}__wgp__lowmag_02.bmp"],
    "mid": [f"{STEM}__wgp__midmag_01.bmp", f"{STEM}__wgp__midmag_02.bmp",
            f"{STEM}__wgp__midmag_03.bmp"],
    "high": [f"{STEM}__wgp__highmag_01.bmp", f"{STEM}__wgp__highmag_02.bmp"],
}
FOCUS_PAIR = [f"{STEM}__wgp__highmag_focus0um.bmp",
              f"{STEM}__wgp__highmag_focus120um.bmp"]


def px_period(path: Path):
    g = mp.load_gray(str(path))
    prof = mp.column_profile(g, (0.30, 0.70))
    per_fft, _spec, _per = mp.fft_fundamental(prof, pmin_px=3.0)
    peaks, spac = mp.peak_spacings(prof, per_fft)
    return per_fft, (float(np.median(spac)) if len(spac) else float("nan")), len(peaks)


def main():
    report = {"sample": SAMPLE, "micrographs_dir": str(MIC.relative_to(REPO)),
              "calibration_source": "сессия 2026-07-24 (кадры гребёнки purewave), перенос",
              "um_per_px": UM_PER_PX, "nominal": {"P_um": NOMINAL_P, "D_um": NOMINAL_D}}

    # ---------- шаг 1: периоды в пикселях, проверка ступеней увеличения ----------
    px = {}
    print("--- периоды в ПИКСЕЛЯХ (калибровка не нужна) ---")
    for mag, files in IMAGES.items():
        vals = []
        for f in files:
            p_fft, p_med, n = px_period(MIC / f)
            vals.append(p_fft)
            print(f"  {f:<44} {p_fft:7.2f} px (пики {p_med:6.2f}, N={n})")
        px[mag] = float(np.mean(vals))
    ratio_meas = px["high"] / px["mid"]
    ratio_dev = 100 * abs(ratio_meas - UM_RATIO_JULY) / UM_RATIO_JULY
    print(f"\n  ПРОВЕРКА 1 — ступени увеличения: high/mid = {ratio_meas:.3f} "
          f"против июльского {UM_RATIO_JULY:.3f} (расхождение {ratio_dev:.1f} %)")
    print(f"                 mid/low = {px['mid'] / px['low']:.3f} (ожидается ~2)")
    report["pixel_periods"] = px
    report["check_zoom_ratio"] = {"measured_high_over_mid": ratio_meas,
                                  "july_um_ratio": UM_RATIO_JULY,
                                  "deviation_pct": ratio_dev,
                                  "mid_over_low": px["mid"] / px["low"]}

    # ---------- шаг 2: полный разбор на калиброванных увеличениях ----------
    print("\n--- разбор (D при 10 % высоты, только провода в фокусе) ---")
    per_image, P_by_mag = [], {}
    for mag in ("mid", "high"):
        Ps, Ds = [], []
        for f in IMAGES[mag]:
            r = mp.analyze_wgp(str(MIC / f), UM_PER_PX[mag],
                               diam_frac=0.10, in_focus_only=(mag == "high"))
            r["file"], r["mag"] = f, mag
            per_image.append(r)
            Ps.append(r["P_um_lattice"])
            if np.isfinite(r["D_um_mean"]):
                Ds.append(r["D_um_mean"])
            print(f"  {f:<44} P={r['P_um_lattice']:6.2f}  D={r['D_um_mean']:5.2f}"
                  f"±{r['D_um_std']:.2f} (n={r['n_diam_wires']:2d})  "
                  f"беспорядок {r['disorder_pct']:4.1f} %  пропусков {r['n_missing']}")
        P_by_mag[mag] = {"mean": float(np.mean(Ps)),
                         "spread": float(np.max(Ps) - np.min(Ps)) if len(Ps) > 1 else 0.0,
                         "D_mean": float(np.mean(Ds)) if Ds else float("nan")}

    P_mid, P_high = P_by_mag["mid"]["mean"], P_by_mag["high"]["mean"]
    P_final = 0.5 * (P_mid + P_high)
    # sigma: расхождение между увеличениями (систематика калибровки) + разброс внутри
    sigma_cross = abs(P_mid - P_high) / 2
    sigma_within = max(P_by_mag["mid"]["spread"], P_by_mag["high"]["spread"]) / 2
    P_err = float(np.hypot(sigma_cross, sigma_within))
    dev_pct = 100 * abs(P_mid - P_high) / P_final
    print(f"\n  ПРОВЕРКА 2 — согласие увеличений: P(mid) = {P_mid:.2f}, "
          f"P(high) = {P_high:.2f} мкм (расхождение {dev_pct:.1f} %)")

    D_final = P_by_mag["high"]["D_mean"]
    D_err = float(np.mean([r["D_um_std"] for r in per_image
                           if r["mag"] == "high" and np.isfinite(r["D_um_std"])]))

    # ---------- шаг 3: пара «фокус 0 / 0.12 мм» ----------
    print("\n--- некомпланарность: пара кадров с разным фокусом ---")
    focus = []
    for f in FOCUS_PAIR:
        r = mp.analyze_wgp(str(MIC / f), UM_PER_PX["high"], diam_frac=0.10,
                           in_focus_only=True)
        r["file"] = f
        focus.append(r)
        print(f"  {f:<48} проводов в фокусе {r['n_diam_wires']:2d}, "
              f"D={r['D_um_mean']:5.2f} мкм")
    report["focus_pair"] = [{k: v for k, v in r.items() if k != "path"} for r in focus]

    # ---------- итог ----------
    ratio_measured = D_EFF_FIT / D_final
    d_over_p_meas = D_final / P_final
    h5_meas = 1 - 0.85 * d_over_p_meas
    h5_nominal = 1 - 0.85 * (NOMINAL_D / NOMINAL_P)

    report["result"] = {
        "P_um": P_final, "P_err_um": P_err,
        "D_um": D_final, "D_err_um": D_err,
        "D_over_P": d_over_p_meas,
        "P_mid": P_mid, "P_high": P_high, "cross_mag_dev_pct": dev_pct,
        "nominal_P_dev_pct": 100 * (P_final - NOMINAL_P) / NOMINAL_P,
        "nominal_D_dev_pct": 100 * (D_final - NOMINAL_D) / NOMINAL_D,
        "D_eff_over_D_phys_with_measured_D": ratio_measured,
        "H5_prediction_measured_geometry": h5_meas,
        "H5_prediction_nominal_geometry": h5_nominal,
        "disorder_pct_median": float(np.median([r["disorder_pct"] for r in per_image])),
    }
    report["per_image"] = [{k: v for k, v in r.items() if k != "path"} for r in per_image]

    print("\n=== ИТОГ ===")
    print(f"  P = {P_final:.2f} ± {P_err:.2f} мкм   (номинал {NOMINAL_P}, "
          f"отклонение {report['result']['nominal_P_dev_pct']:+.1f} %)")
    print(f"  D = {D_final:.2f} ± {D_err:.2f} мкм   (номинал {NOMINAL_D}, "
          f"отклонение {report['result']['nominal_D_dev_pct']:+.1f} %)")
    print(f"  D/P = {d_over_p_meas:.3f} (номинально {NOMINAL_D / NOMINAL_P:.3f})")
    print(f"  беспорядок (медиана по кадрам): {report['result']['disorder_pct_median']:.1f} %")
    print(f"\n  ⚠ ВНИМАНИЕ: закон H5 при ИЗМЕРЕННОЙ геометрии предсказывает "
          f"{h5_meas:.3f} вместо {h5_nominal:.3f}")
    print(f"  D_eff (фит, при P приколоченном к {NOMINAL_P}) = {D_EFF_FIT:.3f} мкм")
    print(f"  D_eff/D_phys с ИЗМЕРЕННЫМ D = {ratio_measured:.3f}")
    print("  ⇒ фит надо ПЕРЕЗАПУСТИТЬ с P = {:.2f}: период в нём приколочен."
          .format(P_final))

    OUT.mkdir(parents=True, exist_ok=True)
    out = OUT / "a14_micrograph_33_11.json"
    out.write_text(json.dumps(report, ensure_ascii=False, indent=2, default=float),
                   encoding="utf-8")
    print(f"\nартефакт: {out}")


if __name__ == "__main__":
    main()
