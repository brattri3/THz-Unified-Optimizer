#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""A16 — единая геометрия всех образцов: период и диаметр с поправкой на расфокусировку.

ПОВОД. A14b показал, что диаметр систематически завышается расфокусировкой, причём
на разную величину у разных образцов (наклон «размытие -> ширина» аппаратный, а
разброс размытия — свойство намотки). Значит опубликованные `D` несопоставимы
между строками, а на этих же `D/P` подбирался эмпирический закон H5. Владелец
поручил пересчитать всё единообразно и поправить закон.

ДВЕ ПРОБЛЕМЫ, ОБНАРУЖЕННЫЕ ПРИ ПОДГОТОВКЕ, И КАК ОНИ РЕШЕНЫ
------------------------------------------------------------
1. МАСШТАБ «ПЛАВАЕТ» МЕЖДУ КАДРАМИ. Отношение периодов в пикселях high/mid должно
   равняться отношению калибровочных масштабов 0.6444/0.1264 = 5.098. Фактически:
   test_grid_33_11 5.02, test_grid_40_20 5.29, purewave 5.39, specac **4.39**.
   То есть зум между съёмками менялся, и применять общий масштаб `high` ко всем
   кадрам большого увеличения НЕЛЬЗЯ.
   РЕШЕНИЕ: каждый кадр большого увеличения калибруется по СОБСТВЕННОМУ периоду
   образца: `мкм/px = P_мкм / P_px(кадра)`, где `P_мкм` взят со среднего
   увеличения. Зум-дрейф сокращается тождественно. Побочно: отношение `D/P`
   вообще не требует калибровки — это отношение пикселей внутри одного кадра.
2. ПЕРИОД БРАЛСЯ ПО-РАЗНОМУ. В `MANIFEST.md` у `purewave` период — смесь среднего
   и большого увеличения (24.7 и 25.8 -> 25.5), у `specac` только среднее.
   РЕШЕНИЕ: у всех образцов период берётся со СРЕДНЕГО увеличения (там гребёнка
   разрешена, периодов в кадре много, а проволоки ещё не сливаются).

ЧЕГО СДЕЛАТЬ НЕЛЬЗЯ. У `att-11-16` (каталог `wgp_grid_D11_P16`) есть только кадры
среднего увеличения, и проволоки на них почти сомкнуты (`D/P ~ 0.69`): выделение
отдельных проволок неустойчиво, FFT цепляется за кластеры яркости (периоды 62,
103, 146 px вместо ожидаемых ~25). Поправка на расфокусировку требует большого
увеличения, поэтому четыре аттенюаторные точки остаются с ПАСПОРТНОЙ геометрией
и в пересчёте закона участвуют отдельной пометкой.

Запуск ИЗ КОРНЯ репозитория:
    .venv/Scripts/python.exe research/two_wgp/a16_geometry_all.py
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "research" / "experiments"))
import micrograph_period as mp                                        # noqa: E402

MIC = REPO / "research" / "validation" / "micrographs"
OUT = REPO / "research" / "results" / "two_wgp"

UM_PER_PX_MID = 0.6444          # калибровка среднего увеличения (гребёнка 2026-07-24)
UM_PER_PX_HIGH = 0.1264         # только для перевода предела резкости в микроны
BLUR_FLOOR_PX_AT_CAL = 8.69     # медиана размытия на калибровочной гребёнке
BLUR_FLOOR_UM = BLUR_FLOOR_PX_AT_CAL * UM_PER_PX_HIGH      # ~1.10 мкм в пространстве объекта

SAMPLES = {
    "purewave": {
        "dir": "purewave",
        "mid": ["purewave__wgp__midmag.bmp"],
        "high": [f"purewave__wgp__highmag_0{i}.bmp" for i in (1, 2, 3)],
        "geometry_now": {"P_um": 25.5, "D_phys_um": 10.6},
    },
    "specac": {
        "dir": "wgp_specac",
        "mid": ["specac__wgp__midmag.bmp"],
        # highmag_04 выпадает по масштабу (198 px против 169.7 у трёх остальных)
        # — снят при другом зуме; на самокалибровке по периоду это не мешает.
        "high": [f"specac__wgp__highmag_0{i}.bmp" for i in (1, 2, 3, 4)],
        "geometry_now": {"P_um": 24.9, "D_phys_um": 14.0},
    },
    "test_grid_40_20": {
        "dir": "wgp_grid_D20_P40",
        "mid": ["wgp_grid_D20_P40__wgp__midmag.bmp"],
        "high": ["wgp_grid_D20_P40__wgp__highmag_focus0um.bmp",
                 "wgp_grid_D20_P40__wgp__highmag_focus200um.bmp"],
        "geometry_now": {"P_um": 38.8, "D_phys_um": 18.0},
    },
    "test_grid_33_11": {
        "dir": "wgp_grid_D11_P33",
        "mid": [f"wgp-grid-D11-P33__wgp__midmag_0{i}.bmp" for i in (1, 2, 3)],
        "high": ["wgp-grid-D11-P33__wgp__highmag_01.bmp",
                 "wgp-grid-D11-P33__wgp__highmag_02.bmp",
                 "wgp-grid-D11-P33__wgp__highmag_focus0um.bmp",
                 "wgp-grid-D11-P33__wgp__highmag_focus120um.bmp"],
        "geometry_now": {"P_um": 31.75, "D_phys_um": 10.91},
    },
}

# Аттенюатор: обрабатывать нечем, но зафиксировать причину в артефакте.
NOT_MEASURABLE = {
    "att-11-16-*": {
        "dir": "wgp_grid_D11_P16",
        "reason": ("только кадры среднего увеличения; проволоки почти сомкнуты (D/P ~ 0.69), "
                   "выделение отдельных проволок неустойчиво — FFT даёт 62/103/146 px вместо "
                   "ожидаемых ~25. Поправка на расфокусировку требует большого увеличения."),
        "geometry_kept": {"P_um": "15.5 (-356) / 16.0 (-s1..-s3)", "D_phys_um": 11.0,
                          "source": "паспорт прибора, ПРИОР, не измерение"},
    }
}


def profile(path: Path):
    return mp.column_profile(mp.load_gray(str(path)), (0.30, 0.70)).astype(float)


def period_px(prof, pmin=3.0):
    per, _spec, _p = mp.fft_fundamental(prof, pmin_px=pmin)
    return per


def wires_of(path: Path, um_per_px: float):
    """Проволоки кадра: ширина следа на 10 % и размытие, обе в МИКРОНАХ."""
    prof = profile(path)
    per = period_px(prof)
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
            rows.append({"file": path.name, "w10_um": float(a * um_per_px),
                         "blur_um": float((p[pk] / steep) * um_per_px),
                         "w10_over_P": float(a / per)})
    return rows, per


def main():
    print(f"предел резкости прибора: {BLUR_FLOOR_PX_AT_CAL:.2f} px при "
          f"{UM_PER_PX_HIGH} мкм/px = **{BLUR_FLOOR_UM:.2f} мкм** в пространстве объекта\n")

    report = {"blur_floor_um": BLUR_FLOOR_UM, "um_per_px_mid": UM_PER_PX_MID,
              "samples": {}, "not_measurable": NOT_MEASURABLE}

    for name, cfg in SAMPLES.items():
        d = MIC / cfg["dir"]
        # --- период со среднего увеличения
        Ps = []
        for f in cfg["mid"]:
            prof = profile(d / f)
            per = period_px(prof)
            peaks, _ = mp.peak_spacings(prof, per)
            P_lat, _rms, _n, _miss = mp.lattice_fit(peaks, per)
            Ps.append(P_lat * UM_PER_PX_MID)
        P_um = float(np.mean(Ps))
        P_err = float(np.std(Ps, ddof=1)) if len(Ps) > 1 else 0.02 * P_um

        # --- проволоки со всех кадров большого увеличения, каждый самокалиброван
        rows, frames = [], []
        for f in cfg["high"]:
            prof = profile(d / f)
            per_px = period_px(prof)
            scale = P_um / per_px                       # мкм/px ДЛЯ ЭТОГО КАДРА
            r, _ = wires_of(d / f, scale)
            rows += r
            frames.append({"file": f, "P_px": float(per_px), "um_per_px": float(scale)})

        b = np.array([x["blur_um"] for x in rows])
        w = np.array([x["w10_um"] for x in rows])
        A = np.vstack([b, np.ones_like(b)]).T
        coef, *_ = np.linalg.lstsq(A, w, rcond=None)
        resid = w - A @ coef
        cov = (float(resid @ resid) / max(len(b) - 2, 1)) * np.linalg.inv(A.T @ A)
        slope, intercept = float(coef[0]), float(coef[1])
        J = np.array([BLUR_FLOOR_UM, 1.0])
        D_um = float(slope * BLUR_FLOOR_UM + intercept)
        D_err = float(np.sqrt(J @ cov @ J))
        r_corr = float(np.corrcoef(b, w)[0, 1])

        old = cfg["geometry_now"]
        rec = {"P_um": P_um, "P_err_um": P_err,
               "D_um": D_um, "D_err_um": D_err,
               "D_over_P": D_um / P_um,
               "D_um_uncorrected_median": float(np.median(w)),
               "slope_um_per_um_blur": slope, "corr_w10_blur": r_corr,
               "n_wires": len(rows), "n_high_frames": len(cfg["high"]),
               "blur_range_um": [float(b.min()), float(b.max())],
               "geometry_previous": old,
               "P_change_pct": 100 * (P_um - old["P_um"]) / old["P_um"],
               "D_change_pct": 100 * (D_um - old["D_phys_um"]) / old["D_phys_um"],
               "frames": frames}
        report["samples"][name] = rec

        print(f"--- {name} ---")
        print(f"  период (среднее увеличение): P = {P_um:6.2f} ± {P_err:.2f} мкм "
              f"(было {old['P_um']}, {rec['P_change_pct']:+.1f} %)")
        print(f"  масштаб кадров high: " +
              ", ".join(f"{fr['um_per_px']:.4f}" for fr in frames) + " мкм/px "
              f"(общая калибровка дала бы {UM_PER_PX_HIGH})")
        print(f"  проволок {len(rows)}, размытие {b.min():.2f}…{b.max():.2f} мкм, "
              f"r = {r_corr:+.3f}, наклон {slope:+.3f}")
        print(f"  диаметр: без поправки медиана {np.median(w):5.2f} -> "
              f"**{D_um:5.2f} ± {D_err:.2f} мкм** (было {old['D_phys_um']}, "
              f"{rec['D_change_pct']:+.1f} %)")
        print(f"  D/P = {D_um / P_um:.3f} (было {old['D_phys_um'] / old['P_um']:.3f})\n")

    print("--- НЕ измеряется ---")
    for k, v in NOT_MEASURABLE.items():
        print(f"  {k}: {v['reason']}")

    out = OUT / "a16_geometry_all.json"
    out.write_text(json.dumps(report, ensure_ascii=False, indent=2, default=float),
                   encoding="utf-8")
    print(f"\nартефакт: {out}")


if __name__ == "__main__":
    main()
