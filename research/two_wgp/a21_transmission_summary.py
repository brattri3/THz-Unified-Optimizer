#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""A21 — максимальное пропускание и его спад по частоте, по всем образцам.

Для каждого набора берётся ЯРКОЕ положение углового скана (оси пропускания
сонаправлены) и считаются три величины:

  T_time   интегральное пропускание во временной области:
           ∫E_sig² dt / ∫E_bg² dt — усреднение по полосе не требуется вовсе
           (по Парсевалю это тот же интеграл по спектру, но без выбора границ)
  T_band   медиана |t(ν)|² по полосе 0.3–1.5 ТГц — спектральный аналог для сверки
  spad     наклон 10·lg|t(ν)|² по частоте, дБ/ТГц: насколько пропускание
           падает с ростом частоты (то, что владелец назвал дисперсией)

Запуск ИЗ КОРНЯ репозитория:
    .venv/Scripts/python.exe research/two_wgp/a21_transmission_summary.py
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[2]
HERE = REPO / "research" / "two_wgp"
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(REPO / "research" / "experiments"))
sys.path.insert(0, str(REPO))

import a3_eps_stability as a3                                        # noqa: E402

OUT = REPO / "research" / "results" / "a18_microscopy" / "a21_transmission.json"
BAND = (0.3, 1.5)
DATASETS = ["purewave", "specac", "test_grid_33_11", "test_grid_40_20",
            "att-11-16-356", "att-11-16-s1", "att-11-16-s2", "att-11-16-s3"]


def per_angle(data, angle):
    t_s, E_s, t_b, E_b = data[angle]
    f, s_s, s_b, tr = a3.get_transmission_spectra(t_s, E_s, t_b, E_b)
    T = np.abs(tr) ** 2
    # временная область: отношение энергий импульсов, без выбора полосы
    T_time = float(np.trapezoid(np.asarray(E_s) ** 2, t_s)
                   / np.trapezoid(np.asarray(E_b) ** 2, t_b))
    return np.asarray(f), T, T_time


def analyse(ds):
    data = a3.DataManager(a3.config.DATA_DIR).get_data_for_dataset(ds)
    angles = a3.exp_builder.select_angles(data, None)

    best = None
    for a in angles:
        f, T, T_time = per_angle(data, a)
        if best is None or T_time > best[3]:
            best = (a, f, T, T_time)
    a_bright, f, T, T_time = best

    m = (f >= BAND[0]) & (f <= BAND[1]) & np.isfinite(T) & (T > 0)
    fb, Tb = f[m], T[m]
    db = 10.0 * np.log10(Tb)
    slope, icept = np.polyfit(fb, db, 1)
    resid = db - (icept + slope * fb)
    # стандартная ошибка наклона: наклоны здесь малы, без неё таблица нечитаема
    s_res = float(np.sqrt(np.sum(resid ** 2) / (len(fb) - 2)))
    slope_err = float(s_res / np.sqrt(np.sum((fb - fb.mean()) ** 2)))

    return {
        "dataset": ds, "angle_bright_deg": float(a_bright),
        "T_time": T_time, "T_time_pct": 100.0 * T_time,
        "T_band_median": float(np.median(Tb)), "T_band_median_pct": float(100 * np.median(Tb)),
        "T_at_0p5": float(np.interp(0.5, fb, Tb)), "T_at_1p0": float(np.interp(1.0, fb, Tb)),
        "slope_db_per_thz": float(slope), "slope_err_db_per_thz": slope_err,
        "intercept_db": float(icept),
        "resid_rms_db": float(np.sqrt(np.mean(resid ** 2))),
        "std_db_in_band": float(np.std(db)),
        "band_thz": list(BAND), "n_freqs": int(m.sum()),
    }


def main():
    rows = []
    for ds in DATASETS:
        try:
            rows.append(analyse(ds))
        except Exception as exc:                                     # noqa: BLE001
            print(f"{ds}: пропущен ({type(exc).__name__}: {exc})")

    print("\nПРОПУСКАНИЕ В ЯРКОМ ПОЛОЖЕНИИ\n")
    print(f"{'образец':<17}{'угол':>6}{'T врем.,%':>11}{'T полоса,%':>12}"
          f"{'T(0.5)':>9}{'T(1.0)':>9}{'спад, дБ/ТГц':>18}{'СКО,дБ':>9}")
    for r in rows:
        sl = f"{r['slope_db_per_thz']:+.2f} ± {r['slope_err_db_per_thz']:.2f}"
        print(f"{r['dataset']:<17}{r['angle_bright_deg']:>6.0f}{r['T_time_pct']:>11.1f}"
              f"{r['T_band_median_pct']:>12.1f}{r['T_at_0p5']:>9.3f}{r['T_at_1p0']:>9.3f}"
              f"{sl:>18}{r['std_db_in_band']:>9.2f}")

    OUT.write_text(json.dumps({
        "definitions": {
            "T_time": "интеграл E^2 по времени: сигнал/фон, без выбора полосы",
            "T_band_median": "медиана |t(nu)|^2 по полосе 0.3-1.5 ТГц",
            "slope_db_per_thz": "наклон 10lg|t|^2 по частоте в той же полосе",
        },
        "band_thz": list(BAND), "rows": rows,
    }, ensure_ascii=False, indent=2, default=float), encoding="utf-8")
    print(f"\nартефакт: {OUT}")


if __name__ == "__main__":
    main()
