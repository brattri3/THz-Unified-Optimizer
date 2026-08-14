#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""A20 — сводка коэффициентов экстинкции по всем измеренным поляризаторам.

Экстинкция считается ПРЯМО ИЗ ДАННЫХ, без модели: для каждой частоты берётся
отношение максимума углового скана к минимуму,

    ER(nu) = |T(theta_max, nu)|^2 / |T(theta_min, nu)|^2,

то есть яркое положение против скрещенного. Величина относится ко ВСЕЙ схеме
(два реальных WGP), а не к одному элементу.

Сетка эксперимента берётся тем же билдером, что и все фиты зоны A
(`a3_eps_stability.build_masked`, полный угловой диапазон), поэтому маски воды и
динамического диапазона — те же, что в подгонках.

Запуск ИЗ КОРНЯ репозитория:
    .venv/Scripts/python.exe research/two_wgp/a20_extinction_summary.py
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

OUT = REPO / "research" / "results" / "a18_microscopy" / "a20_extinction.json"
PROBE_THZ = (0.5, 1.0, 1.5, 2.0)
BAND = (0.3, 2.0)
DATASETS = ["purewave", "specac", "test_grid_33_11", "test_grid_40_20",
            "att-11-16-356", "att-11-16-s1", "att-11-16-s2", "att-11-16-s3"]


def db(x):
    return 10.0 * np.log10(x)


def analyse(ds):
    ang, freqs, exp, valid, meta = a3.build_masked(ds, angles_limit=None)
    amp2 = np.abs(exp) ** 2                       # интенсивность, углы x частоты
    amp2 = np.where(valid, amp2, np.nan)

    with np.errstate(invalid="ignore"):
        i_max = np.nanargmax(np.nanmean(amp2, axis=1))
        i_min = np.nanargmin(np.nanmean(amp2, axis=1))
    er = amp2[i_max, :] / amp2[i_min, :]

    inband = (freqs >= BAND[0]) & (freqs <= BAND[1]) & np.isfinite(er)
    # зонд берём только если рядом с запрошенной частотой ЕСТЬ валидная точка;
    # иначе argmin молча возвращает край сетки и одно значение попадает в два столбца
    TOL_THZ = 0.05
    probes = {}
    for f0 in PROBE_THZ:
        near = np.abs(freqs - f0) <= TOL_THZ
        ok = near & np.isfinite(er)
        probes[f"{f0:.1f}"] = (float(db(np.nanmedian(er[ok]))) if ok.any() else None)

    j_best = int(np.nanargmax(np.where(inband, er, np.nan)))
    return {
        "dataset": ds, "n_angles": int(len(ang)), "n_freqs": int(len(freqs)),
        "angle_bright_deg": float(ang[i_max]), "angle_dark_deg": float(ang[i_min]),
        "ER_db_probes": probes,
        "ER_db_median_band": float(np.nanmedian(db(er[inband]))),
        "ER_db_max": float(db(er[j_best])), "f_at_max_thz": float(freqs[j_best]),
        "f_range_thz": [float(freqs[inband][0]), float(freqs[inband][-1])],
        "ER_db_at_band_edges": [float(db(er[inband][0])), float(db(er[inband][-1]))],
    }


def main():
    rows = []
    for ds in DATASETS:
        try:
            rows.append(analyse(ds))
        except Exception as exc:                                     # noqa: BLE001
            print(f"{ds}: пропущен ({type(exc).__name__}: {exc})")

    head = (f"{'образец':<17}{'углы':>6}{'ярк.':>7}{'тёмн.':>7}"
            + "".join(f"{f'{f:.1f} ТГц':>10}" for f in PROBE_THZ)
            + f"{'медиана':>10}{'макс':>9}{'на ν':>7}")
    print("\nКОЭФФИЦИЕНТ ЭКСТИНКЦИИ, дБ (схема из двух WGP, прямо из данных)\n")
    print(head)
    print("-" * len(head))
    for r in rows:
        cells = "".join(
            f"{r['ER_db_probes'][f'{f:.1f}']:>10.1f}"
            if r["ER_db_probes"][f"{f:.1f}"] is not None else f"{'—':>10}"
            for f in PROBE_THZ)
        print(f"{r['dataset']:<17}{r['n_angles']:>6}{r['angle_bright_deg']:>7.0f}"
              f"{r['angle_dark_deg']:>7.0f}{cells}"
              f"{r['ER_db_median_band']:>10.1f}{r['ER_db_max']:>9.1f}"
              f"{r['f_at_max_thz']:>7.2f}"
              f"   {r['f_range_thz'][0]:.2f}–{r['f_range_thz'][1]:.2f}")

    OUT.write_text(json.dumps({
        "definition": "ER(nu) = |T(theta_bright)|^2 / |T(theta_dark)|^2 по угловому скану; "
                      "относится ко всей схеме из двух WGP",
        "band_thz": BAND, "probes_thz": list(PROBE_THZ), "rows": rows,
    }, ensure_ascii=False, indent=2, default=float), encoding="utf-8")
    print(f"\nартефакт: {OUT}")


if __name__ == "__main__":
    main()
