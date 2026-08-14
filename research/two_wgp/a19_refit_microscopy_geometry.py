#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""A19 — перефит M5 на геометрии из микрофото A18 (2026-08-10).

ЗАЧЕМ. Период входит в модель Бланко приколоченным, поэтому при его сдвиге `D_eff`
нельзя пересчитать арифметически — фит надо повторить. Микрофото A18 сдвинули период
у `purewave` и `test_grid_40_20` на −4.2 %; у `specac` и `test_grid_33_11` сдвиг
в пределах погрешности, и эти два образца работают КОНТРОЛЕМ: их `D_eff` обязан
воспроизвести прежнее значение из `a16_h5_recompute.json:D_eff_used`.

Методика повторяет A14b (`a14_refit_measured_geometry.py`): подмена `fit_lib.GEOMETRY`
и прогон `a3_eps_stability.analyse` на ПОЛНОМ угловом диапазоне (`angles_limit=None`).
Меняется только геометрия — код модели, данные и стартовая схема те же.

Запуск ИЗ КОРНЯ репозитория:
    .venv/Scripts/python.exe research/two_wgp/a19_refit_microscopy_geometry.py
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
HERE = REPO / "research" / "two_wgp"
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(REPO / "research" / "experiments"))
sys.path.insert(0, str(REPO))

import fit_lib                                                        # noqa: E402

A19 = REPO / "research" / "results" / "a18_microscopy" / "a19_h5_dp.json"
A16 = REPO / "research" / "results" / "two_wgp" / "a16_h5_recompute.json"
OUT = REPO / "research" / "results" / "a18_microscopy" / "a19_refit.json"

# Контроль: период сдвинулся в пределах погрешности, D_eff обязан воспроизвестись
CONTROL = {"specac", "test_grid_33_11"}
ORDER = ["specac", "test_grid_33_11", "purewave", "test_grid_40_20"]


def run(dataset, geom):
    fit_lib.GEOMETRY[dataset] = dict(geom)
    sys.modules.pop("a3_eps_stability", None)
    import a3_eps_stability as a3
    if dataset not in a3.GROUP:
        a3.GROUP[dataset] = "single_wgp"
    tag = "КОНТРОЛЬ" if dataset in CONTROL else "ПЕРЕФИТ"
    print(f"\n{'#' * 78}\n### {tag} {dataset}: P = {geom['P_um']:.3f}, "
          f"D_phys = {geom['D_phys_um']:.3f}\n{'#' * 78}", flush=True)
    r = a3.analyse(dataset, angles_limit=None, verbose=True)
    r.pop("starts", None)
    return r


def main():
    a19 = json.loads(A19.read_text(encoding="utf-8"))
    a16 = json.loads(A16.read_text(encoding="utf-8"))
    new_geo = {p["sample"]: p for p in a19["points"]}
    old_deff = a16["D_eff_used"]

    runs, rows = {}, []
    for ds in ORDER:
        g = new_geo[ds]
        geom = {"P_um": g["P_um"], "D_phys_um": g["D_um"]}
        r = run(ds, geom)
        runs[ds] = r
        b = r["best"]
        rows.append({
            "sample": ds,
            "control": ds in CONTROL,
            "P_um": geom["P_um"], "D_phys_um": geom["D_phys_um"],
            "D_over_P": geom["D_phys_um"] / geom["P_um"],
            "D_eff_um": b["D_um"],
            "D_eff_old_um": old_deff[ds],
            "dD_eff_um": b["D_um"] - old_deff[ds],
            "ratio_new": b["D_um"] / geom["D_phys_um"],
            "ratio_prefit": old_deff[ds] / geom["D_phys_um"],
            "aic": b["aic"], "rmse_amp_db": b["rmse_amp_db"],
            "deepshadow_rmse_db": b.get("deepshadow_rmse_db"),
            "durbin_watson": b.get("durbin_watson"),
        })

    # закон H5 на перефиченных значениях: y = 1 - k*(D/P)
    import numpy as np
    x = np.array([r["D_over_P"] for r in rows])
    y = np.array([r["ratio_new"] for r in rows])
    k = float(np.sum(x * (1 - y)) / np.sum(x * x))
    res = y - (1 - k * x)
    rms = float(np.sqrt(np.mean(res ** 2)))
    k_err = float(np.sqrt(np.sum(res ** 2) / (len(x) - 1)) / np.sqrt(np.sum(x * x)))

    print(f"\n{'=' * 78}\n=== ИТОГ: D_eff до и после перефита ===")
    print(f"{'образец':<18}{'D_eff до':>10}{'после':>9}{'сдвиг':>9}"
          f"{'ratio до':>10}{'после':>9}{'AIC':>11}")
    for r in rows:
        mark = " К" if r["control"] else ""
        print(f"{r['sample']:<18}{r['D_eff_old_um']:>10.3f}{r['D_eff_um']:>9.3f}"
              f"{r['dD_eff_um']:>+9.3f}{r['ratio_prefit']:>10.3f}"
              f"{r['ratio_new']:>9.3f}{r['aic']:>11.1f}{mark}")
    print(f"\nзакон H5 после перефита: k = {k:.4f} ± {k_err:.4f}, rms = {rms:.4f}")
    print("(до перефита, только геометрия: k = 0.8825 ± 0.0236, rms = 0.0176)")

    OUT.write_text(json.dumps({
        "source": "перефит M5 на геометрии микрофото A18; методика A14b",
        "control_samples": sorted(CONTROL),
        "rows": rows,
        "h5_fit_after_refit": {"form": "1 - k*(D/P)", "k": k, "k_err": k_err, "rms": rms},
        "runs": runs,
    }, ensure_ascii=False, indent=2, default=float), encoding="utf-8")
    print(f"\nартефакт: {OUT}")


if __name__ == "__main__":
    main()
