#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""A14b — перефит `test_grid_33_11` на ИЗМЕРЕННОЙ геометрии (A14) вместо номинальной.

ЗАЧЕМ. В A13 период был приколочен к номинальным `P = 33` мкм со слов владельца.
Микрофотографии (A14) дали `P = 31.75 ± 0.65` и `D = 12.55 ± 0.60`. Период входит
в модель Бланко как ФИКСИРОВАННЫЙ параметр, поэтому `D_eff` из старого фита к новой
геометрии не пересчитывается арифметически — фит надо повторить.

Сравниваются два прогона на одном и том же полном угловом диапазоне, отличающиеся
ТОЛЬКО геометрией.

Запуск ИЗ КОРНЯ репозитория:
    .venv/Scripts/python.exe research/two_wgp/a14_refit_measured_geometry.py
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

DATASET = "test_grid_33_11"
OUT = REPO / "research" / "results" / "two_wgp"

NOMINAL = {"P_um": 33.0, "D_phys_um": 11.0}
MEASURED = {"P_um": 31.75, "D_phys_um": 12.55}


def run(geom, tag):
    fit_lib.GEOMETRY[DATASET] = dict(geom)
    for m in ("a3_eps_stability",):
        sys.modules.pop(m, None)
    import a3_eps_stability as a3
    print(f"\n########## {tag}: P = {geom['P_um']}, D_phys = {geom['D_phys_um']} ##########")
    r = a3.analyse(DATASET, angles_limit=None, verbose=True)
    r.pop("starts", None)
    return r


def main():
    import a3_eps_stability as a3                                     # noqa: F401
    if DATASET not in a3.GROUP:
        a3.GROUP[DATASET] = "single_wgp"

    res = {"nominal": run(NOMINAL, "НОМИНАЛЬНАЯ геометрия (A13)"),
           "measured": run(MEASURED, "ИЗМЕРЕННАЯ геометрия (A14)")}

    print("\n=== СРАВНЕНИЕ ===")
    print(f"{'':<28}{'номинал':>12}{'измерено':>12}")
    for lab, key in (("P (пин), мкм", None), ("D_phys, мкм", None)):
        pass
    n, m = res["nominal"]["best"], res["measured"]["best"]
    print(f"{'P (приколот), мкм':<28}{NOMINAL['P_um']:>12.2f}{MEASURED['P_um']:>12.2f}")
    print(f"{'D_phys, мкм':<28}{NOMINAL['D_phys_um']:>12.2f}{MEASURED['D_phys_um']:>12.2f}")
    print(f"{'D_eff, мкм':<28}{n['D_um']:>12.3f}{m['D_um']:>12.3f}")
    print(f"{'D_eff/D_phys':<28}{n['D_over_Dphys']:>12.3f}{m['D_over_Dphys']:>12.3f}")
    dp_n = NOMINAL["D_phys_um"] / NOMINAL["P_um"]
    dp_m = MEASURED["D_phys_um"] / MEASURED["P_um"]
    print(f"{'D/P':<28}{dp_n:>12.3f}{dp_m:>12.3f}")
    print(f"{'закон H5 = 1-0.85*D/P':<28}{1 - 0.85 * dp_n:>12.3f}{1 - 0.85 * dp_m:>12.3f}")
    print(f"{'отклонение от H5':<28}"
          f"{n['D_over_Dphys'] - (1 - 0.85 * dp_n):>+12.3f}"
          f"{m['D_over_Dphys'] - (1 - 0.85 * dp_m):>+12.3f}")
    print(f"{'AIC':<28}{n['aic']:>12.1f}{m['aic']:>12.1f}")
    print(f"{'rmse, дБ':<28}{n['rmse_amp_db']:>12.3f}{m['rmse_amp_db']:>12.3f}")
    print(f"{'DS_rmse, дБ':<28}{n['deepshadow_rmse_db']:>12.3f}{m['deepshadow_rmse_db']:>12.3f}")
    print(f"{'Дарбин-Уотсон':<28}{n['durbin_watson']:>12.3f}{m['durbin_watson']:>12.3f}")

    payload = {"dataset": DATASET, "nominal_geometry": NOMINAL,
               "measured_geometry": MEASURED,
               "H5_nominal": 1 - 0.85 * dp_n, "H5_measured": 1 - 0.85 * dp_m,
               "runs": res}
    out = OUT / "a14_refit_measured_geometry.json"
    out.write_text(json.dumps(payload, ensure_ascii=False, indent=2, default=float),
                   encoding="utf-8")
    print(f"\nартефакт: {out}")


if __name__ == "__main__":
    main()
