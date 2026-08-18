"""Строит зашитую конфигурацию устройства для `attenuator_app/tools/service_calc.py`.

Однократный офлайн-шаг: гоняет полную Джонс-матричную модель `measured_curve.py`
(два идентичных WGP, схема S1) по дорожной серии измерений
`ATT-11-16-CA85_02721_s2` (`data_pool/two_wgp_attenuator/ATT-11-16-CA85_02721_s2/`
-- НЕ путать с плоской конвенцией `att-11-16-s2` в корне `data_pool`, это другой,
вырожденный на границе поиска датасет) и сохраняет результат (P, D, потери,
gamma, референсный спектр веса) в JSON рядом с этим файлом.

`service_calc.py` читает ТОЛЬКО этот JSON в рантайме -- data_pool/ приложением
не читается (тот же принцип, что у клиента v0.2: рантайм не зависит от научного
стека и лабораторных данных).

Перезапускать вручную при новой калибровке прибора:
    .venv\\Scripts\\python.exe -m attenuator_app.tools.calibration.build_service_calibration
"""
from __future__ import annotations

import json
import sys
from datetime import date
from pathlib import Path

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[2]
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))

from attenuator_app.tools.measured_curve import DATA_ROOT, analyze_road  # noqa: E402

DEVICE_ID = "ATT-11-16-CA85_02721"
DATASET_FOLDER = "ATT-11-16-CA85_02721_s2"
OUT_PATH = HERE / "att_11_16_ca85_02721.json"


def main() -> int:
    r = analyze_road(DATA_ROOT / DATASET_FOLDER)
    res = r["blanco_res"]
    if res.at_bound:
        print("фит упёрся в границу поиска -- калибровка невалидна, JSON не сохранён")
        return 1

    payload = {
        "device_id": DEVICE_ID,
        "calibration_dataset": DATASET_FOLDER,
        "calibration_convention": "road (NNN_aXX.txt/NNN_bg.txt)",
        "generated": date.today().isoformat(),
        "n_angles": int(len(r["angles"])),
        "band_thz": list(r["band"]),
        "P_um": float(r["P_fit"]),
        "D_um": float(r["D_fit"]),
        "loss_db_per_thz_gamma": float(r["loss_fit"]),
        "gamma": float(r["gamma_fit"]),
        "theta0_calibration_deg": float(r["theta0_fit"]),
        "at_bound": bool(res.at_bound),
        "fit_cost": float(res.cost),
        "freqs_ref_thz": [float(x) for x in r["freqs_ref"]],
        "power_ref": [float(x) for x in r["power_ref"]],
        "note": (
            "Полная Джонс-матричная модель схемы S1, оба WGP идентичны "
            "(P,D,потери,gamma -- параметры ОДНОГО WGP, не суммы двух). Без "
            "eta0/tau_leak (спорная феноменология зоны A, сюда не берётся). "
            "theta0_calibration_deg -- офсет ротатора В ЭТОЙ калибровочной "
            "серии, СПРАВОЧНО; приложение не использует его как рабочий ноль -- "
            "оператор обязан выполнить SET ZERO заново в своей сессии, т.к. "
            "прибор мог быть переустановлен в роторе между калибровкой и "
            "эксплуатацией. freqs_ref_thz/power_ref -- средний спектр мощности "
            "референсных (bg) импульсов калибровочной серии, используется как "
            "вес |E_ref(nu)|^2 при интегральном усреднении по ВСЕЙ записанной "
            "полосе (согласовано с теоремой Парсеваля, "
            "FINDINGS_measured_curve_2026-08-19.md п.4)."
        ),
    }
    OUT_PATH.write_text(json.dumps(payload, ensure_ascii=False, indent=2), encoding="utf-8")
    print(f"записано: {OUT_PATH}")
    print(f"P={payload['P_um']:.3f} D={payload['D_um']:.3f} мкм, "
          f"потери={payload['loss_db_per_thz_gamma']:.4f} дБ/ТГц^{payload['gamma']:.3f}, "
          f"theta0_calibration={payload['theta0_calibration_deg']:+.3f} град, "
          f"cost={payload['fit_cost']:.4g}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
