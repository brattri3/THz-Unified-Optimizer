"""Библиотека + CLI: калькулятор аттенюатора для обслуживания в THz-TDS
спектрометре (задача C9, санкция владельца 2026-08-19, начата после C8).

Один сценарий эксплуатации, двунаправленный:
  (а) дБ -> угол:  оператор вводит желаемое затухание, получает угол(ы),
      на которые поставить WGP1;
  (б) угол -> дБ:  оператор вводит текущий угол WGP1 (по шкале ротатора),
      получает предсказанное затухание.
Оба направления -- ОТНОСИТЕЛЬНО точки SET ZERO, которую оператор задаёт сам
в начале сессии вручную (владелец 2026-08-19: алгоритм автоматической
калибровки офсета -- следующая версия, не входит в этот минимальный релиз).
Метрика затухания -- интегральная (широкополосная, по умолчанию) и/или на
конкретной частоте, если она указана. Только точечное предсказание, БЕЗ
доверительного интервала (санкция владельца 2026-08-19: минимальный масштаб
задачи).

Физика -- полная Джонс-матричная модель схемы S1 (два ИДЕНТИЧНЫХ WGP) из
`attenuator_app/tools/measured_curve.py`; вывод формулы и обоснование --
`FINDINGS_measured_curve_2026-08-19.md`, п.1 (I_perp=|t_perp|^4, не |t_perp|^2
одного WGP). Параметры устройства ЗАШИТЫ в
`calibration/att_11_16_ca85_02721.json` (как получено --
`calibration/build_service_calibration.py`); рантайм `data_pool/` не читает и
зависит только от numpy + `attenuator_app.core.blanco` (тот же принцип
изоляции от научного стека, что у клиента v0.2).

НЕ часть клиентского `attenuator_app.gui`/`cli` (v0.2/v0.3, отдельный трек) --
самостоятельный инструмент поверх модели C8.

Запуск из корня репозитория:
    .venv\\Scripts\\python.exe -m attenuator_app.tools.service_calc --zero -0.4 --to-db 12
    .venv\\Scripts\\python.exe -m attenuator_app.tools.service_calc --zero -0.4 --from-angle 35 --freq 0.8
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))

from attenuator_app.core.blanco import dressed_t          # noqa: E402

DEFAULT_CALIBRATION_PATH = HERE / "calibration" / "att_11_16_ca85_02721.json"


class Calibration:
    """Зашитая конфигурация устройства -- см. `calibration/*.json`."""

    def __init__(self, data: dict):
        self.device_id = data["device_id"]
        self.dataset = data["calibration_dataset"]
        self.generated = data.get("generated", "?")
        self.P_um = float(data["P_um"])
        self.D_um = float(data["D_um"])
        self.loss_db = float(data["loss_db_per_thz_gamma"])
        self.gamma = float(data["gamma"])
        self.band_thz = tuple(data["band_thz"])
        self.theta0_calibration_deg = float(data["theta0_calibration_deg"])
        self.at_bound = bool(data["at_bound"])
        self.freqs_ref = np.array(data["freqs_ref_thz"], dtype=float)
        self.power_ref = np.array(data["power_ref"], dtype=float)


def load_calibration(path: Path = DEFAULT_CALIBRATION_PATH) -> Calibration:
    with open(path, encoding="utf-8") as f:
        return Calibration(json.load(f))


def band_warning(freq_thz: float | None, cal: Calibration) -> str | None:
    if freq_thz is None:
        return None
    lo, hi = cal.band_thz
    if lo <= freq_thz <= hi:
        return None
    return (f"[!] {freq_thz:.3f} ТГц вне откалиброванной полосы "
            f"{lo:.2f}-{hi:.2f} ТГц -- значение экстраполировано, не измерено")


def _transmission_ratio(theta_deg: np.ndarray, theta0_deg: float, cal: Calibration,
                        freq_thz: float | None) -> np.ndarray:
    """T(theta)/T(theta0) для схемы S1 (два идентичных WGP), см.
    `FINDINGS_measured_curve_2026-08-19.md` п.1 и
    `measured_curve.blanco_angular_curve` (тот же вывод, продублирован здесь
    в минимальном виде, чтобы не тянуть импорт `track_viewer` в рантайм этого
    инструмента).

    `theta_deg` -- массив любой формы, градусы. `freq_thz=None` -> интегрально
    по калибровочному весу `|E_ref(nu)|^2` (широкая записанная полоса,
    теорема Парсеваля -- FINDINGS п.4); иначе -- одна частота (плоское
    "среднее" по одному элементу, т.е. само значение в этой точке).
    """
    if freq_thz is None:
        freqs, weight = cal.freqs_ref, cal.power_ref
    else:
        freqs, weight = np.array([float(freq_thz)]), None

    tp, ta, _ = dressed_t(freqs, cal.P_um, cal.D_um, loss_factor=cal.loss_db, gamma=cal.gamma)
    d = np.deg2rad(np.asarray(theta_deg, dtype=float) - theta0_deg)
    c2, s2 = np.cos(d) ** 2, np.sin(d) ** 2
    E1 = tp[None, :] * c2[..., None] + ta[None, :] * s2[..., None]

    def wavg(x):
        return np.mean(x, axis=-1) if weight is None else np.average(x, axis=-1, weights=weight)

    U, U0 = wavg(np.abs(E1) ** 2), wavg(np.abs(tp) ** 2)
    return U / U0


def attenuation_db(theta_deg: float, theta0_deg: float, cal: Calibration,
                    freq_thz: float | None = None) -> float:
    """Прямая задача: предсказанное затухание (дБ) на угле `theta_deg`."""
    T = float(_transmission_ratio(np.array([theta_deg]), theta0_deg, cal, freq_thz)[0])
    return -10.0 * np.log10(max(T, 1e-300))


def angle_for_db(target_db: float, theta0_deg: float, cal: Calibration,
                 freq_thz: float | None = None, n: int = 901) -> dict:
    """Обратная задача: угол(ы) WGP1 для желаемого затухания.

    Кривая T(delta), delta=theta-theta0 in [0,90], монотонно убывает от 1
    (delta=0, положение SET ZERO) до минимума в скрещенном положении
    (delta=90). Симметрична по знаку delta (cos^2/sin^2 -- чётные функции),
    поэтому решения ДВА: theta0+delta и theta0-delta, физически равнозначны --
    какое ближе к текущему положению ротатора, решает оператор (моторизации
    нет, `attenuator_app` C4_motor -- todo).
    """
    if target_db < 0:
        raise ValueError("затухание не может быть отрицательным")
    delta = np.linspace(0.0, 90.0, n)
    atten = -10.0 * np.log10(np.maximum(
        _transmission_ratio(theta0_deg + delta, theta0_deg, cal, freq_thz), 1e-300))
    floor = float(atten[-1])
    if target_db > floor + 1e-6:
        raise ValueError(f"недостижимо на этой калибровке: максимум {floor:.2f} дБ "
                         f"(скрещенное положение, {theta0_deg + 90.0:+.3f} град)")
    atten_mono = np.maximum.accumulate(atten)      # защита от численного шума
    delta_sol = float(np.interp(target_db, atten_mono, delta))
    return {"theta_plus_deg": theta0_deg + delta_sol,
            "theta_minus_deg": theta0_deg - delta_sol,
            "delta_deg": delta_sol, "floor_db": floor}


# --- CLI --------------------------------------------------------------
def _print_forward(theta_deg: float, theta0_deg: float, cal: Calibration,
                   freq_thz: float | None) -> None:
    db_int = attenuation_db(theta_deg, theta0_deg, cal, None)
    print(f"угол WGP1 = {theta_deg:+.3f} град (delta от нуля = {theta_deg - theta0_deg:+.3f})")
    print(f"  интегральное затухание = {db_int:.2f} дБ")
    if freq_thz is not None:
        w = band_warning(freq_thz, cal)
        if w:
            print(f"  {w}")
        db_f = attenuation_db(theta_deg, theta0_deg, cal, freq_thz)
        print(f"  затухание на {freq_thz:.3f} ТГц = {db_f:.2f} дБ")


def _print_inverse(target_db: float, theta0_deg: float, cal: Calibration,
                   freq_thz: float | None) -> None:
    w = band_warning(freq_thz, cal)
    if w:
        print(w)
    sol = angle_for_db(target_db, theta0_deg, cal, freq_thz)
    metric = f"на {freq_thz:.3f} ТГц" if freq_thz is not None else "интегрально"
    print(f"желаемое затухание {target_db:.2f} дБ ({metric}), "
         f"максимум на этой калибровке {sol['floor_db']:.2f} дБ")
    print(f"  угол WGP1 = {sol['theta_plus_deg']:+.3f} град  (delta={sol['delta_deg']:+.3f})")
    print(f"  или       = {sol['theta_minus_deg']:+.3f} град  (delta={-sol['delta_deg']:+.3f})")
    print("  -- выбрать вариант ближе к текущему положению ротатора")


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--zero", type=float, required=True,
                    help="SET ZERO: показание шкалы WGP1 в СОВМЕЩЁННОМ положении "
                         "(минимум затухания), град")
    g = ap.add_mutually_exclusive_group(required=True)
    g.add_argument("--to-db", type=float, help="желаемое затухание (дБ) -> угол")
    g.add_argument("--from-angle", type=float, help="текущий угол WGP1 (град) -> затухание")
    ap.add_argument("--freq", type=float, default=None,
                    help="частота, ТГц (иначе -- интегрально по всей записанной полосе)")
    ap.add_argument("--calibration", default=None, help="путь к JSON калибровки устройства")
    args = ap.parse_args()

    cal = load_calibration(Path(args.calibration)) if args.calibration else load_calibration()
    print(f"устройство {cal.device_id}, калибровка {cal.dataset} ({cal.generated}), "
         f"P={cal.P_um:.2f} D={cal.D_um:.2f} мкм, "
         f"потери={cal.loss_db:.3f} дБ/ТГц^{cal.gamma:.2f}")
    print(f"SET ZERO: theta0 = {args.zero:+.3f} град\n")

    try:
        if args.to_db is not None:
            _print_inverse(args.to_db, args.zero, cal, args.freq)
        else:
            _print_forward(args.from_angle, args.zero, cal, args.freq)
    except ValueError as e:
        print(f"ошибка: {e}")
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
