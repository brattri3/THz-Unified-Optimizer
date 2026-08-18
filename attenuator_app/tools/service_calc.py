"""Библиотека + CLI: калькулятор аттенюатора для обслуживания в THz-TDS
спектрометре (задача C9, санкция владельца 2026-08-19, начата после C8).

Один сценарий эксплуатации, двунаправленный:
  (а) дБ -> угол:  оператор вводит желаемое затухание, получает угол(ы),
      на которые поставить WGP1;
  (б) угол -> дБ:  оператор вводит текущий угол WGP1 (по шкале ротатора),
      получает предсказанное затухание.
Метрика затухания -- интегральная (широкополосная, по умолчанию) и/или на
конкретной частоте, если она указана. Только точечное предсказание, БЕЗ
доверительного интервала (санкция владельца 2026-08-19: минимальный масштаб
задачи).

Затухание считается ПО МОЩНОСТИ (не по амплитуде поля): `attenuation_db`
возвращает `-10*log10(P/P_ref)`. `E1` ниже -- комплексная амплитуда поля
(Джонс-компонента), `|E1|^2` -- уже интенсивность/мощность; кто хочет
амплитудное отношение поля -- `power_field_ratios()` возвращает обе величины.

Два режима нормировки `mode` (владелец 2026-08-19, как в прежней работе над
моделью и в `measured_curve.blanco_angular_curve(relative=...)`):
  'relative' (по умолчанию) -- нормировка на пропускание В ТОЧКЕ theta0
      (SET ZERO): T(theta0)=1 тождественно на любой частоте (см. вывод в
      docstring `transmission_array`). Это ФИЗИКА, а не удобство -- общий
      множитель потерь и |t_perp|^4 сокращаются точно, поэтому этот режим
      устойчив к экстраполяции (см. `attenuator_app/STATE.md`, «Два ключевых
      решения»).
  'absolute' -- доля МОЩНОСТИ, падающей на аттенюатор (до WGP1), которая
      доходит до детектора. Включает СОБСТВЕННУЮ вносимую потерю пары WGP
      даже в совмещённом положении (при theta=theta0 T = |t_perp|^4-усреднённое
      < 1, а не 1) -- то, что раньше в модели называли абсолютным затуханием.

theta0 (SET ZERO) нужен в ОБОИХ режимах одинаково -- это не точка нормировки,
а координата совмещения осей WGP1/WGP2 (`d = theta_reading - theta0` в формуле
Джонса); в 'relative' она дополнительно служит точкой нормировки, в
'absolute' -- нет.

Физика -- полная Джонс-матричная модель схемы S1 (два ИДЕНТИЧНЫХ WGP) из
`attenuator_app/tools/measured_curve.py`; вывод формулы и обоснование --
`FINDINGS_measured_curve_2026-08-19.md`, п.1 (I_perp=|t_perp|^4, не |t_perp|^2
одного WGP). Параметры устройства ЗАШИТЫ в
`calibration/att_11_16_ca85_02721.json` (как получено --
`calibration/build_service_calibration.py`), включая theta0 калибровочной
серии -- он же дефолт SET ZERO ниже. Рантайм `data_pool/` не читает и зависит
только от numpy + `attenuator_app.core.blanco` (тот же принцип изоляции от
научного стека, что у клиента v0.2).

НЕ часть клиентского `attenuator_app.gui`/`cli` (v0.2/v0.3, отдельный трек) --
самостоятельный инструмент поверх модели C8.

Запуск из корня репозитория:
    .venv\\Scripts\\python.exe -m attenuator_app.tools.service_calc --to-db 12
    .venv\\Scripts\\python.exe -m attenuator_app.tools.service_calc --from-angle 35 --freq 0.8
    .venv\\Scripts\\python.exe -m attenuator_app.tools.service_calc --from-angle 35 --mode absolute
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
MODES = ("relative", "absolute")


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


def power_field_ratios(att_db: float) -> tuple[float, float]:
    """(P/P_ref, E/E_ref) по затуханию в дБ ПО МОЩНОСТИ:
    P/P_ref = 10^(-dB/10), E/E_ref = 10^(-dB/20) = sqrt(P/P_ref)."""
    return 10.0 ** (-att_db / 10.0), 10.0 ** (-att_db / 20.0)


def transmission_array(theta_deg, theta0_deg: float, cal: Calibration,
                       freq_thz: float | None = None, mode: str = "relative") -> np.ndarray:
    """Отношение МОЩНОСТЕЙ T(theta) для схемы S1 (два идентичных WGP), см.
    docstring модуля и `FINDINGS_measured_curve_2026-08-19.md` п.1 /
    `measured_curve.blanco_angular_curve` (тот же вывод, продублирован здесь в
    минимальном виде, чтобы не тянуть импорт `track_viewer` в рантайм этого
    инструмента).

        E1(theta,nu) = t_perp(nu)*cos^2(theta-theta0) + t_par(nu)*sin^2(theta-theta0)
        'relative': T = <|E1|^2> / <|t_perp|^2>            (T(theta0) == 1 тождественно)
        'absolute': T = <|t_perp|^2 * |E1|^2>               (доля входной мощности)

    `<x>` -- среднее по частоте: `freq_thz=None` -> взвешенное калибровочным
    референсным спектром `|E_ref(nu)|^2` по всей записанной полосе (теорема
    Парсеваля, FINDINGS п.4); `freq_thz=<число>` -- значение на одной частоте
    (среднее по одному элементу = само значение).

    `theta_deg` -- массив любой формы, градусы. Возвращает массив той же формы.
    """
    if mode not in MODES:
        raise ValueError(f"неизвестный режим {mode!r}, ожидается один из {MODES}")
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

    if mode == "relative":
        return wavg(np.abs(E1) ** 2) / wavg(np.abs(tp) ** 2)
    return wavg(np.abs(tp[None, :]) ** 2 * np.abs(E1) ** 2)


def attenuation_db_array(theta_deg, theta0_deg: float, cal: Calibration,
                         freq_thz: float | None = None, mode: str = "relative") -> np.ndarray:
    """Затухание в ДЕЦИБЕЛАХ ПО МОЩНОСТИ: -10*log10(T), см. `transmission_array`."""
    T = transmission_array(theta_deg, theta0_deg, cal, freq_thz, mode)
    return -10.0 * np.log10(np.maximum(T, 1e-300))


def transmission(theta_deg: float, theta0_deg: float, cal: Calibration,
                 freq_thz: float | None = None, mode: str = "relative") -> float:
    return float(transmission_array(np.array([theta_deg]), theta0_deg, cal, freq_thz, mode)[0])


def attenuation_db(theta_deg: float, theta0_deg: float, cal: Calibration,
                   freq_thz: float | None = None, mode: str = "relative") -> float:
    """Прямая задача: предсказанное затухание (дБ, по мощности) на угле `theta_deg`."""
    return float(attenuation_db_array(np.array([theta_deg]), theta0_deg, cal, freq_thz, mode)[0])


def angle_for_db(target_db: float, theta0_deg: float, cal: Calibration,
                 freq_thz: float | None = None, mode: str = "relative", n: int = 901) -> dict:
    """Обратная задача: угол(ы) WGP1 для желаемого затухания (дБ, по мощности).

    Кривая затухания(delta), delta=theta-theta0 in [0,90], монотонно растёт от
    `floor_min_db` (delta=0, точка совмещения) до `floor_max_db` (delta=90,
    скрещенное положение) -- в 'relative' `floor_min_db` == 0 тождественно, в
    'absolute' это собственная вносимая потеря пары WGP даже при совмещении.
    Симметрична по знаку delta (cos^2/sin^2 -- чётные функции), поэтому решения
    ДВА: theta0+delta и theta0-delta, физически равнозначны -- какое ближе к
    текущему положению ротатора, решает оператор (моторизации нет,
    `attenuator_app` C4_motor -- todo).
    """
    if target_db < 0:
        raise ValueError("затухание не может быть отрицательным")
    delta = np.linspace(0.0, 90.0, n)
    atten = attenuation_db_array(theta0_deg + delta, theta0_deg, cal, freq_thz, mode)
    floor_min, floor_max = float(atten[0]), float(atten[-1])
    if target_db < floor_min - 1e-6:
        raise ValueError(f"ниже минимума на этой калибровке: {floor_min:.2f} дБ "
                         f"(затухание даже в совмещённом положении WGP1||WGP2, "
                         f"режим={mode})")
    if target_db > floor_max + 1e-6:
        raise ValueError(f"недостижимо на этой калибровке: максимум {floor_max:.2f} дБ "
                         f"(скрещенное положение, {theta0_deg + 90.0:+.3f} град)")
    atten_mono = np.maximum.accumulate(atten)      # защита от численного шума
    delta_sol = float(np.interp(target_db, atten_mono, delta))
    return {"theta_plus_deg": theta0_deg + delta_sol,
            "theta_minus_deg": theta0_deg - delta_sol,
            "delta_deg": delta_sol, "floor_min_db": floor_min, "floor_max_db": floor_max}


# --- CLI --------------------------------------------------------------
_MODE_LABEL = {"relative": "относительный (к theta0, T(theta0)=1)",
              "absolute": "абсолютный (доля входной мощности)"}


def _print_forward(theta_deg: float, theta0_deg: float, cal: Calibration,
                   freq_thz: float | None, mode: str) -> None:
    db_int = attenuation_db(theta_deg, theta0_deg, cal, None, mode)
    p_r, f_r = power_field_ratios(db_int)
    print(f"угол WGP1 = {theta_deg:+.3f} град (delta от theta0 = {theta_deg - theta0_deg:+.3f})")
    print(f"  интегральное затухание = {db_int:.2f} дБ по мощности "
         f"(P/P_ref={p_r:.4g}, поле E/E_ref={f_r:.4g})")
    if freq_thz is not None:
        w = band_warning(freq_thz, cal)
        if w:
            print(f"  {w}")
        db_f = attenuation_db(theta_deg, theta0_deg, cal, freq_thz, mode)
        p_rf, f_rf = power_field_ratios(db_f)
        print(f"  затухание на {freq_thz:.3f} ТГц = {db_f:.2f} дБ по мощности "
             f"(P/P_ref={p_rf:.4g}, поле E/E_ref={f_rf:.4g})")


def _print_inverse(target_db: float, theta0_deg: float, cal: Calibration,
                   freq_thz: float | None, mode: str) -> None:
    w = band_warning(freq_thz, cal)
    if w:
        print(w)
    sol = angle_for_db(target_db, theta0_deg, cal, freq_thz, mode)
    metric = f"на {freq_thz:.3f} ТГц" if freq_thz is not None else "интегрально"
    print(f"желаемое затухание {target_db:.2f} дБ по мощности ({metric}), "
         f"диапазон на этой калибровке [{sol['floor_min_db']:.2f}, {sol['floor_max_db']:.2f}] дБ")
    print(f"  угол WGP1 = {sol['theta_plus_deg']:+.3f} град  (delta={sol['delta_deg']:+.3f})")
    print(f"  или       = {sol['theta_minus_deg']:+.3f} град  (delta={-sol['delta_deg']:+.3f})")
    print("  -- выбрать вариант ближе к текущему положению ротатора")


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--zero", type=float, default=None,
                    help="SET ZERO: показание шкалы WGP1 в СОВМЕЩЁННОМ положении "
                         "(минимум затухания), град. По умолчанию -- значение из "
                         "калибровки прибора (theta0_calibration_deg в JSON); "
                         "указать явно, если прибор переустанавливали в роторе")
    ap.add_argument("--mode", choices=MODES, default="relative",
                    help="relative (по умолчанию) -- нормировка на T(theta0)=1; "
                         "absolute -- доля входной мощности, включает собственную "
                         "потерю пары WGP даже в совмещённом положении")
    g = ap.add_mutually_exclusive_group(required=True)
    g.add_argument("--to-db", type=float, help="желаемое затухание (дБ по мощности) -> угол")
    g.add_argument("--from-angle", type=float, help="текущий угол WGP1 (град) -> затухание")
    ap.add_argument("--freq", type=float, default=None,
                    help="частота, ТГц (иначе -- интегрально по всей записанной полосе)")
    ap.add_argument("--calibration", default=None, help="путь к JSON калибровки устройства")
    args = ap.parse_args()

    cal = load_calibration(Path(args.calibration)) if args.calibration else load_calibration()
    theta0 = args.zero if args.zero is not None else cal.theta0_calibration_deg
    print(f"устройство {cal.device_id}, калибровка {cal.dataset} ({cal.generated}), "
         f"P={cal.P_um:.2f} D={cal.D_um:.2f} мкм, "
         f"потери={cal.loss_db:.3f} дБ/ТГц^{cal.gamma:.2f}")
    print(f"theta0 = {theta0:+.3f} град "
         f"({'из калибровки' if args.zero is None else 'задан вручную (--zero)'})")
    print(f"режим: {_MODE_LABEL[args.mode]}\n")

    try:
        if args.to_db is not None:
            _print_inverse(args.to_db, theta0, cal, args.freq, args.mode)
        else:
            _print_forward(args.from_angle, theta0, cal, args.freq, args.mode)
    except ValueError as e:
        print(f"ошибка: {e}")
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
