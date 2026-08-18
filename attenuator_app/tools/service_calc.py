"""Библиотека + CLI: калькулятор аттенюатора для обслуживания в THz-TDS
спектрометре (задача C9, санкция владельца 2026-08-19, начата после C8).

Один сценарий эксплуатации, двунаправленный:
  (а) дБ -> угол:  оператор вводит желаемое затухание, получает угол(ы),
      на которые поставить WGP1;
  (б) угол -> дБ:  оператор вводит текущий угол WGP1 (по шкале ротатора),
      получает предсказанное затухание.
Только точечное предсказание, БЕЗ доверительного интервала (санкция владельца
2026-08-19: минимальный масштаб задачи).

ЗНАК. Затухание -- ОТРИЦАТЕЛЬНЫЕ децибелы (владелец 2026-08-18):
`attenuation_db = 10*log10(T)`, т.е. T=1 -> 0 дБ, T=0.5 -> -3.01 дБ. Величина
положительна только в режиме усиления относительно рабочей точки (см. SET ZERO
ниже). Децибелы -- ПО МОЩНОСТИ (10*log10), не по амплитуде поля (20*log10);
`power_field_ratios()` возвращает обе величины сразу.

ДВА РАЗНЫХ НУЛЯ -- не путать (владелец 2026-08-18):
  SET OFFSET (`offset_deg`, он же theta0) -- ФИЗИЧЕСКИЙ офсет: показание шкалы
      ротатора, при котором оси WGP1 и WGP2 совмещены. Входит в саму формулу
      Джонса как `d = theta_reading - offset`. Определяется ПОДБОРОМ ПАРАМЕТРОВ
      модели и ЗАШИТ в калибровку (`theta0_calibration_deg` в JSON). Оператор
      его обычно не трогает; алгоритм автоматической калибровки офсета -- задача
      следующей версии.
  SET ZERO (`zero_deg`) -- РАБОЧАЯ ТОЧКА ОТСЧЁТА: просто точка на угловой
      кривой, относительно которой оператор считает ДОБАВОЧНОЕ затухание (или
      усиление, если двигаться к совмещённому положению). В формулу физики не
      входит. По умолчанию совпадает с SET OFFSET.

Два режима нормировки `mode`:
  'relative' (по умолчанию) -- нормировка на пропускание В ТОЧКЕ SET ZERO:
      T(zero)=1 (0 дБ) тождественно на любой частоте. При zero==offset это в
      точности прежняя нормировка «на максимум пропускания в совмещённом
      положении»: общий множитель потерь и |t_perp|^4 сокращаются ТОЧНО,
      поэтому режим устойчив к экстраполяции (см. `attenuator_app/STATE.md`,
      «Два ключевых решения»). При сдвинутом SET ZERO кривая может выходить
      выше 1 (0 дБ) -- это и есть «усиление» относительно рабочей точки.
  'absolute' -- доля мощности `P_0`, падающей на аттенюатор (до WGP1), которая
      доходит до детектора. Включает СОБСТВЕННУЮ вносимую потерю пары WGP даже
      в совмещённом положении (T(offset) < 1, порядка -0.4 дБ), поэтому 0 дБ в
      этом режиме недостижим.

Метрика `Metric` -- по какой полосе усредняется пропускание, 4 варианта:
  full         -- полная мощность: вся записанная полоса с весом |E_ref(nu)|^2
                  (теорема Парсеваля, FINDINGS п.4);
  single       -- одна частота;
  band_cw      -- полоса, заданная центром и шириной;
  band_minmax  -- полоса, заданная минимумом и максимумом.
Полосы усредняются тем же весом |E_ref(nu)|^2, но только по точкам сетки,
попавшим внутрь полосы.

Физика -- полная Джонс-матричная модель схемы S1 (два ИДЕНТИЧНЫХ WGP) из
`attenuator_app/tools/measured_curve.py`; вывод формулы и обоснование --
`FINDINGS_measured_curve_2026-08-19.md`, п.1 (I_perp=|t_perp|^4, не |t_perp|^2
одного WGP). Параметры устройства ЗАШИТЫ в
`calibration/att_11_16_ca85_02721.json` (как получено --
`calibration/build_service_calibration.py`). Рантайм `data_pool/` не читает и
зависит только от numpy + `attenuator_app.core.blanco` (тот же принцип изоляции
от научного стека, что у клиента v0.2).

НЕ часть клиентского `attenuator_app.gui`/`cli` (v0.2/v0.3, отдельный трек) --
самостоятельный инструмент поверх модели C8.

Запуск из корня репозитория:
    .venv\\Scripts\\python.exe -m attenuator_app.tools.service_calc --to-db -12
    .venv\\Scripts\\python.exe -m attenuator_app.tools.service_calc --from-angle 35 --freq 0.8
    .venv\\Scripts\\python.exe -m attenuator_app.tools.service_calc --from-angle 35 --band 0.4 1.2
    .venv\\Scripts\\python.exe -m attenuator_app.tools.service_calc --from-angle 35 --mode absolute
"""
from __future__ import annotations

import argparse
import json
import sys
from dataclasses import dataclass
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))

from attenuator_app.core.blanco import dressed_t          # noqa: E402

DEFAULT_CALIBRATION_PATH = HERE / "calibration" / "att_11_16_ca85_02721.json"

MODES = ("relative", "absolute")
MODE_LABEL = {
    "relative": "относительный (нормировка на SET ZERO, там 0 дБ)",
    "absolute": "абсолютный (доля входной мощности P_0)",
}
#: чему равна опорная мощность в знаменателе T в каждом режиме
REF_POWER_LABEL = {
    "relative": "P_zero -- мощность в точке SET ZERO",
    "absolute": "P_0 -- мощность на входе аттенюатора",
}
REF_POWER_SHORT = {"relative": "P_zero", "absolute": "P_0"}

METRIC_KINDS = ("full", "single", "band_cw", "band_minmax")


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
        #: SET OFFSET -- физический офсет совмещения WGP1/WGP2, подгоночный
        #: параметр модели, зашит; НЕ рабочая точка отсчёта (это SET ZERO)
        self.theta0_calibration_deg = float(data["theta0_calibration_deg"])
        self.at_bound = bool(data["at_bound"])
        self.freqs_ref = np.array(data["freqs_ref_thz"], dtype=float)
        self.power_ref = np.array(data["power_ref"], dtype=float)


def load_calibration(path: Path = DEFAULT_CALIBRATION_PATH) -> Calibration:
    with open(path, encoding="utf-8") as f:
        return Calibration(json.load(f))


# --- метрика: по какой полосе усредняем ------------------------------
@dataclass(frozen=True)
class Metric:
    """Спецификация полосы усреднения. `a`/`b` трактуются по `kind`:
    full -- не используются; single -- a=частота; band_cw -- a=центр, b=ширина;
    band_minmax -- a=f_min, b=f_max. Всё в ТГц."""

    kind: str = "full"
    a: float | None = None
    b: float | None = None

    def __post_init__(self):
        if self.kind not in METRIC_KINDS:
            raise ValueError(f"неизвестная метрика {self.kind!r}, ожидается одна из {METRIC_KINDS}")
        if self.kind == "single" and self.a is None:
            raise ValueError("для метрики 'одна частота' нужна частота")
        if self.kind in ("band_cw", "band_minmax") and (self.a is None or self.b is None):
            raise ValueError("для полосы нужны оба значения")
        if self.kind == "band_cw" and float(self.b) <= 0:
            raise ValueError("ширина полосы должна быть больше нуля")
        if self.kind == "band_minmax" and float(self.b) <= float(self.a):
            raise ValueError("f_max должна быть больше f_min")

    @property
    def limits(self) -> tuple[float, float] | None:
        """(lo, hi) для полосовых метрик, иначе None."""
        if self.kind == "band_cw":
            return (self.a - self.b / 2.0, self.a + self.b / 2.0)
        if self.kind == "band_minmax":
            return (float(self.a), float(self.b))
        return None

    @property
    def label(self) -> str:
        if self.kind == "full":
            return "полная мощность (вся записанная полоса)"
        if self.kind == "single":
            return f"на {self.a:.3f} ТГц"
        lo, hi = self.limits
        if self.kind == "band_cw":
            return f"полоса {self.a:.3f} +- {self.b / 2.0:.3f} ТГц ({lo:.3f}-{hi:.3f})"
        return f"полоса {lo:.3f}-{hi:.3f} ТГц"

    def resolve(self, cal: Calibration) -> tuple[np.ndarray, np.ndarray | None]:
        """(частоты, веса) для усреднения. Вес None = равномерное среднее."""
        if self.kind == "full":
            return cal.freqs_ref, cal.power_ref
        if self.kind == "single":
            return np.array([float(self.a)]), None
        lo, hi = self.limits
        m = (cal.freqs_ref >= lo) & (cal.freqs_ref <= hi)
        if not m.any():
            raise ValueError(
                f"в полосе {lo:.3f}-{hi:.3f} ТГц нет ни одной точки спектра "
                f"(сетка {cal.freqs_ref[0]:.3f}-{cal.freqs_ref[-1]:.3f} ТГц, "
                f"шаг {cal.freqs_ref[1] - cal.freqs_ref[0]:.4f} ТГц)")
        return cal.freqs_ref[m], cal.power_ref[m]

    def warning(self, cal: Calibration) -> str | None:
        """Предупреждение об экстраполяции за откалиброванную полосу."""
        lo_c, hi_c = cal.band_thz
        if self.kind == "single":
            if lo_c <= self.a <= hi_c:
                return None
            return (f"[!] {self.a:.3f} ТГц вне откалиброванной полосы "
                    f"{lo_c:.2f}-{hi_c:.2f} ТГц -- значение экстраполировано, не измерено")
        lim = self.limits
        if lim is None:
            return None
        lo, hi = lim
        if lo >= lo_c and hi <= hi_c:
            return None
        return (f"[!] полоса {lo:.3f}-{hi:.3f} ТГц частично вне откалиброванной "
                f"{lo_c:.2f}-{hi_c:.2f} ТГц -- край экстраполирован, не измерен")


FULL = Metric("full")


def band_warning(metric: Metric, cal: Calibration) -> str | None:
    """Обёртка для совместимости вызова из GUI."""
    return metric.warning(cal)


def power_field_ratios(att_db: float) -> tuple[float, float]:
    """(P/P_ref, E/E_ref) по затуханию в дБ ПО МОЩНОСТИ (отрицательному):
    P/P_ref = 10^(dB/10), E/E_ref = 10^(dB/20) = sqrt(P/P_ref)."""
    return 10.0 ** (att_db / 10.0), 10.0 ** (att_db / 20.0)


# --- физика ------------------------------------------------------------
def transmission_array(theta_deg, offset_deg: float, cal: Calibration,
                       metric: Metric = FULL, mode: str = "relative",
                       zero_deg: float | None = None) -> np.ndarray:
    """Отношение МОЩНОСТЕЙ T(theta) для схемы S1 (два идентичных WGP), см.
    docstring модуля и `FINDINGS_measured_curve_2026-08-19.md` п.1 /
    `measured_curve.blanco_angular_curve` (тот же вывод, продублирован здесь в
    минимальном виде, чтобы не тянуть импорт `track_viewer` в рантайм).

        E1(theta,nu) = t_perp(nu)*cos^2(d) + t_par(nu)*sin^2(d),  d = theta - offset
        'relative': T = <|E1(theta)|^2> / <|E1(zero)|^2>     (T(zero) == 1)
        'absolute': T = <|t_perp|^2 * |E1(theta)|^2>          (доля P_0)

    `<x>` -- среднее по частоте с весом |E_ref(nu)|^2 по точкам, отобранным
    метрикой (для одной частоты -- само значение).

    `offset_deg` -- SET OFFSET (физика), `zero_deg` -- SET ZERO (рабочая точка
    отсчёта, по умолчанию = offset). `theta_deg` -- массив любой формы, градусы;
    возвращает массив той же формы.
    """
    if mode not in MODES:
        raise ValueError(f"неизвестный режим {mode!r}, ожидается один из {MODES}")
    if zero_deg is None:
        zero_deg = offset_deg
    freqs, weight = metric.resolve(cal)
    tp, ta, _ = dressed_t(freqs, cal.P_um, cal.D_um, loss_factor=cal.loss_db, gamma=cal.gamma)

    def field(th) -> np.ndarray:
        d = np.deg2rad(np.asarray(th, dtype=float) - offset_deg)
        c2, s2 = np.cos(d) ** 2, np.sin(d) ** 2
        return tp[None, :] * c2[..., None] + ta[None, :] * s2[..., None]

    def wavg(x):
        return np.mean(x, axis=-1) if weight is None else np.average(x, axis=-1, weights=weight)

    e1 = field(theta_deg)
    if mode == "relative":
        p_zero = float(wavg(np.abs(field(np.array([zero_deg]))) ** 2)[0])
        return wavg(np.abs(e1) ** 2) / p_zero
    return wavg(np.abs(tp[None, :]) ** 2 * np.abs(e1) ** 2)


def attenuation_db_array(theta_deg, offset_deg: float, cal: Calibration,
                         metric: Metric = FULL, mode: str = "relative",
                         zero_deg: float | None = None) -> np.ndarray:
    """Затухание в ДЕЦИБЕЛАХ ПО МОЩНОСТИ, ОТРИЦАТЕЛЬНЫХ: 10*log10(T)."""
    T = transmission_array(theta_deg, offset_deg, cal, metric, mode, zero_deg)
    return 10.0 * np.log10(np.maximum(T, 1e-300))


def transmission(theta_deg: float, offset_deg: float, cal: Calibration,
                 metric: Metric = FULL, mode: str = "relative",
                 zero_deg: float | None = None) -> float:
    return float(transmission_array(np.array([theta_deg]), offset_deg, cal,
                                    metric, mode, zero_deg)[0])


def attenuation_db(theta_deg: float, offset_deg: float, cal: Calibration,
                   metric: Metric = FULL, mode: str = "relative",
                   zero_deg: float | None = None) -> float:
    """Прямая задача: предсказанное затухание (отрицательные дБ по мощности)."""
    return float(attenuation_db_array(np.array([theta_deg]), offset_deg, cal,
                                      metric, mode, zero_deg)[0])


def angle_for_db(target_db: float, offset_deg: float, cal: Calibration,
                 metric: Metric = FULL, mode: str = "relative",
                 zero_deg: float | None = None, n: int = 901) -> dict:
    """Обратная задача: угол(ы) WGP1 для желаемого затухания (отрицательные дБ).

    Кривая затухания(delta), delta=theta-offset in [0,90], монотонно УБЫВАЕТ от
    `db_max` (delta=0, совмещённое положение) до `db_min` (delta=90,
    скрещенное). Симметрична по знаку delta (cos^2/sin^2 -- чётные), поэтому
    решения ДВА: offset+delta и offset-delta, физически равнозначны -- какое
    ближе к текущему положению ротатора, решает оператор (моторизации нет,
    `attenuator_app` C4_motor -- todo).
    """
    delta = np.linspace(0.0, 90.0, n)
    atten = attenuation_db_array(offset_deg + delta, offset_deg, cal, metric, mode, zero_deg)
    db_max, db_min = float(atten[0]), float(atten[-1])
    if target_db > db_max + 1e-6:
        hint = ""
        if target_db > 0 and db_min - 1e-6 <= -target_db <= db_max + 1e-6:
            hint = f"; затухание задаётся ОТРИЦАТЕЛЬНЫМ числом -- возможно, нужно {-target_db:.2f} дБ"
        raise ValueError(f"выше максимума на этой калибровке: {db_max:.2f} дБ "
                         f"(совмещённое положение WGP1||WGP2, режим={mode}){hint}")
    if target_db < db_min - 1e-6:
        raise ValueError(f"недостижимо на этой калибровке: минимум {db_min:.2f} дБ "
                         f"(скрещенное положение, {offset_deg + 90.0:+.3f} град)")
    atten_mono = np.minimum.accumulate(atten)      # защита от численного шума
    delta_sol = float(np.interp(target_db, atten_mono[::-1], delta[::-1]))
    return {"theta_plus_deg": offset_deg + delta_sol,
            "theta_minus_deg": offset_deg - delta_sol,
            "delta_deg": delta_sol, "db_max": db_max, "db_min": db_min}


# --- CLI --------------------------------------------------------------
def _print_line(prefix: str, db: float, mode: str) -> None:
    p_r, f_r = power_field_ratios(db)
    print(f"{prefix} = {db:+.2f} дБ по мощности  "
          f"(T = P/{REF_POWER_SHORT[mode]} = {p_r * 100:.3f} %, "
          f"поле E/E_ref = {f_r:.4g})")


def _print_forward(theta_deg: float, offset_deg: float, zero_deg: float,
                   cal: Calibration, metric: Metric, mode: str) -> None:
    w = metric.warning(cal)
    if w:
        print(w)
    db = attenuation_db(theta_deg, offset_deg, cal, metric, mode, zero_deg)
    print(f"угол WGP1 = {theta_deg:+.3f} град "
          f"(delta от SET OFFSET = {theta_deg - offset_deg:+.3f}, "
          f"от SET ZERO = {theta_deg - zero_deg:+.3f})")
    _print_line(f"  затухание, {metric.label}", db, mode)
    db_zero = attenuation_db(zero_deg, offset_deg, cal, metric, mode, zero_deg)
    d_zero = db - db_zero
    kind = "усиление" if d_zero > 0 else "затухание"
    print(f"  относительно SET ZERO ({zero_deg:+.3f} град): {d_zero:+.2f} дБ ({kind})")


def _print_inverse(target_db: float, offset_deg: float, zero_deg: float,
                   cal: Calibration, metric: Metric, mode: str) -> None:
    w = metric.warning(cal)
    if w:
        print(w)
    sol = angle_for_db(target_db, offset_deg, cal, metric, mode, zero_deg)
    print(f"желаемое затухание {target_db:+.2f} дБ по мощности ({metric.label}), "
          f"диапазон на этой калибровке [{sol['db_min']:.2f}, {sol['db_max']:.2f}] дБ")
    print(f"  угол WGP1 = {sol['theta_plus_deg']:+.3f} град  (delta={sol['delta_deg']:+.3f})")
    print(f"  или       = {sol['theta_minus_deg']:+.3f} град  (delta={-sol['delta_deg']:+.3f})")
    print("  -- выбрать вариант ближе к текущему положению ротатора")


def metric_from_args(args) -> Metric:
    given = [args.freq is not None, args.band is not None,
             args.band_center is not None or args.band_width is not None]
    if sum(given) > 1:
        raise ValueError("--freq, --band и --band-center/--band-width взаимоисключающи")
    if args.freq is not None:
        return Metric("single", args.freq)
    if args.band is not None:
        return Metric("band_minmax", args.band[0], args.band[1])
    if args.band_center is not None or args.band_width is not None:
        if args.band_center is None or args.band_width is None:
            raise ValueError("--band-center и --band-width задаются вместе")
        return Metric("band_cw", args.band_center, args.band_width)
    return FULL


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--offset", type=float, default=None,
                    help="SET OFFSET (theta0): показание шкалы WGP1 в СОВМЕЩЁННОМ "
                         "положении, град. Физический параметр модели, по умолчанию "
                         "берётся из зашитой калибровки прибора")
    ap.add_argument("--zero", type=float, default=None,
                    help="SET ZERO: рабочая точка отсчёта, град. Относительно неё "
                         "считается добавочное затухание/усиление; в режиме relative "
                         "она же точка нормировки. По умолчанию = SET OFFSET")
    ap.add_argument("--mode", choices=MODES, default="relative",
                    help="relative (по умолчанию) -- нормировка на SET ZERO; "
                         "absolute -- доля входной мощности P_0, включает собственную "
                         "потерю пары WGP даже в совмещённом положении")
    g = ap.add_mutually_exclusive_group(required=True)
    g.add_argument("--to-db", type=float,
                   help="желаемое затухание (ОТРИЦАТЕЛЬНЫЕ дБ по мощности) -> угол")
    g.add_argument("--from-angle", type=float, help="текущий угол WGP1 (град) -> затухание")
    ap.add_argument("--freq", type=float, default=None, help="метрика: одна частота, ТГц")
    ap.add_argument("--band", type=float, nargs=2, metavar=("FMIN", "FMAX"),
                    help="метрика: полоса по минимуму и максимуму, ТГц")
    ap.add_argument("--band-center", type=float, default=None,
                    help="метрика: центр полосы, ТГц (вместе с --band-width)")
    ap.add_argument("--band-width", type=float, default=None,
                    help="метрика: ширина полосы, ТГц (вместе с --band-center)")
    ap.add_argument("--calibration", default=None, help="путь к JSON калибровки устройства")
    args = ap.parse_args()

    cal = load_calibration(Path(args.calibration)) if args.calibration else load_calibration()
    offset = args.offset if args.offset is not None else cal.theta0_calibration_deg
    zero = args.zero if args.zero is not None else offset

    try:
        metric = metric_from_args(args)
    except ValueError as e:
        print(f"ошибка: {e}")
        return 1

    print(f"устройство {cal.device_id}, калибровка {cal.dataset} ({cal.generated}), "
          f"P={cal.P_um:.2f} D={cal.D_um:.2f} мкм, "
          f"потери={cal.loss_db:.3f} дБ/ТГц^{cal.gamma:.2f}")
    print(f"SET OFFSET = {offset:+.3f} град "
          f"({'из калибровки' if args.offset is None else 'задан вручную'})")
    print(f"SET ZERO   = {zero:+.3f} град "
          f"({'= SET OFFSET' if args.zero is None else 'задан вручную'})")
    print(f"режим: {MODE_LABEL[args.mode]}; метрика: {metric.label}\n")

    try:
        if args.to_db is not None:
            _print_inverse(args.to_db, offset, zero, cal, metric, args.mode)
        else:
            _print_forward(args.from_angle, offset, zero, cal, metric, args.mode)
    except ValueError as e:
        print(f"ошибка: {e}")
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
