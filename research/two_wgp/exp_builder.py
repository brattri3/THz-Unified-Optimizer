"""
Приватный билдер экспериментальной сетки Сессии A (зона `research/two_wgp/`).

ЗАЧЕМ. `fit_lib.build_experiment` жёстко применяет `config.ANGLES_LIMIT_2D=(0,90)`, из-за
чего все фиты Фаз 1-3 и A2 шли на 10 углах из 19: отрицательная полуось 356att
(-90..-10) не использовалась НИ РАЗУ. Задача A6 (тест симметрии U(theta) vs U(-theta))
требует полного диапазона. Продакшн-конфиг read-only (гардрейл RESEARCH_PLAN §6),
`fit_lib.py` -- зона Сессии B, поэтому билдер живёт здесь.

ОТЛИЧИЕ ОТ `fit_lib.build_experiment` РОВНО ОДНО (плюс один опциональный флаг):
  * `angles_limit` -- параметр, а не константа конфига (`"config"` = взять из конфига,
    `None` = все углы, кортеж = свой диапазон);
  * `mask_angles` -- по каким углам считать водяную маску (см. ВАЖНО ниже).
Всё остальное -- те же продакшн-функции (`get_transmission_spectra`,
`find_auto_water_mask`) и те же границы `config.F_MIN` / `config.F_MAX_2D`, в том же
порядке операций, чтобы числа оставались сравнимыми с Фазами 1-3.

ВАЖНО (находка A6, передана B в HANDOFFS). Водяная маска считается по среднему ПО УГЛАМ
(`fit_lib.py:81`), поэтому она ЗАВИСИТ ОТ НАБОРА УГЛОВ: маска, посчитанная по 19 углам,
не равна маске по 10. Значит «тот же датасет на полном диапазоне» молча меняет ещё и
набор частот. Чтобы отделить эффект углов от эффекта маски, `mask_angles` позволяет
считать маску по фиксированному подмножеству:
    mask_angles=None      -- как в fit_lib: по выбранным углам (маска «плывёт»);
    mask_angles=(0., 90.) -- маска ровно как в Фазах 1-3, углы полные (нужный режим A6).

ГАРАНТИЯ. `python research/two_wgp/exp_builder.py --verify` проверяет БИТ-В-БИТ
совпадение с `fit_lib.build_experiment` при `angles_limit="config"`, `mask_angles=None`
на всех доступных датасетах (допуск РОВНО 0). Расхождение -> ненулевой код возврата.

Сессия A. Ответ на хэндофф B «под A6 не форкать билдер»: предложение параметризовать
общий `build_experiment(..., angles_limit=None)` принято; когда B его сделает, этот
модуль становится тонкой обёрткой или удаляется.
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[2]  # THz-Unified-Optimizer
sys.path.insert(0, str(REPO))
sys.path.insert(0, str(REPO / "research" / "experiments"))

from unified_optimizer import config                                   # noqa: E402
from unified_optimizer.data_manager import DataManager                 # noqa: E402
from unified_optimizer.optimizer_2d import get_transmission_spectra    # noqa: E402
from unified_optimizer.utils import find_auto_water_mask               # noqa: E402

_CONFIG = "config"


def select_angles(data_dict, angles_limit=_CONFIG):
    """Отсортированный список углов с применённым фильтром.

    angles_limit: "config" -> config.ANGLES_LIMIT_2D (поведение fit_lib);
                  None     -> все углы;
                  (lo, hi) -> свой диапазон (включительно).
    Углы приходят уже нормализованными DataManager'ом в (-180, 180]:
    файл `..._270deg_...` -> -90, `..._350deg_...` -> -10.
    """
    angles = sorted(data_dict.keys())
    if angles_limit == _CONFIG:
        angles_limit = config.ANGLES_LIMIT_2D
    if angles_limit is not None:
        lo, hi = angles_limit
        angles = [a for a in angles if lo <= a <= hi]
    return angles


def build_experiment(data_dict, angles_limit=_CONFIG, mask_angles=None):
    """2D-сетка эксперимента + водяная маска. См. отличия от fit_lib в докстринге модуля.

    Возвращает (angles_val, analysis_freqs, exp, valid, meta).
    Первые четыре элемента совпадают с `fit_lib.build_experiment` при значениях
    по умолчанию (проверяется `--verify`).
    """
    angles = select_angles(data_dict, angles_limit)
    if not angles:
        raise ValueError("после фильтра angles_limit не осталось ни одного угла")

    freqs_common, individual = None, {}
    for a in angles:
        t_s, E_s, t_b, E_b = data_dict[a]
        f, s_s, s_b, trans = get_transmission_spectra(t_s, E_s, t_b, E_b)
        individual[a] = trans
        if freqs_common is None:
            freqs_common = f

    f_max = getattr(config, "F_MAX_2D", config.F_MAX)
    idx = np.where((freqs_common >= config.F_MIN) & (freqs_common <= f_max))[0]
    analysis_freqs = freqs_common[idx]
    angles_val = np.array(angles)
    exp = np.zeros((len(angles_val), len(analysis_freqs)), dtype=complex)
    for i, a in enumerate(angles):
        exp[i, :] = individual[a][idx]

    # --- водяная маска ---
    if mask_angles is None:
        mask_rows = np.arange(len(angles_val))          # поведение fit_lib
    else:
        lo, hi = mask_angles
        mask_rows = np.array([i for i, a in enumerate(angles) if lo <= a <= hi])
        if mask_rows.size == 0:
            raise ValueError("mask_angles не покрывает ни одного из выбранных углов")
    Ymean_db = np.mean(20 * np.log10(np.maximum(np.abs(exp[mask_rows, :]), 1e-12)), axis=0)
    water_mask, intervals = find_auto_water_mask(analysis_freqs, Ymean_db)

    valid = np.ones_like(exp, dtype=bool)
    for i in range(len(angles_val)):
        valid[i, ~water_mask] = False

    meta = {
        "n_angles": int(len(angles_val)),
        "angles": [float(a) for a in angles_val],
        "n_freqs": int(len(analysis_freqs)),
        "f_min_thz": float(analysis_freqs[0]),
        "f_max_thz": float(analysis_freqs[-1]),
        "n_valid_freqs": int(np.sum(water_mask)),
        "n_valid_points": int(np.sum(valid)),
        "water_mask": water_mask,
        "masked_intervals": [(float(a), float(b)) for a, b in intervals],
        "mask_from_angles": [float(angles_val[i]) for i in mask_rows],
        "angles_limit": (None if angles_limit is None else
                         (list(config.ANGLES_LIMIT_2D) if angles_limit == _CONFIG
                          else list(angles_limit))),
    }
    return angles_val, analysis_freqs, exp, valid, meta


def load(dataset, angles_limit=_CONFIG, mask_angles=None, dm=None):
    """Удобная обёртка: датасет по имени -> build_experiment."""
    dm = dm or DataManager(config.DATA_DIR)
    data = dm.get_data_for_dataset(dataset)
    if not data:
        raise ValueError(f"датасет {dataset!r} не найден в {config.DATA_DIR}")
    return build_experiment(data, angles_limit=angles_limit, mask_angles=mask_angles)


def verify(datasets=None, verbose=True):
    """БИТ-В-БИТ сверка с fit_lib.build_experiment на дефолтных настройках. True = ок."""
    import fit_lib  # зона B, только чтение

    dm = DataManager(config.DATA_DIR)
    datasets = datasets or dm.get_datasets()
    ok = True
    for ds in datasets:
        data = dm.get_data_for_dataset(ds)
        if not data:
            continue
        try:
            ref = fit_lib.build_experiment(data)
        except Exception as exc:                     # датасет с одним углом и т.п.
            if verbose:
                print(f"  {ds:<18} SKIP ({type(exc).__name__}: {exc})")
            continue
        mine = build_experiment(data)[:4]
        same = all(np.array_equal(a, b) for a, b in zip(ref, mine))
        ok &= same
        if verbose:
            n_ang, n_f = ref[2].shape
            print(f"  {ds:<18} {'OK  ' if same else 'FAIL'} "
                  f"angles={n_ang:>3} freqs={n_f:>4} valid={int(np.sum(ref[3])):>5}")
    return ok


if __name__ == "__main__":
    import argparse

    ap = argparse.ArgumentParser(description=__doc__.split("\n")[1])
    ap.add_argument("--verify", action="store_true",
                    help="бит-в-бит сверка с fit_lib.build_experiment (допуск 0)")
    ap.add_argument("--coverage", action="store_true",
                    help="угловое покрытие датасетов: сколько углов теряет ANGLES_LIMIT_2D")
    args = ap.parse_args()

    if args.coverage:
        dm = DataManager(config.DATA_DIR)
        print(f"ANGLES_LIMIT_2D = {config.ANGLES_LIMIT_2D}")
        print(f"{'dataset':<18} {'всего':>6} {'в фите':>7} {'потеряно':>9}  углы (норм., град)")
        for ds in dm.get_datasets():
            data = dm.get_data_for_dataset(ds)
            full = select_angles(data, None)
            lim = select_angles(data, _CONFIG)
            print(f"{ds:<18} {len(full):>6} {len(lim):>7} {len(full)-len(lim):>9}  "
                  f"{[int(a) for a in full]}")

    if args.verify or not args.coverage:
        print("Бит-в-бит сверка exp_builder.build_experiment vs fit_lib.build_experiment:")
        good = verify()
        print("ИТОГ:", "совпадение бит-в-бит (допуск 0)" if good else "РАСХОЖДЕНИЕ!")
        sys.exit(0 if good else 1)
