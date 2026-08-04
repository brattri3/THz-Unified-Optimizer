"""A7 — ввод датасета `purewave` в оборот: качество данных, сшивка двух дней, фит.

ЗАЧЕМ. `purewave` числился BLOCKED (A3, 2026-08-04): в `data_pool/purewave/` лежали одни
заглушки нулевого размера. Владелец домерил его за два дня (28.07 и 03.08.2026), имена
приведены к конвенции `normalize_purewave_data.py`. Датасет ценен тем, что это ТРЕТИЙ
образец категории «одиночный WGP между плёнками TYDEX» и единственный с угловым покрытием
-95..110 град (24 угла) — шире, чем у specac (13) и att-11-16 (19).

ОСОБОЕ ОБСТОЯТЕЛЬСТВО (метаданные владельца 2026-08-04). Этот поляризатор ВОССТАНОВЛЕН:
центральная зона потеряла натяжение (отслоение клея от оправы), провода провисли; владелец
переклеил и деформацией оправы подтянул центральные провода. Периодичность в центре
восстановлена НЕ полностью, по краям остались заводские зоны с хорошей периодичностью.
⇒ `purewave` — натурный тест на связь «регулярность решётки ↔ загадка D_eff»: если
D_eff/D_phys у него выпадает из коридора 0.37..0.59, найденного A3 на шести образцах,
это первое прямое указание, что нерегулярность структуры вносит вклад (ср. HN4, отвергнутую
на СИНТЕТИЧЕСКОЙ нерегулярности — здесь нерегулярность настоящая и известна).

ЧТО ДЕЛАЕТ СКРИПТ
-----------------
1. КАЧЕСТВО (`--quality`). Дрейф фона по 10 трассам bg (пиковое поле, положение пика,
   спектральная амплитуда) — внутри дня и МЕЖДУ днями; динамический диапазон по углам;
   положение реального минимума пропускания.
2. СШИВКА ДНЕЙ (входит в `--quality`). Пропускание π-периодично по углу, поэтому
   theta = -95 град (день 1) и theta = +85 град (день 2) — ФИЗИЧЕСКИ ОДНА И ТА ЖЕ
   ориентация. Угол +85 не измерялся, но измерены 84 и 86 ⇒ сравнение с интерполяцией.
   То же для -85 ↔ +95 (интерполяция 94/96). Это ПРЯМАЯ мера межднев­ного дрейфа, не
   требующая никакой модели.
3. ГАРМОНИКИ A0 (`--harmonics`). Модельно-независимая утечка eta и офсет theta0 —
   на полном диапазоне и на продакшн-обрезке (0,90), чтобы число было сравнимо с A3.
4. ФИТ M5 (`--fit`). Мультистарт по сетке через `a3_eps_stability.analyse` — те же 21
   старт, тот же тест бассейнов; ключевой выход — D_eff/D_phys и eta0.

Гардрейлы: `data_pool/` и `unified_optimizer/` — чтение; ядро (`fit_lib`, `model_core`,
`t6_hn*`) — чтение; артефакты — в `research/results/two_wgp/`.

Запуск ИЗ КОРНЯ репозитория:
    .venv/Scripts/python.exe research/two_wgp/a7_purewave.py --quality --harmonics
    .venv/Scripts/python.exe research/two_wgp/a7_purewave.py            # всё, включая фит
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(REPO / "research" / "experiments"))
sys.path.insert(0, str(REPO))

import exp_builder                                                    # noqa: E402
import malus_linearization as ml                                      # noqa: E402
import fit_lib                                                        # noqa: E402  зона B, ЧТЕНИЕ
from unified_optimizer import config                                  # noqa: E402
from unified_optimizer.data_manager import DataManager                # noqa: E402
from unified_optimizer.optimizer_2d import get_transmission_spectra   # noqa: E402

DATASET = "purewave"
OUT = REPO / "research" / "results" / "two_wgp"
OUT.mkdir(parents=True, exist_ok=True)

# Разбиение по дням — из `normalize_purewave_data.MAPPING` (018_META: граница дней).
DAY1_ANGLES = [-95, -85, -60, -40, -20, 0, 20, 40, 55, 65, 72, 80]
DAY2_ANGLES = [82, 84, 86, 88, 90, 92, 94, 96, 98, 100, 105, 110]
# Пары «физически один и тот же угол» через pi-периодичность: (день1, день2-эквивалент)
PI_PAIRS = [(-95.0, 85.0), (-85.0, 95.0)]
REP2_ANGLES = [94, 96, 98, 100, 105, 110]


# ------------------------------------------------------------------ утилиты
def load_raw():
    """Сырые пары (t,E) по (угол, rep) — минуя усреднение реплик в DataManager."""
    dm = DataManager(config.DATA_DIR)
    return {k[1:]: v for k, v in dm.raw_data.items() if k[0] == DATASET}


def trans_by_angle(rep_filter=None):
    """{угол: (freqs, комплексное пропускание)} с выбором реплики."""
    raw = load_raw()
    out, freqs = {}, None
    for (angle, rep), d in sorted(raw.items()):
        if rep_filter is not None and rep != rep_filter:
            continue
        if "sig" not in d or "bg" not in d:
            continue
        (t_s, E_s), (t_b, E_b) = d["sig"], d["bg"]
        f, _s, _b, tr = get_transmission_spectra(t_s, E_s, t_b, E_b)
        freqs = f if freqs is None else freqs
        out[(angle, rep)] = tr
    return freqs, out


def band(freqs):
    """Рабочая полоса анализа (та же, что у билдера): F_MIN..F_MAX_2D."""
    fmax = getattr(config, "F_MAX_2D", config.F_MAX)
    return np.where((freqs >= config.F_MIN) & (freqs <= fmax))[0]


def interp_complex(angles, values, target):
    """Линейная интерполяция комплексного спектра по углу (по Re и Im отдельно)."""
    a = np.asarray(angles, float)
    order = np.argsort(a)
    a = a[order]
    V = np.asarray(values)[order]
    re = np.array([np.interp(target, a, V[:, j].real) for j in range(V.shape[1])])
    im = np.array([np.interp(target, a, V[:, j].imag) for j in range(V.shape[1])])
    return re + 1j * im


# --------------------------------------------------------------- 1. качество
def quality(verbose=True):
    raw = load_raw()
    # --- дрейф фона -------------------------------------------------------
    # У каждого угла свой файл bg, но физически их всего 10 (копии одного фона).
    # Группируем по хешу содержимого — это восстанавливает исходные 10 трасс.
    seen, bgs = {}, []
    for (angle, rep), d in sorted(raw.items()):
        if "bg" not in d:
            continue
        t_b, E_b = d["bg"]
        key = (float(E_b[0]), float(E_b[-1]), float(np.sum(E_b)))
        if key in seen:
            continue
        seen[key] = True
        day = 1 if angle in DAY1_ANGLES else 2
        i_pk = int(np.argmax(np.abs(E_b)))
        bgs.append({"day": day, "first_angle": float(angle),
                    "peak_abs": float(np.abs(E_b[i_pk])),
                    "t_peak_ps": float(t_b[i_pk]),
                    "rms": float(np.sqrt(np.mean(E_b ** 2)))})
    d1 = [b for b in bgs if b["day"] == 1]
    d2 = [b for b in bgs if b["day"] == 2]

    def stat(rows, key):
        v = np.array([r[key] for r in rows], float)
        return {"n": len(v), "mean": float(v.mean()), "std": float(v.std(ddof=1)) if len(v) > 1 else 0.0,
                "min": float(v.min()), "max": float(v.max()),
                "spread_pct": float(100 * (v.max() - v.min()) / abs(v.mean())) if v.mean() else None}

    bg_drift = {
        "n_unique_backgrounds": len(bgs),
        "day1": {"peak_abs": stat(d1, "peak_abs"), "t_peak_ps": stat(d1, "t_peak_ps")},
        "day2": {"peak_abs": stat(d2, "peak_abs"), "t_peak_ps": stat(d2, "t_peak_ps")},
        "day2_over_day1_peak": (float(np.mean([b["peak_abs"] for b in d2])
                                      / np.mean([b["peak_abs"] for b in d1]))
                                if d1 and d2 else None),
        "day2_minus_day1_tpeak_ps": (float(np.mean([b["t_peak_ps"] for b in d2])
                                           - np.mean([b["t_peak_ps"] for b in d1]))
                                     if d1 and d2 else None),
    }

    # --- динамический диапазон и положение тени --------------------------
    freqs, tr = trans_by_angle(rep_filter=1)
    idx = band(freqs)
    per_angle = []
    for (angle, rep), t in sorted(tr.items()):
        u = np.abs(t[idx]) ** 2
        per_angle.append({"angle": float(angle),
                          "U_parseval": float(u.sum()),
                          "amp_median": float(np.median(np.abs(t[idx])))})
    U = np.array([p["U_parseval"] for p in per_angle])
    A = np.array([p["angle"] for p in per_angle])
    i_max, i_min = int(np.argmax(U)), int(np.argmin(U))
    # параболическая интерполяция вершины по трём точкам вокруг минимума
    theta_min_parab = None
    if 0 < i_min < len(U) - 1:
        x = A[i_min - 1:i_min + 2]
        y = np.log(U[i_min - 1:i_min + 2])
        c = np.polyfit(x, y, 2)
        if c[0] > 0:
            theta_min_parab = float(-c[1] / (2 * c[0]))
    dyn = {
        "U_max": float(U[i_max]), "angle_U_max": float(A[i_max]),
        "U_min": float(U[i_min]), "angle_U_min": float(A[i_min]),
        "extinction_ratio": float(U[i_max] / U[i_min]),
        "extinction_db": float(10 * np.log10(U[i_max] / U[i_min])),
        "angle_U_min_parabolic": theta_min_parab,
        "per_angle": per_angle,
    }

    # --- сшивка дней через pi-периодичность -------------------------------
    day2_pairs = [(a, r) for (a, r) in tr if a in DAY2_ANGLES]
    ang2 = [a for a, _ in day2_pairs]
    V2 = np.array([tr[(a, r)] for a, r in day2_pairs])
    stitch = []
    for a1, a2_equiv in PI_PAIRS:
        key1 = (a1, 1)
        if key1 not in tr:
            continue
        if not (min(ang2) <= a2_equiv <= max(ang2)):
            continue
        v1 = tr[key1][idx]
        v2 = interp_complex(ang2, V2[:, idx], a2_equiv)
        u1, u2 = np.abs(v1) ** 2, np.abs(v2) ** 2
        stitch.append({
            "angle_day1": a1, "angle_day2_equiv": a2_equiv,
            "note": "pi-периодичность: одна и та же ориентация проводов",
            "U_day1": float(u1.sum()), "U_day2": float(u2.sum()),
            "ratio_U": float(u1.sum() / u2.sum()),
            "ratio_db": float(10 * np.log10(u1.sum() / u2.sum())),
            "amp_ratio_median": float(np.median(np.abs(v1) / np.abs(v2))),
            "phase_diff_median_deg": float(np.degrees(np.median(
                np.angle(v1 / v2)))),
        })

    # --- воспроизводимость реплик (обратный ход того же дня) --------------
    _f, tr2 = trans_by_angle(rep_filter=2)
    reps = []
    for a in REP2_ANGLES:
        k1, k2 = (float(a), 1), (float(a), 2)
        if k1 not in tr or k2 not in tr2:
            continue
        v1, v2 = tr[k1][idx], tr2[k2][idx]
        u1, u2 = float((np.abs(v1) ** 2).sum()), float((np.abs(v2) ** 2).sum())
        reps.append({"angle": float(a), "U_rep1": u1, "U_rep2": u2,
                     "ratio_db": float(10 * np.log10(u1 / u2)),
                     "amp_ratio_median": float(np.median(np.abs(v1) / np.abs(v2)))})
    rep_spread_db = (float(np.std([r["ratio_db"] for r in reps], ddof=1))
                     if len(reps) > 1 else None)

    out = {"background_drift": bg_drift, "dynamic_range": dyn,
           "day_stitching": stitch,
           "replicates": {"per_angle": reps, "std_ratio_db": rep_spread_db}}

    if verbose:
        print("\n--- КАЧЕСТВО ДАННЫХ purewave ---")
        print(f"  уникальных фонов: {bg_drift['n_unique_backgrounds']} "
              f"(день1 {len(d1)}, день2 {len(d2)})")
        for tag, dd in (("день1", bg_drift["day1"]), ("день2", bg_drift["day2"])):
            print(f"    {tag}: пик |E| = {dd['peak_abs']['mean']:.4f} "
                  f"(разброс {dd['peak_abs']['spread_pct']:.1f}%), "
                  f"t_пика = {dd['t_peak_ps']['mean']:.3f} пс "
                  f"(разброс {dd['t_peak_ps']['max'] - dd['t_peak_ps']['min']:.3f} пс)")
        print(f"    день2/день1 по амплитуде фона: {bg_drift['day2_over_day1_peak']:.4f} "
              f"({20 * np.log10(bg_drift['day2_over_day1_peak']):+.2f} дБ), "
              f"сдвиг пика {bg_drift['day2_minus_day1_tpeak_ps']:+.3f} пс")
        print(f"  динамический диапазон: max при {dyn['angle_U_max']:+.0f} град, "
              f"min при {dyn['angle_U_min']:+.0f} град "
              f"(парабола: {dyn['angle_U_min_parabolic']:+.2f} град), "
              f"экстинкция {dyn['extinction_db']:.1f} дБ")
        print("  СШИВКА ДНЕЙ (pi-периодичность, модельно-независимо):")
        for s in stitch:
            print(f"    U({s['angle_day1']:+.0f}) / U({s['angle_day2_equiv']:+.0f}) = "
                  f"{s['ratio_U']:.4f}  ({s['ratio_db']:+.2f} дБ), "
                  f"фаза {s['phase_diff_median_deg']:+.1f} град")
        if reps:
            print(f"  реплики rep1/rep2 ({len(reps)} углов): "
                  f"разброс отношения {rep_spread_db:.3f} дБ, "
                  f"диапазон {min(r['ratio_db'] for r in reps):+.2f}.."
                  f"{max(r['ratio_db'] for r in reps):+.2f} дБ")
    return out


# Наборы углов для раздельного анализа. Границы подобраны так, чтобы РОВНО совпасть
# с днями: день 1 = -95..80, день 2 = 82..110 (см. MAPPING). Продакшн-диапазон (0,90)
# СМЕШИВАЕТ дни (0..80 из дня 1, 82..90 из дня 2) — это отдельно отмечается в выводе.
RANGES = (("full_range", None),
          ("day1_-95_80", (-95.0, 80.0)),
          ("day2_82_110", (82.0, 110.0)),
          ("config_0_90", "config"))


# ------------------------------------------------------------- 2. гармоники
def harmonics(verbose=True):
    out = {}
    for tag, limit in RANGES:
        try:
            h = ml.analyse_dataset(DATASET, angles_limit=limit, verbose=False)
        except Exception as exc:                                   # noqa: BLE001
            out[tag] = {"error": f"{type(exc).__name__}: {exc}"}
            continue
        ps = h["parseval_order4"]
        out[tag] = {
            "eta_amplitude": ps["eta_amplitude"],
            "eta_per_freq_median": h["per_freq_summary"]["eta_median"],
            "theta0_deg": ps["theta0_deg_from_h2"],
            "theta0_err_deg": ps["theta0_err_deg"],
            "A4_over_A2": ps.get("A4_over_A2"),
            "n_angles": h.get("n_angles"),
        }
        if verbose:
            o = out[tag]
            print(f"  [{tag:<12}] eta = {o['eta_amplitude']:.5f} "
                  f"(медиана по частотам {o['eta_per_freq_median']:.5f}), "
                  f"theta0 = {o['theta0_deg']:+.3f} +- {o['theta0_err_deg']:.3f} град")
    return out


# ------------------------------------------------------------------ 3. фит
def fit(verbose=True, ranges=RANGES):
    import a3_eps_stability as a3
    res = {}
    for tag, limit in ranges:
        try:
            r = a3.analyse(DATASET, angles_limit=limit, verbose=verbose)
        except Exception as exc:                                   # noqa: BLE001
            res[tag] = {"error": f"{type(exc).__name__}: {exc}"}
            if verbose:
                print(f"  [{tag}] фит не выполнен: {type(exc).__name__}: {exc}")
            continue
        r.pop("starts", None)                    # 21 старт на диапазон — в артефакт не нужны
        res[tag] = r
    return res


# ------------------------------------------------------------------- main
def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--quality", action="store_true")
    ap.add_argument("--harmonics", action="store_true")
    ap.add_argument("--fit", action="store_true")
    ap.add_argument("--out", default=str(OUT / "a7_purewave.json"))
    args = ap.parse_args()
    do_all = not (args.quality or args.harmonics or args.fit)

    dm = DataManager(config.DATA_DIR)
    data = dm.get_data_for_dataset(DATASET)
    if not data:
        print(f"ОСТАНОВ: датасет {DATASET} не виден DataManager "
              f"(запусти normalize_purewave_data.py --apply)")
        return 1
    geo = fit_lib.GEOMETRY[DATASET]
    print(f"purewave: {len(data)} углов {min(data):.0f}..{max(data):.0f} град, "
          f"геометрия по микрофото P = {geo['P_um']} мкм, D_phys = {geo['D_phys_um']} мкм "
          f"(D/P = {geo['D_phys_um'] / geo['P_um']:.2f})")

    payload = {"task": "A7_purewave", "session": "A", "nick": "A3", "dataset": DATASET,
               "geometry": geo}
    if args.quality or do_all:
        payload["quality"] = quality()
    if args.harmonics or do_all:
        print("\n--- ГАРМОНИКИ A0 (модельно-независимо) ---")
        payload["harmonics"] = harmonics()
    if args.fit or do_all:
        print("\n--- ФИТ M5 (мультистарт по сетке A3) ---")
        payload["fit"] = fit()

    Path(args.out).write_text(json.dumps(payload, ensure_ascii=False, indent=2,
                                         default=float), encoding="utf-8")
    print(f"\nАртефакт: {Path(args.out).relative_to(REPO)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
