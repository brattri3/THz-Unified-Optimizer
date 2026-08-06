#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""A13 — ввод датасета `test_grid_33_11` в оборот: качество, два контроля, гармоники, фит.

ЗАЧЕМ. Восьмой образец проекта и **самый разреженный**: `P = 33 мкм`, `D = 11 мкм`,
`D/P = 0.333` против прежнего диапазона 0.416…0.710 по семи образцам. Поэтому это
первая ЭКСТРАПОЛЯЦИОННАЯ проверка эмпирического закона H5 `D_eff/D_phys = 1 − 0.85·D/P`,
а не ещё одна точка внутри интервала.

ПРЕДСКАЗАНИЯ, ЗАФИКСИРОВАННЫЕ ДО ИЗМЕРЕНИЯ (`data_pool/test_grid_33_11/README.md`,
коммит ad9818b от 2026-08-06, данные сняты позже в тот же день):
    D_eff/D_phys = 0.717 (D_eff ~ 7.9 мкм) | eta = 0.04 +- 0.02
    экстинкция 24…30 дБ                    | минимум не обязательно при 90 град
Скрипт печатает предсказание рядом с измеренным — без ручного сравнения постфактум.

ДВА КОНТРОЛЯ, ЗАЛОЖЕННЫЕ В ПОРЯДОК СЪЁМКИ (маршрутный лист §2)
--------------------------------------------------------------
Порядок блоков подобран так, чтобы дрейф и физику можно было РАЗДЕЛИТЬ:

  * ДРЕЙФ:   тень B rep1 (блок 1, начало сессии) против тени B rep2 (блок 4, конец).
             Один и тот же угол, разное время ⇒ различие = дрейф установки.
  * ФИЗИКА:  тень A (блок 3) против тени B rep2 (блок 4) — сняты рядом по времени.
             Углы theta и theta+180 суть одна ориентация проводов (pi-периодичность),
             поэтому различие = неоднородность образца, а не дрейф.

На `purewave` эти два эффекта были смешаны (два дня, переустановка) и разделить их
удалось только постфактум. Здесь разделение заложено в план.

ЧТО ДЕЛАЕТ СКРИПТ
  1. `--quality`   — дрейф фона по 17 трассам bg, динамический диапазон, положение
                     реального минимума (парабола по трём точкам), оба контроля выше.
  2. `--harmonics` — модельно-независимые eta, theta0, A4/A2 (линейный МНК, гармоники A0)
                     на полном диапазоне, на каждой половине и на продакшн-обрезке (0,90).
  3. `--fit`       — мультистартовый M5 через `a3_eps_stability.analyse` (21 старт).

Гардрейлы: `data_pool/` и `unified_optimizer/` — чтение; ядро (`fit_lib`, `model_core`) —
чтение. Геометрия нового образца ВРЕМЕННО инжектируется в `fit_lib.GEOMETRY` в памяти
процесса (файл зоны B не правится); постоянная запись — запросом к B.

Запуск ИЗ КОРНЯ репозитория:
    .venv/Scripts/python.exe research/two_wgp/a13_test_grid_33_11.py --quality --harmonics
    .venv/Scripts/python.exe research/two_wgp/a13_test_grid_33_11.py
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

import fit_lib                                                        # noqa: E402  зона B, ЧТЕНИЕ

DATASET = "test_grid_33_11"
P_UM, D_PHYS_UM = 33.0, 11.0

# --- ВРЕМЕННАЯ инжекция геометрии (постоянная запись — зона B, запрос в HANDOFFS).
# setdefault: если B уже внесла запись, она имеет приоритет.
fit_lib.GEOMETRY.setdefault(DATASET, {"P_um": P_UM, "D_phys_um": D_PHYS_UM})

import exp_builder                                                    # noqa: E402
import malus_linearization as ml                                      # noqa: E402
from unified_optimizer import config                                  # noqa: E402
from unified_optimizer.data_manager import DataManager                # noqa: E402
from unified_optimizer.optimizer_2d import get_transmission_spectra   # noqa: E402

OUT = REPO / "research" / "results" / "two_wgp"
OUT.mkdir(parents=True, exist_ok=True)

# Предсказания, зафиксированные ДО измерения (README датасета, коммит ad9818b).
PREREG = {
    "D_over_Dphys": 0.717,          # закон H5: 1 - 0.85 * 0.3333
    "eta": (0.02, 0.06),            # 0.04 +- 0.02
    "extinction_db": (24.0, 30.0),
    "theta_min_deg": "не обязательно 90",
}

SHADOW_B = [80, 82, 84, 86, 88, 90, 92, 94, 96, 98, 100]
SHADOW_A = [-80, -82, -84, -86, -88, -90, -92, -94, -96, -98, -100]
BRIGHT_ANCHORS = [40, 0, -40]

# Диапазоны для сравнения. ВАЖНО: каждый обязан содержать окно (0,90) — внутренняя
# служебная сверка `malus_linearization.analyse_dataset` жёстко считает подфит на
# `0 <= theta <= 90`, и на диапазоне без этого окна падает с «углов 1 < параметров 5».
#   full_range        — все 37 углов, то, ради чего снимался полный диапазон;
#   one_shadow_-70_100 — как если бы тень A НЕ снимали (цена отказа от контрольной группы);
#   config_0_90       — продакшн-обрезка, для сравнимости с числами A3 по другим образцам.
RANGES = (("full_range", None),
          ("one_shadow_-70_100", (-70.0, 100.0)),
          ("config_0_90", "config"))


# ------------------------------------------------------------------ утилиты
def load_raw():
    dm = DataManager(config.DATA_DIR)
    return {k[1:]: v for k, v in dm.raw_data.items() if k[0] == DATASET}


def trans_by_angle(rep_filter=None):
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
    fmax = getattr(config, "F_MAX_2D", config.F_MAX)
    return np.where((freqs >= config.F_MIN) & (freqs <= fmax))[0]


def _stat(rows, key):
    v = np.array([r[key] for r in rows], float)
    return {"n": len(v), "mean": float(v.mean()),
            "std": float(v.std(ddof=1)) if len(v) > 1 else 0.0,
            "min": float(v.min()), "max": float(v.max()),
            "spread_pct": float(100 * (v.max() - v.min()) / abs(v.mean())) if v.mean() else None}


# --------------------------------------------------------------- 1. качество
def quality(verbose=True):
    raw = load_raw()

    # --- дрейф фона: восстанавливаем уникальные трассы bg (они копировались) ---
    seen, bgs = {}, []
    for (angle, rep), d in sorted(raw.items()):
        if "bg" not in d:
            continue
        t_b, E_b = d["bg"]
        key = (float(E_b[0]), float(E_b[-1]), float(np.sum(E_b)))
        if key in seen:
            continue
        seen[key] = True
        i_pk = int(np.argmax(np.abs(E_b)))
        bgs.append({"first_angle": float(angle), "first_rep": int(rep),
                    "peak_abs": float(np.abs(E_b[i_pk])),
                    "t_peak_ps": float(t_b[i_pk]),
                    "rms": float(np.sqrt(np.mean(E_b ** 2)))})
    bg_drift = {"n_unique_backgrounds": len(bgs),
                "peak_abs": _stat(bgs, "peak_abs"),
                "t_peak_ps": _stat(bgs, "t_peak_ps"),
                "per_background": bgs}

    # --- динамический диапазон и положение тени ---
    freqs, tr1 = trans_by_angle(rep_filter=1)
    idx = band(freqs)
    per_angle = []
    for (angle, rep), t in sorted(tr1.items()):
        u = np.abs(t[idx]) ** 2
        per_angle.append({"angle": float(angle), "U_parseval": float(u.sum()),
                          "amp_median": float(np.median(np.abs(t[idx])))})
    U = np.array([p["U_parseval"] for p in per_angle])
    A = np.array([p["angle"] for p in per_angle])
    i_max, i_min = int(np.argmax(U)), int(np.argmin(U))
    theta_min_parab = None
    if 0 < i_min < len(U) - 1:
        x, y = A[i_min - 1:i_min + 2], np.log(U[i_min - 1:i_min + 2])
        c = np.polyfit(x, y, 2)
        if c[0] > 0:
            theta_min_parab = float(-c[1] / (2 * c[0]))
    dyn = {"U_max": float(U[i_max]), "angle_U_max": float(A[i_max]),
           "U_min": float(U[i_min]), "angle_U_min": float(A[i_min]),
           "extinction_ratio": float(U[i_max] / U[i_min]),
           "extinction_db": float(10 * np.log10(U[i_max] / U[i_min])),
           "angle_U_min_parabolic": theta_min_parab,
           "prereg_extinction_db": PREREG["extinction_db"],
           "per_angle": per_angle}

    # --- КОНТРОЛЬ 1: дрейф (один угол, разное время) ---
    _f, tr2 = trans_by_angle(rep_filter=2)
    drift = []
    for a in SHADOW_B + BRIGHT_ANCHORS:
        k1, k2 = (float(a), 1), (float(a), 2)
        if k1 not in tr1 or k2 not in tr2:
            continue
        v1, v2 = tr1[k1][idx], tr2[k2][idx]
        u1, u2 = float((np.abs(v1) ** 2).sum()), float((np.abs(v2) ** 2).sum())
        drift.append({"angle": float(a), "zone": "shadow" if a in SHADOW_B else "bright",
                      "U_rep1": u1, "U_rep2": u2,
                      "ratio_db": float(10 * np.log10(u1 / u2)),
                      "amp_ratio_median": float(np.median(np.abs(v1) / np.abs(v2)))})
    d_sh = [d["ratio_db"] for d in drift if d["zone"] == "shadow"]
    d_br = [d["ratio_db"] for d in drift if d["zone"] == "bright"]
    drift_sum = {
        "per_angle": drift,
        "shadow_mean_db": float(np.mean(d_sh)) if d_sh else None,
        "shadow_std_db": float(np.std(d_sh, ddof=1)) if len(d_sh) > 1 else None,
        "bright_mean_db": float(np.mean(d_br)) if d_br else None,
        "bright_std_db": float(np.std(d_br, ddof=1)) if len(d_br) > 1 else None,
        "note": "rep1 снят в начале сессии, rep2 — в конце; один угол ⇒ чистый дрейф",
    }

    # --- КОНТРОЛЬ 2: pi-периодичность (разные участки, совпадающее время) ---
    # theta и theta-180 — одна ориентация проводов. Тень A (блок 3) и тень B rep2
    # (блок 4) сняты подряд, поэтому дрейф между ними минимален.
    pi_pairs = []
    for a in SHADOW_B:
        a_neg = float(a - 180)                    # 80 -> -100, 100 -> -80
        k_neg, k_pos2 = (a_neg, 1), (float(a), 2)
        if k_neg not in tr1 or k_pos2 not in tr2:
            continue
        v_neg, v_pos = tr1[k_neg][idx], tr2[k_pos2][idx]
        u_neg = float((np.abs(v_neg) ** 2).sum())
        u_pos = float((np.abs(v_pos) ** 2).sum())
        pi_pairs.append({"angle_neg": a_neg, "angle_pos": float(a),
                         "U_neg": u_neg, "U_pos": u_pos,
                         "ratio_db": float(10 * np.log10(u_neg / u_pos)),
                         "amp_ratio_median": float(np.median(np.abs(v_neg) / np.abs(v_pos)))})
    pr = [p["ratio_db"] for p in pi_pairs]
    pi_sum = {"per_pair": pi_pairs,
              "mean_db": float(np.mean(pr)) if pr else None,
              "std_db": float(np.std(pr, ddof=1)) if len(pr) > 1 else None,
              "max_abs_db": float(np.max(np.abs(pr))) if pr else None,
              "note": "тень A (блок 3) против тени B rep2 (блок 4): сняты рядом по времени, "
                      "разные участки решётки ⇒ мера неоднородности образца"}

    out = {"background_drift": bg_drift, "dynamic_range": dyn,
           "control_drift": drift_sum, "control_pi_periodicity": pi_sum}

    if verbose:
        print("\n--- КАЧЕСТВО ДАННЫХ test_grid_33_11 ---")
        print(f"  уникальных фонов: {bg_drift['n_unique_backgrounds']}; "
              f"пик |E| = {bg_drift['peak_abs']['mean']:.4f} "
              f"(разброс {bg_drift['peak_abs']['spread_pct']:.1f}%), "
              f"t_пика разброс {bg_drift['t_peak_ps']['max'] - bg_drift['t_peak_ps']['min']:.3f} пс")
        print(f"  максимум при {dyn['angle_U_max']:+.0f}°, минимум при {dyn['angle_U_min']:+.0f}° "
              f"(парабола {dyn['angle_U_min_parabolic']:+.2f}°), "
              f"экстинкция {dyn['extinction_db']:.2f} дБ "
              f"[предсказано {PREREG['extinction_db'][0]:.0f}…{PREREG['extinction_db'][1]:.0f}]")
        print("  КОНТРОЛЬ 1 — дрейф (один угол, начало против конца сессии):")
        if d_sh:
            print(f"    тень:  {drift_sum['shadow_mean_db']:+.3f} дБ "
                  f"± {drift_sum['shadow_std_db']:.3f} ({len(d_sh)} углов)")
        if d_br:
            print(f"    ярко:  {drift_sum['bright_mean_db']:+.3f} дБ "
                  f"± {drift_sum['bright_std_db']:.3f} ({len(d_br)} углов)")
        print("  КОНТРОЛЬ 2 — π-периодичность (разные участки, совпадающее время):")
        if pr:
            print(f"    среднее {pi_sum['mean_db']:+.3f} дБ ± {pi_sum['std_db']:.3f}, "
                  f"максимум по модулю {pi_sum['max_abs_db']:.3f} дБ ({len(pr)} пар)")
    return out


# ------------------------------------------------------------- 2. гармоники
def harmonics(verbose=True):
    out = {}
    for tag, limit in RANGES:
        try:
            h = ml.analyse_dataset(DATASET, angles_limit=limit, verbose=False)
        except Exception as exc:                                      # noqa: BLE001
            out[tag] = {"error": f"{type(exc).__name__}: {exc}"}
            if verbose:
                print(f"  [{tag:<16}] не выполнено: {type(exc).__name__}: {exc}")
            continue
        ps = h["parseval_order4"]
        out[tag] = {"eta_amplitude": ps["eta_amplitude"],
                    "eta_per_freq_median": h["per_freq_summary"]["eta_median"],
                    "theta0_deg": ps["theta0_deg_from_h2"],
                    "theta0_err_deg": ps["theta0_err_deg"],
                    "A4_over_A2": ps.get("harmonic4_over_harmonic2"),
                    "n_angles": h.get("n_angles")}
        if verbose:
            o = out[tag]
            a42 = f", A4/A2 = {o['A4_over_A2']:.4f}" if o.get("A4_over_A2") is not None else ""
            print(f"  [{tag:<18}] n={o['n_angles']:>2}  eta = {o['eta_amplitude']:.5f} "
                  f"(медиана по частотам {o['eta_per_freq_median']:.5f}), "
                  f"theta0 = {o['theta0_deg']:+.3f} ± {o['theta0_err_deg']:.3f}°{a42}")
    out["prereg_eta_range"] = PREREG["eta"]
    return out


# ------------------------------------------------------------------ 3. фит
def fit(verbose=True, ranges=RANGES):
    import a3_eps_stability as a3
    res = {}
    for tag, limit in ranges:
        try:
            r = a3.analyse(DATASET, angles_limit=limit, verbose=verbose)
        except Exception as exc:                                      # noqa: BLE001
            res[tag] = {"error": f"{type(exc).__name__}: {exc}"}
            if verbose:
                print(f"  [{tag}] фит не выполнен: {type(exc).__name__}: {exc}")
            continue
        r.pop("starts", None)
        res[tag] = r
    res["prereg_D_over_Dphys"] = PREREG["D_over_Dphys"]
    return res


# ------------------------------------------------------------------- main
def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--quality", action="store_true")
    ap.add_argument("--harmonics", action="store_true")
    ap.add_argument("--fit", action="store_true")
    ap.add_argument("--out", default=str(OUT / "a13_test_grid_33_11.json"))
    args = ap.parse_args()
    do_all = not (args.quality or args.harmonics or args.fit)

    dm = DataManager(config.DATA_DIR)
    if not dm.get_data_for_dataset(DATASET):
        print(f"ОСТАНОВ: датасет {DATASET} не виден DataManager "
              f"(запусти normalize_test_grid_33_11.py --apply)")
        return 2

    payload = {"dataset": DATASET, "P_um": P_UM, "D_phys_um": D_PHYS_UM,
               "D_over_P": round(D_PHYS_UM / P_UM, 4), "prereg": PREREG}
    if args.quality or do_all:
        payload["quality"] = quality()
    if args.harmonics or do_all:
        print("\n--- ГАРМОНИКИ A0 (модельно-независимо) ---")
        payload["harmonics"] = harmonics()
    if args.fit or do_all:
        print("\n--- МУЛЬТИСТАРТОВЫЙ ФИТ M5 ---")
        payload["fit"] = fit()

    Path(args.out).write_text(json.dumps(payload, ensure_ascii=False, indent=2,
                                         default=float), encoding="utf-8")
    print(f"\nартефакт: {args.out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
