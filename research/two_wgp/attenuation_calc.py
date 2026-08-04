#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
A4 — калькулятор затухания двух-WGP аттенюатора: A(theta1, theta2, nu).

НАЗНАЧЕНИЕ (разграничение зафиксировано хэндоффом C->A от 2026-07-28)
---------------------------------------------------------------------
Это **физический эталон**, а не продуктовый расчёт. Критерий качества —
точность и идентифицируемость; потребители — сама Сессия A (проверка гипотез)
и Сессия D (статья). Продуктовый расчёт для оператора прибора живёт отдельно в
`attenuator_app/core/` (критерий там другой: воспроизводимость и консервативный
интервал). Две реализации существуют СОЗНАТЕЛЬНО; сверка между ними — задача
`C2_a4_crosscheck`, и режим `--grid-for-C` этого файла готовит для неё таблицу.

ЧТО СЧИТАЕТСЯ
-------------
Поле на детекторе (модель A1, матрицы Джонса, per-element геометрия):

    E(θ1, θ2, ν) = d^T · J_WGP2(θ2, ν) · J_WGP1(θ1, ν) · E_in

Затухание — ОТНОСИТЕЛЬНОЕ, к согласованной конфигурации (оба элемента на 0):

    A(θ1, θ2, ν) = −10·lg [ |E(θ1,θ2,ν)|² / |E(0,0,ν)|² ]      [дБ]

Относительная нормировка выбрана НЕ произвольно: она совпадает с определением,
по которому Сессия C валидировала прибор на 69 точках, и тождественно сокращает
общие множители (потери ℓ(ν), |t⊥|⁴ и т.п.). Поэтому число сравнимо между
реализациями и с измерением.

ОБОЗНАЧЕНИЯ
-----------
θ1, θ2   углы осей ПРОПУСКАНИЯ элементов (⊥ проводам), град, CCW+
d        ось чувствительности детектора; E_in — поляризация излучателя
ν        частота, ТГц; P, D — период и эффективный диаметр проволоки, мкм

Зона A. Общее ядро Бланко зовётся через `model_2wgp.blanco_t` -> `model_core`.

Запуск:
    .venv\\Scripts\\python.exe research/two_wgp/attenuation_calc.py --selftest
    ... --map --p 15.5 --d 11.0
    ... --spectrum --theta1 45 --theta2 0
    ... --grid-for-C
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO))
sys.path.insert(0, str(REPO / "research" / "two_wgp"))

import model_2wgp as m2                                    # noqa: E402
from unified_optimizer import config                       # noqa: E402

RESULTS = REPO / "research" / "results" / "two_wgp"

# Геометрия образцов — единственный источник истины на сегодня (fit_lib.GEOMETRY,
# read-only). Дублируется здесь ТОЛЬКО как значения по умолчанию для CLI.
DEFAULT_GEOM = {"att-11-16-356": (15.5, 11.0), "series": (16.0, 11.0),
                "specac": (24.9, 14.0), "test_grid_40_20": (38.8, 18.0),
                "purewave": (25.5, 10.6)}

_EPS = 1e-300


# --------------------------------------------------------------- сборка модели
def make_wgp(p_um, d_um, **kw):
    """Один элемент. kw — необязательные надстройки (eta0, delta0, tau_leak_ps…)."""
    return m2.WGP(p_um=p_um, d_um=d_um, **kw)


def _field(theta1, theta2, freqs, wgp1, wgp2, **kw):
    return m2.forward_field(theta1, theta2, freqs, wgp1, wgp2, **kw)


# ------------------------------------------------------------------ основное
def attenuation_db(theta1_deg, theta2_deg, freqs_thz, wgp1, wgp2=None,
                   reference=(0.0, 0.0), absolute=False, **kw):
    """Затухание в дБ. Основная точка входа.

    theta1_deg, theta2_deg  скаляры или массивы (broadcast между собой):
        - развёртка прибора:  theta1 = (n,),  theta2 = скаляр
        - карта:              theta1 = (n1,1), theta2 = (1,n2)  -> (n1,n2,nf)
    reference  углы (θ1,θ2) согласованной конфигурации; None вместе с
        absolute=True даёт АБСОЛЮТНОЕ затухание (−10 lg |E|², без нормировки).
    Возвращает ndarray формы broadcast(θ1,θ2) + (nf,), в дБ, >= 0 для reference.
    """
    freqs = np.asarray(freqs_thz, float)
    P = np.abs(_field(theta1_deg, theta2_deg, freqs, wgp1, wgp2, **kw)) ** 2
    if absolute or reference is None:
        return -10.0 * np.log10(np.maximum(P, _EPS))
    P0 = np.abs(_field(reference[0], reference[1], freqs, wgp1, wgp2, **kw)) ** 2
    return -10.0 * np.log10(np.maximum(P, _EPS) / np.maximum(P0, _EPS))


def band_attenuation_db(theta1_deg, theta2_deg, freqs_thz, wgp1, wgp2=None,
                        reference=(0.0, 0.0), **kw):
    """Одно число на конфигурацию: затухание ПОЛНОЙ ЭНЕРГИИ в полосе.

    A = −10·lg [ Σ_ν |E(θ1,θ2,ν)|² / Σ_ν |E(ref,ν)|² ]  (теорема Парсеваля).
    Именно эта величина измеряется прибором и валидировалась Сессией C — не
    среднее по частоте от A(ν), которое было бы другим числом (среднее
    логарифма ≠ логарифм среднего).
    """
    freqs = np.asarray(freqs_thz, float)
    U = np.sum(np.abs(_field(theta1_deg, theta2_deg, freqs, wgp1, wgp2, **kw)) ** 2,
               axis=-1)
    U0 = np.sum(np.abs(_field(reference[0], reference[1], freqs,
                              wgp1, wgp2, **kw)) ** 2, axis=-1)
    return -10.0 * np.log10(np.maximum(U, _EPS) / np.maximum(U0, _EPS))


def attenuation_map(freqs_thz, wgp1, wgp2=None, theta1_grid=None,
                    theta2_grid=None, **kw):
    """Карта A(θ1, θ2) — затухание энергии в полосе. Возвращает (θ1, θ2, A_дБ)."""
    t1 = np.arange(0.0, 90.5, 2.0) if theta1_grid is None else np.asarray(theta1_grid, float)
    t2 = np.array([0.0]) if theta2_grid is None else np.asarray(theta2_grid, float)
    A = band_attenuation_db(t1[:, None], t2[None, :], freqs_thz, wgp1, wgp2, **kw)
    return t1, t2, A


def required_angle(target_db, freqs_thz, wgp1, wgp2=None, theta2_deg=0.0,
                   bracket=(0.0, 90.0), tol=1e-6, **kw):
    """ОБРАТНАЯ задача: угол θ1, дающий заданное затухание (бисекция).

    A(θ1) монотонно возрастает на [0,90] у любого разумного поляризатора, но
    НЕ строго — вблизи 90° кривая выполаживается на «полке» утечки. Поэтому
    если цель выше достижимого максимума, возвращается None вместе с этим
    максимумом: молча упереться в границу — худший из возможных ответов.
    """
    lo, hi = bracket
    f_lo = band_attenuation_db(lo, theta2_deg, freqs_thz, wgp1, wgp2, **kw)
    f_hi = band_attenuation_db(hi, theta2_deg, freqs_thz, wgp1, wgp2, **kw)
    if not (f_lo <= target_db <= f_hi):
        return {"theta1_deg": None, "achievable_min_db": float(f_lo),
                "achievable_max_db": float(f_hi),
                "reason": "цель вне достижимого диапазона затухания"}
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if band_attenuation_db(mid, theta2_deg, freqs_thz, wgp1, wgp2, **kw) < target_db:
            lo = mid
        else:
            hi = mid
        if hi - lo < tol:
            break
    th = 0.5 * (lo + hi)
    return {"theta1_deg": float(th),
            "attained_db": float(band_attenuation_db(th, theta2_deg, freqs_thz,
                                                     wgp1, wgp2, **kw)),
            "achievable_max_db": float(f_hi)}


# ------------------------------------------------------------------ проверки
def selftest(verbose=True):
    """Валидация: идеальный предел, симметрии, согласие с A1, обратная задача."""
    checks = []

    def chk(name, ok, detail):
        checks.append({"check": name, "ok": bool(ok), "detail": detail})
        if verbose:
            print(f"  [{'OK ' if ok else 'FAIL'}] {name}: {detail}")

    freqs = np.arange(0.2, 1.5001, 0.01)
    ideal = m2.WGP(t_override=(1.0, 0.0))
    real = m2.WGP(p_um=15.5, d_um=11.0)

    # 1. Идеальный предел = закон Малюса для ДВУХ поляризаторов.
    #    E = cos(t1)*cos(t1-t2)*cos(t2)  (E_in || d || x)
    t1 = np.array([0.0, 15.0, 30.0, 45.0, 60.0, 75.0])
    t2 = np.array([0.0, 10.0, 25.0])
    A = attenuation_db(t1[:, None], t2[None, :], freqs, ideal, ideal)
    c = (np.cos(np.deg2rad(t1))[:, None] * np.cos(np.deg2rad(t1[:, None] - t2[None, :]))
         * np.cos(np.deg2rad(t2))[None, :])
    A_ref = -20.0 * np.log10(np.abs(c))
    dev = float(np.max(np.abs(A - A_ref[..., None])))
    chk("идеальный предел = закон Малюса для двух элементов", dev < 1e-12,
        f"max|отклонение| = {dev:.2e} дБ")

    # 2. Затухание в опорной точке тождественно нулевое
    z = float(np.max(np.abs(attenuation_db(0.0, 0.0, freqs, real, real))))
    chk("A(0,0) == 0 по определению нормировки", z < 1e-12, f"max|A| = {z:.2e} дБ")

    # 3. Чётность по θ1 у РЕАЛЬНЫХ элементов (нет нечётных членов в модели)
    a_p = attenuation_db(37.0, 0.0, freqs, real, real)
    a_m = attenuation_db(-37.0, 0.0, freqs, real, real)
    dev3 = float(np.max(np.abs(a_p - a_m)))
    chk("A(+θ1) == A(−θ1): модель ЧЁТНАЯ (ср. вердикт HN12 на данных)",
        dev3 < 1e-12, f"max|разность| = {dev3:.2e} дБ")

    # 4. Согласие с A1: тот же |E|^2, что даёт model_2wgp.attenuation
    P = m2.attenuation(40.0, 0.0, freqs, real, real)
    P0 = m2.attenuation(0.0, 0.0, freqs, real, real)
    dev4 = float(np.max(np.abs(attenuation_db(40.0, 0.0, freqs, real, real)
                               - (-10 * np.log10(P / P0)))))
    chk("бит-совпадение с ядром A1 (model_2wgp.attenuation)", dev4 == 0.0,
        f"max|разность| = {dev4:.2e} дБ")

    # 5. Монотонность по θ1 на [0,85]
    th = np.arange(0.0, 85.1, 2.5)
    Ab = band_attenuation_db(th, 0.0, freqs, real, real)
    chk("A(θ1) монотонно возрастает на [0,85]", bool(np.all(np.diff(Ab) > 0)),
        f"мин. приращение = {float(np.min(np.diff(Ab))):+.4f} дБ")

    # 6. Затухание энергии в полосе <= взвешенного среднего A(nu) — НЕРАВЕНСТВО
    #    ЙЕНСЕНА (логарифм вогнут). Знакоопределённо, порог не нужен: считать
    #    надо энергию, «среднее затухание по частоте» её систематически завышает.
    ok6, det6 = True, []
    for th_i in (30.0, 60.0, 85.0):
        Pn = np.abs(_field(th_i, 0.0, freqs, real, real)) ** 2
        P0n = np.abs(_field(0.0, 0.0, freqs, real, real)) ** 2
        w6 = P0n / P0n.sum()
        a_band = float(band_attenuation_db(th_i, 0.0, freqs, real, real))
        a_wmean = float(-10.0 * np.sum(w6 * np.log10(Pn / P0n)))
        ok6 &= a_band <= a_wmean + 1e-12
        det6.append(f"{th_i:.0f}°: полоса {a_band:.3f} <= среднее {a_wmean:.3f}")
    chk("затухание энергии <= взвеш. среднего A(ν) (Йенсен): считать надо энергию",
        ok6, "; ".join(det6))

    # 7. Обратная задача воспроизводит прямую
    tgt = 20.0
    inv = required_angle(tgt, freqs, real, real)
    chk("обратная задача: A(θ1(20 дБ)) == 20 дБ",
        inv["theta1_deg"] is not None and abs(inv["attained_db"] - tgt) < 1e-4,
        f"θ1 = {inv['theta1_deg']:.4f}°, получено {inv['attained_db']:.6f} дБ")

    # 8. Недостижимая цель отвергается явно, а не молча упирается в границу
    inv2 = required_angle(1e3, freqs, real, real)
    chk("недостижимая цель -> явный отказ с указанием максимума",
        inv2["theta1_deg"] is None,
        f"максимум {inv2.get('achievable_max_db', float('nan')):.2f} дБ")

    # 9. Утечка ограничивает динамический диапазон (физическая осмысленность)
    leaky = m2.WGP(p_um=15.5, d_um=11.0, eta0=0.03)
    dr_clean = float(band_attenuation_db(90.0, 0.0, freqs, real, real))
    dr_leaky = float(band_attenuation_db(90.0, 0.0, freqs, leaky, leaky))
    chk("добавленная утечка eta0=0.03 УМЕНЬШАЕТ предельное затухание",
        dr_leaky < dr_clean,
        f"{dr_clean:.2f} дБ -> {dr_leaky:.2f} дБ (потеря {dr_clean - dr_leaky:.2f} дБ)")

    n_ok = sum(c["ok"] for c in checks)
    if verbose:
        print(f"\nИТОГ: {n_ok}/{len(checks)}")
    return {"checks": checks, "n_ok": n_ok, "n_total": len(checks),
            "ALL_PASS": n_ok == len(checks)}


# ----------------------------------------------------- таблица для сверки с C
def grid_for_C(p_um=15.5, d_um=11.0, verbose=True):
    """Таблица A(θ, ν) на сетке 0.2–1.5 ТГц x 0–90° — критерий готовности C->A.

    Формат сознательно плоский и самодостаточный: C читает json и сравнивает со
    своим `attenuator_app/core`, ничего не импортируя из зоны A.
    """
    freqs = np.round(np.arange(0.2, 1.5001, 0.05), 4)
    angles = np.arange(0.0, 90.1, 5.0)
    w = m2.WGP(p_um=p_um, d_um=d_um)
    A = attenuation_db(angles[:, None], 0.0, freqs, w, w)[:, 0, :]
    band = band_attenuation_db(angles, 0.0, freqs, w, w)
    out = {
        "producer": "A4 research/two_wgp/attenuation_calc.py (физический эталон)",
        "definition": "A = -10*lg[ |E(theta,nu)|^2 / |E(0,nu)|^2 ]; "
                      "band = -10*lg[ sum_nu|E(theta)|^2 / sum_nu|E(0)|^2 ]",
        "model": "два реальных WGP, Джонс, per-element геометрия; "
                 "WGP2 == WGP1, theta2 = 0, E_in = d = 0 град",
        "geometry": {"p_um": p_um, "d_um": d_um, "N": 15, "drude": True},
        "caveat": "d_um здесь — ЭФФЕКТИВНЫЙ диаметр (загадка D_eff не решена); "
                  "при сверке брать одинаковый d_um в обеих реализациях, иначе "
                  "сравниваются разные образцы, а не разные калькуляторы.",
        "freqs_thz": freqs.tolist(),
        "angles_deg": angles.tolist(),
        "A_db": [[float(v) for v in row] for row in A],
        "A_band_db": [float(v) for v in band],
        "alarm_threshold_db": 1.39,
        "alarm_threshold_note": "паспортная RMSE прибора; расхождение выше — ошибка "
                                "в одной из реализаций (порог задан Сессией C)",
    }
    if verbose:
        print(f"  сетка: {len(angles)} углов x {len(freqs)} частот, P={p_um} D={d_um}")
        print("  A_band(θ), дБ: " +
              " ".join(f"{a:.0f}:{v:.2f}" for a, v in zip(angles, band)))
    return out


# ---------------------------------------------------------------------- CLI
def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--selftest", action="store_true")
    ap.add_argument("--map", action="store_true", help="карта A(θ1,θ2)")
    ap.add_argument("--spectrum", action="store_true", help="спектр A(ν)")
    ap.add_argument("--grid-for-C", action="store_true",
                    help="таблица для сверки C2_a4_crosscheck")
    ap.add_argument("--target-db", type=float, help="обратная задача: нужный уровень")
    ap.add_argument("--p", type=float, default=15.5, help="период, мкм")
    ap.add_argument("--d", type=float, default=11.0, help="эфф. диаметр, мкм")
    ap.add_argument("--theta1", type=float, default=45.0)
    ap.add_argument("--theta2", type=float, default=0.0)
    ap.add_argument("--out", default=None)
    args = ap.parse_args()

    freqs = np.arange(config.F_MIN, getattr(config, "F_MAX_2D", config.F_MAX)
                      + 1e-12, 0.01)
    w = m2.WGP(p_um=args.p, d_um=args.d)
    RESULTS.mkdir(parents=True, exist_ok=True)

    if args.selftest:
        r = selftest()
        p = Path(args.out or RESULTS / "a4_validation.json")
        p.write_text(json.dumps(r, ensure_ascii=False, indent=2), encoding="utf-8")
        print(f"артефакт -> {p}")
        raise SystemExit(0 if r["ALL_PASS"] else 1)

    if args.grid_for_C:
        r = grid_for_C(args.p, args.d)
        p = Path(args.out or RESULTS / "a4_grid_for_crosscheck.json")
        p.write_text(json.dumps(r, ensure_ascii=False, indent=2), encoding="utf-8")
        print(f"артефакт -> {p}")
        raise SystemExit(0)

    if args.map:
        t1, t2, A = attenuation_map(freqs, w, w,
                                    theta2_grid=np.arange(0.0, 90.1, 5.0))
        print(f"  карта A(θ1,θ2), дБ: {A.shape[0]}x{A.shape[1]}, "
              f"диапазон {A.min():.2f}…{A.max():.2f}")
        p = Path(args.out or RESULTS / "a4_map.json")
        p.write_text(json.dumps({"theta1_deg": t1.tolist(), "theta2_deg": t2.tolist(),
                                 "A_band_db": A.tolist(),
                                 "geometry": {"p_um": args.p, "d_um": args.d}},
                                ensure_ascii=False, indent=2), encoding="utf-8")
        print(f"артефакт -> {p}")
        raise SystemExit(0)

    if args.spectrum:
        A = attenuation_db(args.theta1, args.theta2, freqs, w, w)
        band = float(band_attenuation_db(args.theta1, args.theta2, freqs, w, w))
        print(f"  θ1={args.theta1}°, θ2={args.theta2}°: A(ν) от {A.min():.2f} "
              f"до {A.max():.2f} дБ; энергия в полосе {band:.3f} дБ")
        p = Path(args.out or RESULTS / "a4_spectrum.json")
        p.write_text(json.dumps({"freqs_thz": freqs.tolist(), "A_db": A.tolist(),
                                 "A_band_db": band, "theta1_deg": args.theta1,
                                 "theta2_deg": args.theta2,
                                 "geometry": {"p_um": args.p, "d_um": args.d}},
                                ensure_ascii=False, indent=2), encoding="utf-8")
        print(f"артефакт -> {p}")
        raise SystemExit(0)

    if args.target_db is not None:
        print(f"  обратная задача: {required_angle(args.target_db, freqs, w, w)}")
        raise SystemExit(0)

    ap.print_help()


if __name__ == "__main__":
    main()
