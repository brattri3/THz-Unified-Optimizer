#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
A0_malus_linearization — гармоническая (линейная) форма закона Малюса
для НЕидеального поляризатора с утечкой.

Идея одной фразой
-----------------
Наблюдаемая мощность U(theta) нелинейна по углу, но ЛИНЕЙНА по коэффициентам,
если разложить её по гармоникам угла. Тогда подгонка — обычный линейный МНК
(замкнутая формула, единственный минимум, аналитическая ковариация), а не
нелинейная оптимизация с мультистартом.

Два уровня строгости
--------------------
order=2  «учебный» (интенсивностная / некогерентная композиция):
         U = a0 + a1*cos2t + b1*sin2t                        (3 параметра)
         эквивалентно U = U_min + (U_max-U_min)*cos^2(t-t0).
order=4  ТОЧНЫЙ для когерентной модели этого проекта
         (E = t_perp*cos^2(t) + t_par*sin^2(t), затем U = |E|^2):
         U = c0 + a2*cos2t + b2*sin2t + a4*cos4t + b4*sin4t   (5 параметров)

order=4 — не приближение: |E|^2 при E, линейном по {1, cos2t, sin2t},
ТОЧНО является тригонометрическим многочленом 4-го порядка. Всё, что угловое
сканирование может сказать на ОДНОЙ частоте, — это 5 чисел.

Обозначения
-----------
theta   угол поворота WGP, град (конвенция репозитория: ось ПРОПУСКАНИЯ)
t_perp  комплексный коэффициент пропускания «яркого» канала (E перпенд. проводам)
t_par   то же для «тёмного» канала (утечка)
eta     амплитудная утечка |t_par|/|t_perp|
dphi    относительная фаза arg(t_perp) - arg(t_par)  (в M5: delta0 + 2*pi*nu*tau_leak)
theta0  угловой офсет (нуль шкалы поворотного столика)

Зона A. Общий код (fit_lib/model_core) НЕ трогается — здесь только numpy.

Запуск:  .venv\\Scripts\\python.exe research/two_wgp/malus_linearization.py --selftest
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[2]
RESULTS = REPO / "research" / "results" / "two_wgp"

D2R = np.pi / 180.0
R2D = 180.0 / np.pi


# --------------------------------------------------------------- базис
def harmonic_design(theta_deg, order=4):
    """Матрица плана X: столбцы [1, cos2t, sin2t, (cos4t, sin4t)].

    Модель U = X @ p строго ЛИНЕЙНА по p — отсюда весь замкнутый аппарат МНК.
    """
    t = np.asarray(theta_deg, float) * D2R
    cols = [np.ones_like(t), np.cos(2 * t), np.sin(2 * t)]
    if order >= 4:
        cols += [np.cos(4 * t), np.sin(4 * t)]
    return np.column_stack(cols)


def fit_harmonics(theta_deg, U, sigma=None, order=4):
    """Взвешенный линейный МНК. Возвращает (p, cov, info).

    cov = (X^T W X)^{-1} — ТОЧНАЯ ковариация (модель линейна), бутстрап не нужен.
    Если sigma не задана, она оценивается из невязки: sigma^2 = RSS/(n-k).
    """
    X = harmonic_design(theta_deg, order)
    y = np.asarray(U, float)
    n, k = X.shape
    if n < k:
        raise ValueError(f"углов {n} < параметров {k}: система недоопределена")

    if sigma is None:
        w = np.ones(n)
    else:
        s = np.asarray(sigma, float)
        if np.any(s <= 0):
            raise ValueError("sigma должна быть > 0")
        w = 1.0 / s ** 2

    XtWX = X.T @ (w[:, None] * X)
    XtWy = X.T @ (w * y)
    p = np.linalg.solve(XtWX, XtWy)
    resid = y - X @ p
    dof = n - k
    rss = float(np.sum(w * resid ** 2))

    if sigma is None:
        s2 = rss / dof if dof > 0 else np.nan          # оценка дисперсии из невязки
        cov = s2 * np.linalg.inv(XtWX)
        chi2_red = np.nan
    else:
        cov = np.linalg.inv(XtWX)
        chi2_red = rss / dof if dof > 0 else np.nan

    info = {
        "n_angles": int(n), "n_params": int(k), "dof": int(dof),
        "rss": rss, "chi2_red": float(chi2_red),
        "resid_rms": float(np.sqrt(np.mean(resid ** 2))),
        "cond_XtWX": float(np.linalg.cond(XtWX)),
        "resid": resid,
    }
    return p, cov, info


# ------------------------------------------------- замкнутые формулы (order=2)
def malus_params(p, cov=None):
    """Учебный случай U = a0 + a1*cos2t + b1*sin2t -> контраст, утечка, офсет.

    theta0  = 1/2 * atan2(b1, a1)                      (град, mod 90)
    U_max   = a0 + A,  U_min = a0 - A,  A = hypot(a1,b1)
    ER      = U_max/U_min                              (степень экстинкции)
    eta_pow = U_min/U_max                              (утечка по МОЩНОСТИ)
    eta     = sqrt(eta_pow)                            (утечка по АМПЛИТУДЕ)
    V       = A/a0                                     (глубина модуляции)
    Ошибки — аналитические (дельта-метод), без бутстрапа.
    """
    a0, a1, b1 = float(p[0]), float(p[1]), float(p[2])
    A = float(np.hypot(a1, b1))
    t0 = 0.5 * np.arctan2(b1, a1) * R2D
    umax, umin = a0 + A, a0 - A
    out = {
        "theta0_deg": t0,
        "U_max": umax, "U_min": umin,
        "amplitude_A": A, "visibility": A / a0 if a0 else np.nan,
        "extinction_ratio": umax / umin if umin > 0 else np.inf,
        "extinction_db": 10 * np.log10(umax / umin) if umin > 0 else np.inf,
        "eta_power": umin / umax if umax else np.nan,
        "eta_amplitude": float(np.sqrt(max(umin / umax, 0.0))) if umax else np.nan,
    }
    if cov is not None and A > 0:
        v_a1, v_b1, c_ab = cov[1, 1], cov[2, 2], cov[1, 2]
        # theta0 = 1/2 atan2(b1,a1):  d/da1 = -b1/(2A^2), d/db1 = a1/(2A^2)
        var_t0 = (b1 ** 2 * v_a1 + a1 ** 2 * v_b1 - 2 * a1 * b1 * c_ab) / (4 * A ** 4)
        # A = hypot: dA/da1 = a1/A, dA/db1 = b1/A
        var_A = (a1 ** 2 * v_a1 + b1 ** 2 * v_b1 + 2 * a1 * b1 * c_ab) / A ** 2
        out["theta0_err_deg"] = float(np.sqrt(max(var_t0, 0.0)) * R2D)
        out["amplitude_A_err"] = float(np.sqrt(max(var_A, 0.0)))
        v_a0 = cov[0, 0]
        ca0A = (a1 * cov[0, 1] + b1 * cov[0, 2]) / A
        out["U_max_err"] = float(np.sqrt(max(v_a0 + var_A + 2 * ca0A, 0.0)))
        out["U_min_err"] = float(np.sqrt(max(v_a0 + var_A - 2 * ca0A, 0.0)))
    return out


# --------------------------------------- замкнутые формулы (order=4, когерентно)
def theta0_error(p, cov, chi2_red=None):
    """sigma(theta0) из 2-й гармоники — АНАЛИТИЧЕСКИ (дельта-метод), без бутстрапа.

    Если задан chi2_red, ковариация масштабируется на него (поправка Бирджа):
    это честно, когда абсолютный уровень шума не известен, а известна лишь его
    ФОРМА (у нас sigma_i ~ U_i с неизвестным множителем).
    """
    a2, b2 = float(p[1]), float(p[2])
    A = np.hypot(a2, b2)
    if A <= 0:
        return float("nan")
    v_a, v_b, c_ab = cov[1, 1], cov[2, 2], cov[1, 2]
    var = (b2 ** 2 * v_a + a2 ** 2 * v_b - 2 * a2 * b2 * c_ab) / (4 * A ** 4)
    if chi2_red is not None and np.isfinite(chi2_red):
        var *= chi2_red
    return float(np.sqrt(max(var, 0.0)) * R2D)


def channel_params(p, cov=None):
    """Точный случай: 5 гармонических коэффициентов -> физика каналов.

    В системе с офсетом: U = c0 + c2*cos2(t-t0) + c4*cos4(t-t0), c2,c4 >= 0.
    Тождества (вывод — в NOTE_malus_linearization.md):
        |t_perp|^2             = c0 + c2 + c4
        |t_par|^2              = c0 - c2 + c4
        Re(t_perp*conj(t_par)) = c0 - 3*c4
    Отсюда eta и КОСИНУС относительной фазы — модельно-независимо, из линейного МНК.
    ВАЖНО: фаза сидит ТОЛЬКО в 4-й гармонике. Учебная 3-параметрическая
    линеаризация выбрасывает ровно ту наблюдаемую, что несёт tau_leak.
    """
    c0, a2, b2, a4, b4 = (float(v) for v in p[:5])
    A2, A4 = float(np.hypot(a2, b2)), float(np.hypot(a4, b4))
    t0_h2 = 0.5 * np.arctan2(b2, a2) * R2D                 # офсет из 2-й гармоники
    t0_h4 = 0.25 * np.arctan2(b4, a4) * R2D                # он же из 4-й (mod 45)
    c2 = A2
    ph = 4 * t0_h2 * D2R
    c4 = a4 * np.cos(ph) + b4 * np.sin(ph)                 # знак c4 в системе офсета

    I_perp = c0 + c2 + c4
    I_par = c0 - c2 + c4
    cross = c0 - 3 * c4
    denom = np.sqrt(I_perp * I_par) if I_perp > 0 and I_par > 0 else np.nan
    cos_dphi = cross / denom if denom and np.isfinite(denom) else np.nan

    # согласование двух независимых оценок офсета (mod 45 град)
    d = (t0_h4 - t0_h2 + 22.5) % 45.0 - 22.5

    return {
        "theta0_deg_from_h2": t0_h2,
        "theta0_deg_from_h4_mod45": t0_h4,
        "theta0_h4_minus_h2_deg": float(d),
        "A2": A2, "A4": A4, "c4_signed": float(c4),
        "harmonic4_over_harmonic2": A4 / A2 if A2 else np.nan,
        "I_perp": float(I_perp), "I_par": float(I_par),
        "eta_amplitude": float(np.sqrt(I_par / I_perp)) if I_perp > 0 and I_par > 0 else np.nan,
        "extinction_db": float(10 * np.log10(I_perp / I_par)) if I_perp > 0 and I_par > 0 else np.inf,
        "Re_cross": float(cross),
        "cos_dphi": float(cos_dphi),
        "dphi_deg": float(np.degrees(np.arccos(np.clip(cos_dphi, -1, 1)))) if np.isfinite(cos_dphi) else np.nan,
        "cauchy_schwarz_ok": bool(np.isfinite(cos_dphi) and abs(cos_dphi) <= 1.0 + 1e-9),
    }


# ------------------------------------------------------ диагностика плана
def design_diagnostics(theta_deg, order=4, sigma_rel=0.01, U_scale=1.0, A2_rel=1.0):
    """Что угловой НАБОР (а не данные) делает с точностью офсета.

    Возвращает число обусловленности X^T X и предсказанную sigma(theta0)
    при единичном шуме — чистая геометрия плана эксперимента.
    """
    X = harmonic_design(theta_deg, order)
    XtX = X.T @ X
    cov = np.linalg.inv(XtX) * (sigma_rel * U_scale) ** 2
    a1, b1 = A2_rel * U_scale, 0.0                       # без потери общности t0=0
    A = np.hypot(a1, b1)
    var_t0 = (b1 ** 2 * cov[1, 1] + a1 ** 2 * cov[2, 2]
              - 2 * a1 * b1 * cov[1, 2]) / (4 * A ** 4)
    # неортогональность гармоник 2 и 4 на данном наборе углов
    n = len(theta_deg)
    G = XtX / n
    off24 = float(np.max(np.abs(G[1:3, 3:5]))) if order >= 4 else 0.0
    return {
        "n_angles": int(n),
        "theta_min": float(np.min(theta_deg)), "theta_max": float(np.max(theta_deg)),
        "cond_XtX": float(np.linalg.cond(XtX)),
        "sigma_theta0_deg": float(np.sqrt(max(var_t0, 0.0)) * R2D),
        "max_overlap_h2_h4": off24,
    }


# ------------------------------------------------------------- самопроверка
def _synth(theta_deg, t_perp, t_par, theta0_deg):
    t = (np.asarray(theta_deg, float) - theta0_deg) * D2R
    E = t_perp * np.cos(t) ** 2 + t_par * np.sin(t) ** 2
    return np.abs(E) ** 2


def selftest(verbose=True):
    checks = []

    def chk(name, ok, detail):
        checks.append({"check": name, "ok": bool(ok), "detail": detail})
        if verbose:
            print(f"  [{'OK ' if ok else 'FAIL'}] {name}: {detail}")

    rng = np.random.default_rng(20260729)
    th_full = np.arange(-90.0, 91.0, 10.0)               # 19 углов, как 356att
    th_half = np.arange(0.0, 91.0, 10.0)                 # 10 углов, как Фазы 1-3

    # 1. Точность: order=4 воспроизводит когерентную модель до машинного нуля
    t_perp, t_par, t0 = 0.93 + 0.11j, 0.031 * np.exp(1j * 0.7), 1.35
    U = _synth(th_full, t_perp, t_par, t0)
    p, cov, info = fit_harmonics(th_full, U, order=4)
    chk("order=4 точно описывает |t_perp cos^2 + t_par sin^2|^2",
        info["resid_rms"] < 1e-14, f"rms невязки = {info['resid_rms']:.2e}")

    # 2. order=2 НЕ точен для когерентной модели (есть 4-я гармоника)
    p2, cov2, info2 = fit_harmonics(th_full, U, order=2)
    chk("order=2 оставляет структурную невязку (4-я гармоника реальна)",
        info2["resid_rms"] > 1e3 * max(info["resid_rms"], 1e-16),
        f"rms(order=2) = {info2['resid_rms']:.3e} против {info['resid_rms']:.2e}")

    # 3. Замкнутые формулы возвращают физику канала
    ch = channel_params(p)
    dphi_true = np.angle(t_perp) - np.angle(t_par)
    eta_true = abs(t_par) / abs(t_perp)
    chk("восстановление theta0 (2-я гармоника)", abs(ch["theta0_deg_from_h2"] - t0) < 1e-9,
        f"{ch['theta0_deg_from_h2']:.12f} против {t0}")
    chk("восстановление eta", abs(ch["eta_amplitude"] - eta_true) < 1e-12,
        f"{ch['eta_amplitude']:.10f} против {eta_true:.10f}")
    chk("восстановление относительной фазы", abs(ch["cos_dphi"] - np.cos(dphi_true)) < 1e-10,
        f"cos(dphi) = {ch['cos_dphi']:.10f} против {np.cos(dphi_true):.10f}")
    chk("два независимых офсета (h2 и h4) согласованы",
        abs(ch["theta0_h4_minus_h2_deg"]) < 1e-8,
        f"h4-h2 = {ch['theta0_h4_minus_h2_deg']:.2e} град")

    # 4. Некогерентный (интенсивностный) случай: order=2 ТОЧЕН
    Ui = (abs(t_perp) ** 2 * np.cos((th_full - t0) * D2R) ** 2
          + abs(t_par) ** 2 * np.sin((th_full - t0) * D2R) ** 2)
    pi_, ci_, ii_ = fit_harmonics(th_full, Ui, order=2)
    mp = malus_params(pi_, ci_)
    chk("интенсивностная модель: order=2 точен", ii_["resid_rms"] < 1e-14,
        f"rms = {ii_['resid_rms']:.2e}")
    chk("интенсивностная модель: theta0 и утечка замкнуто",
        abs(mp["theta0_deg"] - t0) < 1e-9 and abs(mp["eta_amplitude"] - eta_true) < 1e-10,
        f"theta0={mp['theta0_deg']:.10f}, eta={mp['eta_amplitude']:.10f}, "
        f"экстинкция={mp['extinction_db']:.3f} дБ")

    # 5. Суммирование по частотам НЕ ломает гармоническую структуру
    nu = np.linspace(0.2, 1.5, 40)
    Usum = np.zeros_like(th_full)
    for f in nu:                                          # eta(nu) и фаза(nu) — произвольные
        tp = (0.95 - 0.1 * f) * np.exp(1j * 0.3 * f)
        tl = 0.02 * (1 + 0.5 * f) * np.exp(1j * (0.4 + 2 * np.pi * f * 0.19))
        Usum += _synth(th_full, tp, tl, t0)
    _, _, isum = fit_harmonics(th_full, Usum, order=4)
    chk("сумма по частотам (Парсеваль) остаётся многочленом 4-го порядка",
        isum["resid_rms"] / Usum.max() < 1e-14,
        f"относит. rms = {isum['resid_rms'] / Usum.max():.2e}")

    # 6. Логарифмическая шкала ЛОМАЕТ линейность
    _, _, ilog = fit_harmonics(th_full, np.log(U), order=4)
    chk("в лог-шкале (дБ) линейность разрушается — фитовать надо U, не ln U",
        ilog["resid_rms"] > 1e-3,
        f"rms(ln U) = {ilog['resid_rms']:.4f} нат = {ilog['resid_rms'] * 10 / np.log(10):.3f} дБ")

    # 7. Обусловленность: полный диапазон против половинного
    d_full = design_diagnostics(th_full, order=4)
    d_half = design_diagnostics(th_half, order=4)
    gain = d_half["sigma_theta0_deg"] / d_full["sigma_theta0_deg"]
    chk("полный диапазон -90..90 обусловлен лучше половинного 0..90",
        d_full["cond_XtX"] < d_half["cond_XtX"] and gain > 1.0,
        f"cond {d_full['cond_XtX']:.2f} против {d_half['cond_XtX']:.2f}; "
        f"sigma(theta0) лучше в {gain:.2f} раза")

    # 8. Смещение theta0 при подгонке order=2 к когерентным данным
    #    (4-я гармоника не ортогональна базису на неполном периоде)
    bias_full = malus_params(fit_harmonics(th_full, U, order=2)[0])["theta0_deg"] - t0
    U_half = _synth(th_half, t_perp, t_par, t0)
    bias_half = malus_params(fit_harmonics(th_half, U_half, order=2)[0])["theta0_deg"] - t0
    chk("пренебрежение 4-й гармоникой смещает theta0 тем сильнее, чем уже диапазон",
        abs(bias_half) > abs(bias_full),
        f"смещение: полный {bias_full:+.4f} град, половинный {bias_half:+.4f} град")

    n_ok = sum(c["ok"] for c in checks)
    if verbose:
        print(f"\nИТОГ самопроверки: {n_ok}/{len(checks)}")
    return {
        "checks": checks, "n_ok": n_ok, "n_total": len(checks),
        "design_full_range": d_full, "design_half_range": d_half,
        "theta0_bias_order2_deg": {"full_range": float(bias_full),
                                   "half_range": float(bias_half)},
    }


# ------------------------------------------------------- применение к данным
def _weights_for(U, mode):
    """Веса МНК. Выбор шкалы шума — НЕ косметика: он решает, слышна ли тень.

    'uniform'  : sigma = const. Вклад точки ~ U^2 -> яркая зона давит тень
                 полностью, eta и экстинкция фактически не определены.
    'relative' : sigma_i = U_i (постоянная ОТНОСИТЕЛЬНАЯ ошибка) -> все углы
                 весят одинаково в лог-смысле; тень начинает влиять.
                 Задача остаётся ЛИНЕЙНОЙ (веса не зависят от параметров).
    """
    if mode == "uniform":
        return None
    if mode == "relative":
        return np.maximum(np.asarray(U, float), 1e-300)
    raise ValueError(mode)


def analyse_dataset(dataset, angles_limit=None, verbose=True):
    """Гармонический разбор РЕАЛЬНЫХ данных: по каждой частоте отдельно + сумма.

    exp[theta, nu] = FFT(sig)/FFT(bg) — комплексное пропускание. Нормировка на
    опорный сигнал — общий по theta множитель, поэтому она сокращается в eta и
    в cos(dphi): обе величины ИНВАРИАНТНЫ к ней (и, по теореме выравнивания T8,
    к скалярному множителю t_perp2(nu) второго WGP).
    """
    import sys
    sys.path.insert(0, str(REPO / "research" / "two_wgp"))
    import exp_builder

    ang, freqs, exp, valid, meta = exp_builder.load(
        dataset, angles_limit=angles_limit, mask_angles=(0.0, 90.0))
    theta = np.asarray(ang, float)
    idx = np.where(valid[0])[0]
    P = np.abs(exp[:, idx]) ** 2                       # (n_ang, n_freq), мощность

    per_nu = []
    for j in range(P.shape[1]):
        col = P[:, j]
        p, _, info = fit_harmonics(theta, col, sigma=_weights_for(col, "relative"),
                                   order=4)
        ch = channel_params(p)
        ch["freq_thz"] = float(freqs[idx[j]])
        ch["resid_rms_rel"] = float(np.sqrt(info["chi2_red"]))   # среднее ОТНОСИТ. отклонение
        per_nu.append(ch)

    # энергия Парсеваля (сумма по частотам) — тоже точный многочлен 4-го порядка
    U = P.sum(axis=1)
    p_sum, cov_sum, info_sum = fit_harmonics(
        theta, U, sigma=_weights_for(U, "relative"), order=4)
    ch_sum = channel_params(p_sum)
    ch_sum["theta0_err_deg"] = theta0_error(p_sum, cov_sum, info_sum["chi2_red"])
    p_sum2, cov_sum2, info_sum2 = fit_harmonics(
        theta, U, sigma=_weights_for(U, "relative"), order=2)
    mp_sum2 = malus_params(p_sum2, cov_sum2)

    # --- две служебные сверки: шкала весов и урезанный диапазон углов
    p_u4, cov_u4, _ = fit_harmonics(theta, U, order=4)
    ch_u4 = channel_params(p_u4)
    half = (theta >= 0) & (theta <= 90)
    p_h4, cov_h4, _ = fit_harmonics(theta[half], U[half],
                                    sigma=_weights_for(U[half], "relative"), order=4)
    ch_h4 = channel_params(p_h4)
    p_h2, cov_h2, _ = fit_harmonics(theta[half], U[half],
                                    sigma=_weights_for(U[half], "relative"), order=2)
    mp_h2 = malus_params(p_h2, cov_h2)
    weight_sensitivity = {
        "uniform_theta0_deg": ch_u4["theta0_deg_from_h2"],
        "relative_theta0_deg": ch_sum["theta0_deg_from_h2"],
        "uniform_eta": ch_u4["eta_amplitude"],
        "relative_eta": ch_sum["eta_amplitude"],
        "uniform_I_par_negative": bool(ch_u4["I_par"] < 0),
        "comment": "невзвешенный МНК по мощности слеп к тени; eta оттуда недостоверна",
    }
    range_sensitivity = {
        "n_angles_half": int(half.sum()),
        "theta0_full_range_deg": ch_sum["theta0_deg_from_h2"],
        "theta0_half_range_deg": ch_h4["theta0_deg_from_h2"],
        "delta_deg": float(ch_h4["theta0_deg_from_h2"] - ch_sum["theta0_deg_from_h2"]),
        "theta0_half_range_order2_deg": mp_h2["theta0_deg"],
        "order2_bias_on_half_range_deg":
            float(mp_h2["theta0_deg"] - ch_h4["theta0_deg_from_h2"]),
        "cond_full": design_diagnostics(theta, 4)["cond_XtX"],
        "cond_half": design_diagnostics(theta[half], 4)["cond_XtX"],
        "sigma_theta0_ratio_half_over_full":
            design_diagnostics(theta[half], 4)["sigma_theta0_deg"]
            / design_diagnostics(theta, 4)["sigma_theta0_deg"],
    }

    nu = np.array([c["freq_thz"] for c in per_nu])
    dphi = np.array([c["dphi_deg"] for c in per_nu])
    ok = np.isfinite(dphi)
    tau_ps = slope = np.nan
    if ok.sum() > 5:
        # dphi(nu) = delta0 + 2*pi*nu*tau_leak  ->  наклон даёт tau_leak (пс)
        sl, b0 = np.polyfit(nu[ok], np.deg2rad(dphi[ok]), 1)
        slope, tau_ps = float(sl), float(sl / (2 * np.pi))   # nu в ТГц -> tau в пс

    out = {
        "dataset": dataset,
        "n_angles": int(len(theta)), "angles_deg": [float(a) for a in theta],
        "n_freqs": int(len(idx)),
        "design": design_diagnostics(theta, order=4),
        "parseval_order4": ch_sum,
        "parseval_order4_resid_rms_rel": float(np.sqrt(info_sum["chi2_red"])),
        "parseval_order2_malus": mp_sum2,
        "parseval_order2_resid_rms_rel": float(np.sqrt(info_sum2["chi2_red"])),
        "weight_sensitivity": weight_sensitivity,
        "range_sensitivity": range_sensitivity,
        "theta0_err_deg_analytic": mp_sum2.get("theta0_err_deg"),
        "theta0_bias_order2_minus_order4_deg":
            float(mp_sum2["theta0_deg"] - ch_sum["theta0_deg_from_h2"]),
        "per_freq_summary": {
            "eta_median": float(np.nanmedian([c["eta_amplitude"] for c in per_nu])),
            "eta_p16_p84": [float(np.nanpercentile([c["eta_amplitude"] for c in per_nu], 16)),
                            float(np.nanpercentile([c["eta_amplitude"] for c in per_nu], 84))],
            "theta0_h2_median_deg": float(np.nanmedian(
                [c["theta0_deg_from_h2"] for c in per_nu])),
            "theta0_h2_mad_deg": float(np.nanmedian(np.abs(
                np.array([c["theta0_deg_from_h2"] for c in per_nu])
                - np.nanmedian([c["theta0_deg_from_h2"] for c in per_nu])))),
            "theta0_h4_minus_h2_median_deg": float(np.nanmedian(
                [c["theta0_h4_minus_h2_deg"] for c in per_nu])),
            "cauchy_schwarz_violations": int(sum(not c["cauchy_schwarz_ok"] for c in per_nu)),
            "dphi_slope_rad_per_thz": slope,
            "tau_leak_ps_model_free": tau_ps,
        },
        "per_freq": per_nu,
    }

    if verbose:
        s = out["per_freq_summary"]
        print(f"\n=== {dataset} ===  углов {len(theta)}, частот {len(idx)}")
        print(f"  план: cond(X^T X) = {out['design']['cond_XtX']:.2f}, "
              f"sigma(theta0) при 1% шума = {out['design']['sigma_theta0_deg']:.4f} град")
        print(f"  Парсеваль, order=4: theta0(h2) = {ch_sum['theta0_deg_from_h2']:+.4f} "
              f"+- {ch_sum['theta0_err_deg']:.4f} град, "
              f"theta0(h4) = {ch_sum['theta0_deg_from_h4_mod45']:+.4f}, "
              f"расхожд. = {ch_sum['theta0_h4_minus_h2_deg']:+.4f}")
        print(f"      A4/A2 = {ch_sum['harmonic4_over_harmonic2']:.4f} "
              f"(идеальный cos^4 даёт ровно 0.25)")
        print(f"      eta = {ch_sum['eta_amplitude']:.5f}, "
              f"экстинкция = {ch_sum['extinction_db']:.2f} дБ, "
              f"cos(dphi) = {ch_sum['cos_dphi']:+.4f}")
        print(f"      отн. rms невязки: order=4 {out['parseval_order4_resid_rms_rel']:.2e}, "
              f"order=2 {out['parseval_order2_resid_rms_rel']:.2e}")
        print(f"  учебная 3-парам. форма даёт theta0 = {mp_sum2['theta0_deg']:+.4f} "
              f"+- {mp_sum2.get('theta0_err_deg', float('nan')):.4f} град "
              f"-> СМЕЩЕНИЕ {out['theta0_bias_order2_minus_order4_deg']:+.4f} град")
        w, r = weight_sensitivity, range_sensitivity
        print(f"  шкала весов: theta0 равномерн. {w['uniform_theta0_deg']:+.4f} / "
              f"относит. {w['relative_theta0_deg']:+.4f} град; "
              f"eta {w['uniform_eta']:.5f} / {w['relative_eta']:.5f}"
              f"{'  [равномерн. дала I_par<0!]' if w['uniform_I_par_negative'] else ''}")
        print(f"  урезание до 0..90 ({r['n_angles_half']} углов): theta0 "
              f"{r['theta0_full_range_deg']:+.4f} -> {r['theta0_half_range_deg']:+.4f} "
              f"(сдвиг {r['delta_deg']:+.4f} град); cond {r['cond_full']:.2f} -> "
              f"{r['cond_half']:.2f}, sigma(theta0) хуже в "
              f"{r['sigma_theta0_ratio_half_over_full']:.2f} раза; "
              f"смещение от отбрасывания 4-й гармоники на этом диапазоне "
              f"{r['order2_bias_on_half_range_deg']:+.3f} град")
        print(f"  по частотам: eta = {s['eta_median']:.5f} "
              f"[{s['eta_p16_p84'][0]:.5f}, {s['eta_p16_p84'][1]:.5f}], "
              f"theta0(h2) = {s['theta0_h2_median_deg']:+.4f} +- {s['theta0_h2_mad_deg']:.4f} (MAD)")
        print(f"      h4-h2 медиана = {s['theta0_h4_minus_h2_median_deg']:+.4f} град, "
              f"нарушений Коши-Буняковского {s['cauchy_schwarz_violations']}/{len(per_nu)}")
        print(f"      tau_leak модельно-независимо = {s['tau_leak_ps_model_free']:.4f} пс "
              f"(M5 даёт ~0.19 пс)")
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--selftest", action="store_true", help="прогнать самопроверку")
    ap.add_argument("--data", nargs="*", help="разобрать реальные датасеты")
    ap.add_argument("--out", default=str(RESULTS / "a0_malus_linearization.json"))
    args = ap.parse_args()

    if args.data is not None:
        ds = args.data or ["att-11-16-356", "att-11-16-s1", "att-11-16-s2", "test_grid_40_20"]
        res = {"datasets": {}, "note": "гармонический разбор, полный угловой диапазон"}
        for d in ds:
            try:
                res["datasets"][d] = analyse_dataset(d)
            except Exception as e:                       # noqa: BLE001
                print(f"  [пропуск] {d}: {type(e).__name__}: {e}")
                res["datasets"][d] = {"error": f"{type(e).__name__}: {e}"}
        out = RESULTS / "a0_malus_harmonics_data.json"
        out.parent.mkdir(parents=True, exist_ok=True)
        out.write_text(json.dumps(res, ensure_ascii=False, indent=2), encoding="utf-8")
        print(f"\nартефакт -> {out}")
        raise SystemExit(0)

    if args.selftest:
        res = selftest()
        Path(args.out).parent.mkdir(parents=True, exist_ok=True)
        Path(args.out).write_text(json.dumps(res, ensure_ascii=False, indent=2),
                                  encoding="utf-8")
        print(f"артефакт -> {args.out}")
        raise SystemExit(0 if res["n_ok"] == res["n_total"] else 1)
    ap.print_help()


if __name__ == "__main__":
    main()
