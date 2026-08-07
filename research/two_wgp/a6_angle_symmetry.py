"""
A6 -- МОДЕЛЬНО-НЕЗАВИСИМЫЙ тест угловой симметрии U(theta) vs U(-theta).

ИДЕЯ. Все фиты Фаз 1-3 и A2 шли на углах 0..90 (config.ANGLES_LIMIT_2D), хотя данные
покрывают -90..+90 (356att, series1/2) и даже -100..+100 (test_grid_40_20). Между тем:
  * угловой офсет theta0 порождает АСИММЕТРИЮ U(+theta) != U(-theta);
  * физическая утечка eta и сниженный D_eff -- СИММЕТРИЧНЫ по theta.
Значит полный диапазон даёт рычаг на вырожденность angle_offset <-> phi2 <-> eps_cross,
которую A2 закрыть не смогла (корр. -0.95, два бассейна с разрывом 110 AIC), причём БЕЗ
всякой модели и без фита.

ЧТО СЧИТАЕТСЯ (U -- Парсеваль по комплексному пропусканию T = FFT(sig)/FFT(bg)):
    U(theta) = sum_{nu in valid} |T(theta,nu)|^2 ,   S(theta) = ln U(theta)
Профиль S симметричен относительно theta0: S(theta) = f(theta - theta0), f чётная.

ТРИ НЕЗАВИСИМЫХ ОЦЕНКИ theta0:
  E1 (по-парная, без интерполяции и без модели). Delta(theta) = S(theta) - S(-theta)
     ~= -2*theta0*f'(theta), где f оценивается симметризацией f_sym=(S(+)+S(-))/2.
     Даёт СВОЙ theta0 на каждой паре -> проверка согласованности: настоящий офсет даёт
     ПЛОСКИЙ профиль theta0(theta), любая другая причина асимметрии -- разброс.
  E2 (глобальная, сплайн). theta0 = argmin int [S(theta0+x) - S(theta0-x)]^2 w(x) dx.
  E3 (параметрическая сверка, HN7-форма). U = U0*cos^4(theta-theta0) + C; побочно даёт
     пол утечки eps = C/U0, сравнимый с eta.
  E4 (положение минимума экстинкции). Параболой по ln U вокруг минимума; ВНУТРИ диапазона
     он лежит только у test_grid_40_20 (+-100) и specac (0..100) -> прямая сверка с
     "минимум Specac ~84 град" (T16).

ДАЛЕЕ: остаточная асимметрия ПОСЛЕ снятия theta0 = верхняя граница на любой истинно
невзаимный / eps_cross-подобный член (он, в отличие от офсета, не сводится к сдвигу оси).

Ошибки -- бутстрап по частотным отсчётам (B выборок с возвратом), пересчёт всей цепочки.

Запуск (из корня репозитория):
    .venv\\Scripts\\python.exe research/two_wgp/a6_angle_symmetry.py
    ... --datasets att-11-16-356 att-11-16-s1 --nboot 800 --no-plot

Сессия A. Билдер -- research/two_wgp/exp_builder.py (бит-в-бит сверен с fit_lib на 0..90).
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
from scipy.interpolate import CubicSpline
from scipy.optimize import minimize_scalar, least_squares

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO))
sys.path.insert(0, str(REPO / "research" / "two_wgp"))

import exp_builder                                        # noqa: E402
import malus_linearization as ml                          # noqa: E402  (A0: линейный оценщик)
from unified_optimizer import config                      # noqa: E402
from unified_optimizer.data_manager import DataManager     # noqa: E402

OUT_DIR = REPO / "research" / "results" / "two_wgp"
NEP = 10.0 / np.log(10.0)      # ln -> dB (для энергии: 10*log10)

# Опубликованные angle_offset из A2 (research/results/two_wgp/a2_356att_ladder.json),
# 356att, фит на 10 углах 0..90 -- то, с чем сверяем модельно-независимый theta0.
M5_ANGLE_OFFSET_356ATT = {
    "M3ref_1wgp": (-1.2939, 0.0513),
    "M5_1wgp": (0.2190, 0.0650),
    "M5_eps": (0.5839, 0.0701),
    "W2_mis_psi_eps_seedEps": (-0.1994, 1957.2),
}


# ---------------------------------------------------------------- энергия U(theta)
def energy_profile(exp, valid, freq_idx=None):
    """U(theta) = sum_nu |T|^2 по валидным частотам (маска одинакова для всех углов)."""
    col = valid[0]
    if not np.all(valid == col[None, :]):
        raise AssertionError("маска воды зависит от угла -- U(theta) несравнимы")
    idx = np.where(col)[0] if freq_idx is None else freq_idx
    return np.sum(np.abs(exp[:, idx]) ** 2, axis=1)


# ---------------------------------------------------------------- E1: по-парная
def pairwise_offset(theta, S, sigma_S=None):
    """theta0 из каждой пары (+theta,-theta): theta0 = -Delta/(2 f'_sym).

    Возвращает список записей на каждую пару. Ключевое -- СОГЛАСОВАННОСТЬ значений.
    """
    th = np.asarray(theta, float)
    S = np.asarray(S, float)
    lut = {round(float(t), 6): s for t, s in zip(th, S)}
    pos = sorted(t for t in th if t > 0 and round(-t, 6) in lut)

    # симметризованный профиль f_sym(|theta|) (чётная часть; ошибка O(theta0^2))
    absth = [0.0] + pos if 0.0 in set(th) else pos
    f_sym = {}
    for t in absth:
        if t == 0.0:
            f_sym[t] = lut[0.0]
        else:
            f_sym[t] = 0.5 * (lut[round(t, 6)] + lut[round(-t, 6)])
    keys = sorted(f_sym)

    out = []
    for t in pos:
        i = keys.index(t)
        if i == 0 or i == len(keys) - 1:
            continue                                     # нет центральной разности
        tm, tp = keys[i - 1], keys[i + 1]
        fprime = (f_sym[tp] - f_sym[tm]) / (tp - tm)      # d f_sym / d theta
        if abs(fprime) < 1e-12:
            continue
        delta = lut[round(t, 6)] - lut[round(-t, 6)]
        t0 = -delta / (2.0 * fprime)
        rec = {"theta_deg": float(t),
               "delta_db": float(delta * NEP),
               "dS_dtheta_per_deg": float(fprime),
               "theta0_deg": float(t0)}
        if sigma_S is not None:
            sd = float(np.hypot(sigma_S[list(th).index(t)],
                                sigma_S[list(th).index(-t)]))
            rec["sigma_delta_db"] = sd * NEP
            rec["sigma_theta0_deg"] = sd / (2.0 * abs(fprime))
        out.append(rec)
    return out


def weighted_mean(vals, errs):
    """Взвешенное среднее + ошибка + chi2_red СОГЛАСОВАННОСТИ.

    chi2_red >> 1 означает, что разброс между парами не объясняется заявленными
    ошибками -> оценки НЕ описывают одну величину (для theta0: асимметрия имеет не
    офсетную природу). Тогда честная ошибка -- инфляция по Бирджу sqrt(chi2_red).
    """
    v = np.asarray(vals, float)
    e = np.asarray(errs, float)
    good = np.isfinite(v) & np.isfinite(e) & (e > 0)
    if good.sum() == 0:
        return float(np.mean(v)), float("nan"), float("nan")
    w = 1.0 / e[good] ** 2
    m = float(np.sum(w * v[good]) / np.sum(w))
    sm = float(np.sqrt(1.0 / np.sum(w)))
    chi2 = float(np.sum(w * (v[good] - m) ** 2) / max(good.sum() - 1, 1))
    return m, sm, chi2


def angle_level_noise(theta, S, c4, bright_max=40.0):
    """ВЕРХНЯЯ граница ошибки на уровне УГЛА (то, чего бутстрап по частотам не видит).

    Бутстрап по частотам ловит только частотный разброс внутри одного скана. Между
    сканами есть дрейф лазера и перепозиционирование ротатора -- на угол приходится
    ОДИН повтор (356att/series/test_grid), независимой оценки нет. Прокси: невязка
    S(theta) к гладкой форме cos^4+пол (E3). Она включает и неадекватность модели,
    поэтому это ВЕРХНЯЯ граница; бутстрап -- НИЖНЯЯ. Истина между ними.
    """
    th = np.asarray(theta, float)
    S = np.asarray(S, float)
    pred = np.log(c4["U0"] * np.cos(np.deg2rad(th - c4["theta0_deg"])) ** 4 + c4["floor"])
    r = (S - pred) * NEP
    bright = np.abs(th) <= bright_max
    shadow = ~bright
    return {
        "rms_db_bright": float(np.sqrt(np.mean(r[bright] ** 2))) if bright.any() else None,
        "rms_db_shadow": float(np.sqrt(np.mean(r[shadow] ** 2))) if shadow.any() else None,
        "bright_max_deg": bright_max,
        "note": "невязка к cos^4+пол; ВЕРХНЯЯ граница ошибки уровня угла",
    }


# ---------------------------------------------------------------- E2: сплайн
def spline_symmetry_offset(theta, S, x_max=None, bound=12.0, dx=0.5):
    """theta0 = argmin int_0^X [S(t0+x) - S(t0-x)]^2 dx (кубический сплайн по S)."""
    th = np.asarray(theta, float)
    order = np.argsort(th)
    sp = CubicSpline(th[order], np.asarray(S, float)[order])
    lo, hi = float(th.min()), float(th.max())
    x_cap = x_max if x_max is not None else min(hi, -lo)

    def J(t0):
        xm = min(x_cap, hi - t0, t0 - lo)
        if xm <= dx:
            return 1e9
        x = np.arange(dx, xm + 1e-9, dx)
        d = sp(t0 + x) - sp(t0 - x)
        return float(np.mean(d ** 2))

    r = minimize_scalar(J, bounds=(-bound, bound), method="bounded",
                        options={"xatol": 1e-4})
    return float(r.x), sp, float(r.fun)


def residual_asymmetry(sp, t0, lo, hi, x_max=None, dx=1.0):
    """Остаточная асимметрия после снятия theta0 (в дБ) -- граница на невзаимный член."""
    xm = min(x_max if x_max is not None else min(hi, -lo), hi - t0, t0 - lo)
    x = np.arange(dx, xm + 1e-9, dx)
    d = (sp(t0 + x) - sp(t0 - x)) * NEP
    return x, d


# ---------------------------------------------------------------- E3: cos^4 + пол
def cos4_fit(theta, U, sigma_U=None):
    """U = U0*cos^4(theta - theta0) + C. Возвращает (U0, theta0, C, eps=C/U0)."""
    th = np.asarray(theta, float)
    U = np.asarray(U, float)
    w = np.ones_like(U) if sigma_U is None else 1.0 / np.maximum(sigma_U, 1e-30)

    def res(p):
        U0, t0, C = p
        return (U0 * np.cos(np.deg2rad(th - t0)) ** 4 + C - U) * w

    p0 = [float(U.max()), 0.0, float(max(U.min(), 1e-12))]
    r = least_squares(res, p0, bounds=([0, -20, 0], [np.inf, 20, np.inf]))
    U0, t0, C = r.x
    return {"U0": float(U0), "theta0_deg": float(t0), "floor": float(C),
            "eps_floor": float(C / U0) if U0 > 0 else float("nan"),
            "cost": float(r.cost)}


# ---------------------------------------------------------------- E4: минимум
def extinction_minimum(theta, S, near, half_window=15.0):
    """Положение минимума ln U параболой по точкам в окне вокруг `near`.

    Возвращает None, если минимум лежит на КРАЮ покрытия (не разрешается).
    """
    th = np.asarray(theta, float)
    S = np.asarray(S, float)
    sel = np.abs(th - near) <= half_window
    if sel.sum() < 3:
        return {"edge_limited": True, "n_points": int(sel.sum()),
                "reason": "покрытие обрывается: <3 точек в окне вокруг минимума"}
    ts, ss = th[sel], S[sel]
    i_min = int(np.argmin(ss))
    if i_min == 0 or i_min == len(ts) - 1:
        return {"edge_limited": True, "argmin_deg": float(ts[i_min])}
    c = np.polyfit(ts, ss, 2)
    if c[0] <= 0:
        return {"edge_limited": True, "argmin_deg": float(ts[i_min])}
    return {"edge_limited": False, "theta_min_deg": float(-c[1] / (2 * c[0])),
            "argmin_deg": float(ts[i_min]), "n_points": int(sel.sum())}


# ---------------------------------------------------------------- анализ датасета
def _offset_by_zone(pairs, bright_max, chi2_global):
    """theta0 отдельно по яркой зоне и по тени — ТЕСТ НА ОФСЕТНУЮ ФОРМУ.

    Угловой офсет — ОДНО число: он обязан давать одинаковый theta0 из любой пары
    углов. Если оценки по яркой зоне и по тени расходятся значимо, значит
    асимметрия НЕ сводится к сдвигу нуля шкалы, и подгонять её одним параметром
    `angle_offset` некорректно. Ошибки инфлируются на sqrt(chi2) (поправка
    Бирджа): разброс между парами сам по себе — свидетельство недомоделирования.
    """
    if not pairs:
        return None
    br = [p for p in pairs if abs(p["theta_deg"]) <= bright_max]
    sh = [p for p in pairs if abs(p["theta_deg"]) > bright_max]

    def agg(ps):
        if len(ps) < 2:
            return None
        m, s, c = weighted_mean([p["theta0_deg"] for p in ps],
                                [p.get("sigma_theta0_deg", np.nan) for p in ps])
        return {"theta0_deg": m, "err_deg": s * max(np.sqrt(c), 1.0),
                "err_raw_deg": s, "chi2_red": c, "n_pairs": len(ps)}

    b, s_ = agg(br), agg(sh)
    out = {"bright": b, "shadow": s_, "bright_max_deg": bright_max,
           "chi2_red_global": chi2_global}
    if b and s_:
        d = s_["theta0_deg"] - b["theta0_deg"]
        sd = float(np.hypot(b["err_deg"], s_["err_deg"]))
        out["shadow_minus_bright_deg"] = float(d)
        out["shadow_minus_bright_err_deg"] = sd
        out["n_sigma"] = float(abs(d) / sd) if sd > 0 else None
        out["pure_offset_form"] = bool(sd > 0 and abs(d) / sd < 3.0)
        out["interpretation"] = (
            "theta0 из яркой зоны и из тени СОВПАДАЮТ -> асимметрия совместима "
            "с чистым угловым офсетом"
            if out["pure_offset_form"] else
            "theta0 из яркой зоны и из тени РАСХОДЯТСЯ -> асимметрия НЕ сводится "
            "к одному angle_offset; подгонка её одним параметром некорректна")
    return out


def _asym_verdict(xr, dres, boot_rms, noise, bright_max):
    """Остаточная асимметрия против ДВУХ границ шума — нижней и верхней.

    Урок A2: одна нижняя граница (бутстрап по частотам) завышает значимость,
    потому что не видит шума УРОВНЯ УГЛА. Вердикт выносится по ВЕРХНЕЙ границе:
    значимо только то, что её превышает. Яркая зона и тень разделены — там
    разные механизмы шума и разные амплитуды эффекта.
    """
    xr = np.asarray(xr, float)
    dres = np.asarray(dres, float)
    br = xr <= bright_max
    sh = ~br

    def blk(sel, up):
        if not np.any(sel):
            return None
        rms = float(np.sqrt(np.mean(dres[sel] ** 2)))
        out = {"rms_db": rms, "max_abs_db": float(np.max(np.abs(dres[sel]))),
               "n": int(sel.sum()), "noise_upper_db": up}
        out["significance_vs_upper"] = float(rms / up) if up else None
        out["significant"] = bool(up and rms > up)
        return out

    up_b = noise.get("rms_db_bright")
    up_s = noise.get("rms_db_shadow")
    rms_all = float(np.sqrt(np.mean(dres ** 2)))
    return {
        "x_deg": [float(v) for v in xr],
        "delta_db": [float(v) for v in dres],
        "max_abs_db": float(np.max(np.abs(dres))),
        "rms_db": rms_all,
        "noise_lower_boot_db": boot_rms,
        "significance_vs_lower_boot": (float(rms_all / boot_rms) if boot_rms else None),
        "bright": blk(br, up_b),
        "shadow": blk(sh, up_s),
        "verdict_rule": "значимо ТОЛЬКО если rms превышает ВЕРХНЮЮ границу шума "
                        "(невязка к cos^4+пол). Бутстрап по частотам — нижняя граница, "
                        "по ней значимость систематически завышена (урок A2).",
    }


def analyse(dataset, dm, nboot=600, x_max=70.0, bright_max=40.0,
            seed=20260728, verbose=True):
    ang_lim, freqs, exp, valid, meta = exp_builder.load(
        dataset, angles_limit=None, mask_angles=(0.0, 90.0), dm=dm)
    # контрольный прогон: маска, посчитанная по ВСЕМ углам (проверка чувствительности)
    _, _, _, valid_full, meta_full = exp_builder.load(
        dataset, angles_limit=None, mask_angles=None, dm=dm)

    theta = np.asarray(ang_lim, float)
    U = energy_profile(exp, valid)
    S = np.log(U)
    idx = np.where(valid[0])[0]

    # --- бутстрап по частотам
    rng = np.random.default_rng(seed)
    P = np.abs(exp[:, idx]) ** 2                          # (n_ang, n_freq)
    nf = P.shape[1]
    boot_S = np.empty((nboot, len(theta)))
    boot_t0_sp, boot_t0_c4, boot_eps = [], [], []
    for b in range(nboot):
        take = rng.integers(0, nf, nf)
        Ub = P[:, take].sum(axis=1) * 1.0
        boot_S[b] = np.log(Ub)
        t0b, _, _ = spline_symmetry_offset(theta, boot_S[b], x_max=x_max)
        boot_t0_sp.append(t0b)
        c4b = cos4_fit(theta, Ub)
        boot_t0_c4.append(c4b["theta0_deg"])
        boot_eps.append(c4b["eps_floor"])
    sigma_S = boot_S.std(axis=0, ddof=1)
    sigma_U = U * sigma_S                                  # dU = U * dlnU

    # --- E1
    pairs = pairwise_offset(theta, S, sigma_S)
    if pairs:
        m1, s1, chi1 = weighted_mean([p["theta0_deg"] for p in pairs],
                                     [p.get("sigma_theta0_deg", np.nan) for p in pairs])
    else:
        m1 = s1 = chi1 = float("nan")

    # РАЗДЕЛЕНИЕ по зонам: если асимметрия — чистый угловой офсет, оба значения
    # обязаны совпасть. Расхождение = асимметрия НЕ офсетной формы.
    e1_zones = _offset_by_zone(pairs, bright_max, chi1)

    # --- E2
    t0_sp, sp, Jmin = spline_symmetry_offset(theta, S, x_max=x_max)
    s_sp = float(np.std(boot_t0_sp, ddof=1)) if nboot > 2 else float("nan")
    xr, dres = residual_asymmetry(sp, t0_sp, theta.min(), theta.max(), x_max=x_max)

    # шум остаточной асимметрии: та же величина, посчитанная по бутстрап-репликам
    boot_res = []
    for b in range(min(nboot, 200)):
        _, spb, _ = spline_symmetry_offset(theta, boot_S[b], x_max=x_max)
        _, db = residual_asymmetry(spb, boot_t0_sp[b], theta.min(), theta.max(), x_max=x_max)
        boot_res.append(db)
    boot_res = np.array(boot_res)
    res_noise_rms = float(np.sqrt(np.mean(boot_res.std(axis=0, ddof=1) ** 2)))

    # --- E3 / E4
    c4_full = cos4_fit(theta, U, sigma_U)
    half = theta >= 0
    c4_half = cos4_fit(theta[half], U[half], sigma_U[half])
    min_pos = extinction_minimum(theta, S, near=+90.0)
    min_neg = extinction_minimum(theta, S, near=-90.0)

    # --- E5: гармонический (ЛИНЕЙНЫЙ) оценщик из A0 — замкнутая формула,
    #     единственный минимум, аналитическая ошибка. Опорная точка для E1-E4.
    p5, cov5, info5 = ml.fit_harmonics(theta, U, sigma=np.maximum(U, 1e-300), order=4)
    ch5 = ml.channel_params(p5)
    ch5["theta0_err_deg"] = ml.theta0_error(p5, cov5, info5["chi2_red"])
    ch5["resid_rms_rel"] = float(np.sqrt(info5["chi2_red"]))
    ch5["design"] = ml.design_diagnostics(theta, order=4)

    # --- шум уровня УГЛА: ВЕРХНЯЯ граница (бутстрап по частотам даёт НИЖНЮЮ)
    noise = angle_level_noise(theta, S, c4_full, bright_max)

    out = {
        "dataset": dataset,
        "coverage": {"angles": [float(a) for a in theta],
                     "n_angles": int(len(theta)),
                     "n_angles_in_phase123": int(np.sum((theta >= 0) & (theta <= 90))),
                     "n_freqs_valid": int(len(idx)),
                     "f_range_thz": [float(freqs[0]), float(freqs[-1])],
                     "mask_from": "angles 0..90 (как в Фазах 1-3)",
                     "n_valid_freqs_mask_all_angles": int(np.sum(valid_full[0])),
                     "mask_shift_freqs": int(np.sum(valid[0] != valid_full[0]))},
        "U_theta": {f"{a:+.0f}": float(u) for a, u in zip(theta, U)},
        "U_db_rel_max": {f"{a:+.0f}": float(10 * np.log10(u / U.max()))
                         for a, u in zip(theta, U)},
        "sigma_S_boot": {f"{a:+.0f}": float(s) for a, s in zip(theta, sigma_S)},
        "E1_pairwise": {"pairs": pairs, "theta0_wmean_deg": m1,
                        "theta0_wmean_err_deg": s1, "chi2_red_consistency": chi1,
                        "by_zone": e1_zones},
        "E2_spline": {"theta0_deg": t0_sp, "theta0_boot_err_deg": s_sp,
                      "J_min": Jmin, "x_max_deg": x_max},
        "E3_cos4_full": c4_full,
        "E3_cos4_halfrange_0_90": c4_half,
        "E3_theta0_boot_err_deg": (float(np.std(boot_t0_c4, ddof=1)) if nboot > 2 else None),
        "E3_eps_boot_err": (float(np.std(boot_eps, ddof=1)) if nboot > 2 else None),
        "E4_extinction_min_plus": min_pos,
        "E4_extinction_min_minus": min_neg,
        "E5_harmonic_linear": ch5,
        "angle_level_noise": noise,
        "residual_asymmetry_after_offset": _asym_verdict(
            xr, dres, res_noise_rms, noise, bright_max),
        "raw_asymmetry_before_offset_db": {
            f"{p['theta_deg']:.0f}": p["delta_db"] for p in pairs},
        "nboot": nboot,
    }

    if verbose:
        print(f"\n=== {dataset} ===  углов {len(theta)} "
              f"(в Фазах 1-3 использовалось {out['coverage']['n_angles_in_phase123']}), "
              f"частот {len(idx)}")
        print(f"  U(theta), дБ отн. макс.: " +
              " ".join(f"{a:+.0f}:{10*np.log10(u/U.max()):+.2f}" for a, u in zip(theta, U)))
        print(f"  СЫРАЯ асимметрия 10log10(U(+t)/U(-t)), дБ: " +
              " ".join(f"{p['theta_deg']:.0f}:{p['delta_db']:+.3f}" for p in pairs))
        print(f"  E1 по-парно theta0: " +
              " ".join(f"{p['theta_deg']:.0f}:{p['theta0_deg']:+.3f}" for p in pairs))
        print(f"  E1 взвеш.среднее theta0 = {m1:+.4f} +- {s1:.4f} град "
              f"(chi2_red согласованности {chi1:.2f})")
        if e1_zones and e1_zones.get("bright") and e1_zones.get("shadow"):
            z = e1_zones
            print(f"  E1 ПО ЗОНАМ: яркая {z['bright']['theta0_deg']:+.4f}"
                  f"+-{z['bright']['err_deg']:.4f} против тень "
                  f"{z['shadow']['theta0_deg']:+.4f}+-{z['shadow']['err_deg']:.4f} "
                  f"-> разность {z['shadow_minus_bright_deg']:+.4f} "
                  f"({z['n_sigma']:.1f} sigma) "
                  f"{'ОФСЕТНАЯ ФОРМА' if z['pure_offset_form'] else 'НЕ ОФСЕТНАЯ ФОРМА'}")
        print(f"  E2 сплайн         theta0 = {t0_sp:+.4f} +- {s_sp:.4f} град")
        print(f"  E3 cos^4 полный   theta0 = {c4_full['theta0_deg']:+.4f} град, "
              f"пол eps = {c4_full['eps_floor']:.5f}")
        print(f"  E3 cos^4 ТОЛЬКО 0..90  theta0 = {c4_half['theta0_deg']:+.4f} град, "
              f"пол eps = {c4_half['eps_floor']:.5f}")
        print(f"  E4 минимум экстинкции: +{90}: {min_pos}, -{90}: {min_neg}")
        print(f"  E5 гармонич.(лин.) theta0 = {ch5['theta0_deg_from_h2']:+.4f} "
              f"+- {ch5['theta0_err_deg']:.4f} град  [замкнутая формула, ОПОРНАЯ ТОЧКА]; "
              f"сверка h4-h2 = {ch5['theta0_h4_minus_h2_deg']:+.4f}, "
              f"A4/A2 = {ch5['harmonic4_over_harmonic2']:.4f}, eta = {ch5['eta_amplitude']:.5f}")
        print(f"     план: cond = {ch5['design']['cond_XtX']:.2f}, "
              f"sigma(theta0)|1% = {ch5['design']['sigma_theta0_deg']:.4f} град")
        print(f"  ШУМ УРОВНЯ УГЛА (верхняя граница): яркая зона "
              f"{noise['rms_db_bright']:.4f} дБ, тень {noise['rms_db_shadow']:.4f} дБ; "
              f"бутстрап по частотам (нижняя) {res_noise_rms:.4f} дБ")
        r = out["residual_asymmetry_after_offset"]
        print(f"  ОСТАТОЧНАЯ асимметрия после снятия theta0: "
              f"rms={r['rms_db']:.4f} дБ, max={r['max_abs_db']:.4f} дБ")
        for nm in ("bright", "shadow"):
            b = r[nm]
            if b:
                print(f"     {nm:>6}: rms={b['rms_db']:.4f} дБ против верхней границы "
                      f"{b['noise_upper_db']:.4f} дБ -> {b['significance_vs_upper']:.2f}x "
                      f"{'ЗНАЧИМО' if b['significant'] else 'не значимо'}")
        print(f"     (по нижней границе бутстрапа было бы "
              f"{r['significance_vs_lower_boot']:.2f}x — ЗАВЫШЕНО, урок A2)")
    return out


# ---------------------------------------------------------------- график
def make_plot(results, path):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    n = len(results)
    fig, axes = plt.subplots(2, n, figsize=(4.2 * n, 7.0), squeeze=False)
    for j, r in enumerate(results):
        th = np.array(r["coverage"]["angles"])
        udb = np.array([r["U_db_rel_max"][f"{a:+.0f}"] for a in th])
        ax = axes[0][j]
        ax.plot(th, udb, "o-", lw=1.2, ms=4, color="#2b6cb0")
        ax.axvline(0, color="#999", lw=0.7)
        ax.axvspan(0, 90, color="#f0ad4e", alpha=0.13,
                   label="использовалось в Фазах 1-3")
        t0 = r["E2_spline"]["theta0_deg"]
        ax.axvline(t0, color="#c53030", ls="--", lw=1.1,
                   label=f"theta0(E2) = {t0:+.2f}°")
        ax.set_title(f"{r['dataset']}: U(θ), {len(th)} углов")
        ax.set_xlabel("θ, град"); ax.set_ylabel("10 lg U/U_max, дБ")
        ax.legend(fontsize=7); ax.grid(alpha=0.25)

        ax = axes[1][j]
        pairs = r["E1_pairwise"]["pairs"]
        xs = [p["theta_deg"] for p in pairs]
        ys = [p["delta_db"] for p in pairs]
        es = [p.get("sigma_delta_db", 0.0) for p in pairs]
        ax.errorbar(xs, ys, yerr=es, fmt="s-", ms=4, lw=1.1, color="#276749",
                    label="сырая асимметрия")
        rr = r["residual_asymmetry_after_offset"]
        ax.plot(rr["x_deg"], rr["delta_db"], "-", lw=1.1, color="#c53030",
                label="после снятия θ0")
        ax.axhline(0, color="#999", lw=0.7)
        ax.set_xlabel("|θ|, град"); ax.set_ylabel("10 lg U(+θ)/U(−θ), дБ")
        ax.set_title("асимметрия: до и после снятия офсета")
        ax.legend(fontsize=7); ax.grid(alpha=0.25)
    fig.tight_layout()
    fig.savefig(path, dpi=140)
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser(description="A6: тест симметрии U(theta) vs U(-theta)")
    ap.add_argument("--datasets", nargs="*",
                    default=["att-11-16-356", "att-11-16-s1", "att-11-16-s2", "test_grid_40_20"])
    ap.add_argument("--nboot", type=int, default=600)
    ap.add_argument("--x-max", type=float, default=70.0,
                    help="макс. плечо |x| в градусах для сплайновой оценки (глубокая тень шумна)")
    ap.add_argument("--no-plot", action="store_true")
    ap.add_argument("--out", default=str(OUT_DIR / "a6_angle_symmetry.json"))
    args = ap.parse_args()

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    dm = DataManager(config.DATA_DIR)

    results = []
    for ds in args.datasets:
        try:
            results.append(analyse(ds, dm, nboot=args.nboot, x_max=args.x_max))
        except Exception as exc:
            print(f"\n=== {ds} === ПРОПУЩЕН: {type(exc).__name__}: {exc}")

    payload = {
        "task": "A6_angle_symmetry",
        "session": "A",
        "date": "2026-07-28",
        "method": "U(theta)=sum_nu |T|^2 (Парсеваль по комплексному пропусканию); "
                  "маска воды взята по углам 0..90 (как в Фазах 1-3), углы -- полные; "
                  "ошибки -- бутстрап по частотам",
        "builder": "research/two_wgp/exp_builder.py (бит-в-бит == fit_lib.build_experiment на 0..90)",
        "reference_M5_angle_offset_356att": M5_ANGLE_OFFSET_356ATT,
        "results": results,
    }
    with open(args.out, "w", encoding="utf-8") as fh:
        json.dump(payload, fh, ensure_ascii=False, indent=2)
    print(f"\nJSON -> {args.out}")

    if not args.no_plot and results:
        png = str(Path(args.out).with_suffix(".png"))
        try:
            make_plot(results, png)
            print(f"PNG  -> {png}")
        except Exception as exc:
            print(f"график не построен: {type(exc).__name__}: {exc}")


if __name__ == "__main__":
    main()
