#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""A22 — восстановление спектров каналов: ПОЛЕ против МОЩНОСТИ, на данных.

Этап 1 двухэтапной схемы: модельно-независимо восстановить T_perp(nu) и
T_par(nu), не подгоняя ни одного параметра модели Бланко. Обе ветки считаются
на ОДНИХ данных, чтобы разница между областями была измерена, а не обсуждена.

Что выдаётся (всё модельно-независимо)
--------------------------------------
theta0 + sigma           замкнутой формулой и профилем chi^2
lambda2/lambda1          тест единого theta0 = модельно-независимый тест HN13
T_perp(nu), T_par(nu)    комплексные, с ковариацией
eta(nu), eta0, eta_exp   прямое измерение HN2 без модели
dphi(nu), tau_leak       развёрнутая ЗНАКОВАЯ фаза; не требует пина D_eff (HN8)
RSS_perp(nu)             невязка, недостижимая никакой моделью cos^2*T+sin^2*T
4-я гармоника в поле     тест выхода за схему «один WGP x скаляр»
дрейф уровня угла        комплексный множитель на трассу (HN6)
шум по повторам          прямая оценка sigma уровня угла (там, где есть повторы)

Гардрейлы: `data_pool/**`, `unified_optimizer/**` — только чтение; ядро зоны B
(`model_core`, `fit_lib`) не правится; `malus_linearization` импортируется, а
не копируется (урок C-1).

Запуск из КОРНЯ репозитория:
    .venv\\Scripts\\python.exe research/two_wgp/a22_field_vs_power.py --dataset test_grid_33_11
    .venv\\Scripts\\python.exe research/two_wgp/a22_field_vs_power.py --dataset specac-2
    .venv\\Scripts\\python.exe research/two_wgp/a22_field_vs_power.py --all
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
for p in (str(REPO), str(HERE)):
    if p not in sys.path:
        sys.path.insert(0, p)

import field_linearization as fl                                   # noqa: E402
import scan_loader                                                 # noqa: E402
import malus_linearization as ml                                   # noqa: E402

RESULTS = REPO / "research" / "results" / "two_wgp"

#: датасеты, читаемые каталогом (журнальный формат), а не через DataManager
DIR_DATASETS = {
    "specac-2": REPO / "data_pool" / "specac" / "exp_data_2",
}

#: известные числа проекта для сверки (модельно-независимые)
REFERENCE = {
    "test_grid_33_11": {"theta0_deg": -1.2770, "sigma_theta0_deg": 0.0354,
                        "A4_over_A2": 0.25183, "eta_parseval": 0.05281,
                        "source": "A13, research/results/two_wgp/a13_test_grid_33_11.json"},
    "att-11-16-356": {"tau_leak_ps": 0.195, "cs_violations": "12/151",
                      "source": "NOTE_malus_linearization.md §6.3"},
    "att-11-16-s1": {"tau_leak_ps": 0.269, "cs_violations": "43/94",
                     "source": "NOTE_malus_linearization.md §6.3"},
}


# ------------------------------------------------------------------ загрузка
def load_dataset(name, angles_limit=None, mask_angles=(0.0, 90.0),
                 average_reps=False, policy=scan_loader.BG_TWO_ABS):
    """-> (angles, freqs, exp, valid, meta, tnoise). Два пути в одном месте.

    Политика фона по умолчанию — `two_abs` (приём владельца против дрейфа:
    амплитуда как среднее двух ближайших фонов, фаза у ближайшего). Она
    измеримо лучше прочих на specac-2: избыточная невязка 29.2 против 41.7 у
    «предшествующего» прецедента, калиброванный дрейф 2.82 % против 3.31 %.
    Все три политики считаются и кладутся в артефакт — разница между ними есть
    прямая мера дрейфа.
    """
    if name in DIR_DATASETS:
        ang, freqs, exp, valid, meta = scan_loader.load_scan(
            DIR_DATASETS[name], average_reps=average_reps, policy=policy,
            mask_angles=mask_angles, angles_limit=angles_limit)
        tnoise = _tnoise_from_dir(DIR_DATASETS[name], freqs)
        meta["source"] = "scan_loader"
    else:
        import exp_builder                                          # noqa: E402
        from a3_eps_stability import noise_floor                    # noqa: E402
        ang, freqs, exp, valid, meta = exp_builder.load(
            name, angles_limit=angles_limit, mask_angles=mask_angles)
        tnoise = noise_floor(name, angles_limit)
        meta["source"] = "exp_builder"
    return ang, freqs, exp, valid, meta, tnoise


def _tnoise_from_dir(directory, freqs):
    """Шумовой пол по ВЧ-хвосту для каталожного датасета.

    Повторяет логику `a3_eps_stability.noise_floor`: полка спектра ФОНА выше
    3.0 ТГц принимается за уровень шума, tnoise = (шум/сигнал_фона)^2.
    """
    traces = scan_loader.scan_dir(directory)
    bgs = [tr for tr in traces if tr.kind == "bg"]
    spec, fc = [], None
    for tr in bgs:
        t, E = scan_loader._read(tr.path)
        dt = float(t[1] - t[0])
        spec.append(np.abs(np.fft.rfft(E)))
        fc = np.fft.rfftfreq(len(E), dt) if fc is None else fc
    bg_avg = np.mean(spec, axis=0)
    tail = np.where(fc >= 3.0)[0]
    if not len(tail):
        tail = np.where(fc >= 2.5)[0]
    nfa = float(np.mean(bg_avg[tail]))
    tn = 1.0 / ((bg_avg / nfa) ** 2)
    idx = np.array([int(np.argmin(np.abs(fc - f))) for f in freqs])
    return tn[idx]


# ------------------------------------------------------- дрейф уровня угла
def fit_drift(theta, Y, sigma2, theta0_deg, n_iter=6):
    """Комплексный множитель на трассу (HN6: мультипликативный дрейф).

    Один множитель на ВСЮ трассу (общий по частотам), поэтому он не может
    подменить модель: у той 3 комплексных коэффициента НА КАЖДЫЙ бин.
    Физически это дрейф мощности лазера и перепозиционирование между трассами.
    """
    Yc = Y.copy()
    g_tot = np.ones(Y.shape[0], dtype=complex)
    w = 1.0 / sigma2
    for _ in range(n_iter):
        fb = fl.solve_bins(theta, Yc, sigma2, theta0_deg=theta0_deg)
        M = Yc - fb.resid
        num = np.sum(w * np.conj(M) * Yc, axis=1)
        den = np.maximum(np.sum(w * np.abs(M) ** 2, axis=1), 1e-300)
        g = num / den
        g = g / np.abs(np.mean(g))            # общий масштаб не трогаем
        Yc = Yc / g[:, None]
        g_tot = g_tot * g
    return Yc, g_tot


def noise_from_reps(theta, Y, sigma2):
    """Прямая оценка шума УРОВНЯ УГЛА по повторным трассам.

    Для каждого угла с >= 2 повторами дисперсия между повторами оценивает
    ПОЛНЫЙ шум трассы (включая дрейф лазера и перепозиционирование), тогда как
    tnoise видит только ВЧ-полку спектра. Отношение — во сколько раз модель
    шума занижена. До `specac/exp_data_2` этой оценки в проекте не было:
    A6/HN13 выносили вердикт по консервативной верхней границе.
    """
    uniq, cnt = np.unique(theta, return_counts=True)
    rep_ang = uniq[cnt >= 2]
    if not len(rep_ang):
        return None
    var_emp, var_mod, n_pairs = [], [], 0
    for a in rep_ang:
        rows = np.where(theta == a)[0]
        block = Y[rows]
        n = len(rows)
        v = np.sum(np.abs(block - block.mean(axis=0)) ** 2, axis=0) / (n - 1)
        var_emp.append(v)
        var_mod.append(np.mean(sigma2[rows], axis=0))
        n_pairs += n
    var_emp = np.mean(var_emp, axis=0)
    var_mod = np.mean(var_mod, axis=0)
    ratio = var_emp / np.maximum(var_mod, 1e-300)
    return {"n_angles_with_reps": int(len(rep_ang)), "n_traces": int(n_pairs),
            "var_emp": var_emp, "var_model": var_mod,
            "ratio_median": float(np.median(ratio)),
            "ratio_per_bin": ratio}


# ------------------------------------------------------- статистика остатков
def residual_stats(resid, sigma2):
    """Нормальность и структура остатков — гардрейл CLAUDE.md.

    Вывод о качестве подгонки делается по СТАТИСТИКЕ остатков, а не по одному
    числу невязки. Re и Im имеют дисперсию sigma^2/2 каждая.
    """
    from scipy import stats
    r = np.concatenate([np.real(resid).ravel(), np.imag(resid).ravel()])
    s = np.concatenate([np.sqrt(sigma2 / 2.0).ravel()] * 2)
    z = r / np.maximum(s, 1e-300)
    z = z[np.isfinite(z)]
    sub = z if len(z) <= 5000 else z[np.linspace(0, len(z) - 1, 5000).astype(int)]
    out = {"n": int(len(z)), "mean": float(np.mean(z)), "std": float(np.std(z)),
           "skew": float(stats.skew(z)), "kurtosis": float(stats.kurtosis(z))}
    try:
        out["shapiro_p"] = float(stats.shapiro(sub).pvalue)
    except Exception:                                    # noqa: BLE001
        out["shapiro_p"] = float("nan")
    jb = stats.jarque_bera(z)
    out["jarque_bera_p"] = float(jb.pvalue)
    return out


# ------------------------------------------------------------------ ветки
def run_field(theta, freqs, Y, tnoise, kappa=fl.KAPPA_ABS_TO_VAR, drift=True,
              eps_mult=None):
    """Полевая ветка целиком. -> словарь результатов.

    `eps_mult=None` — калибровать мультипликативную компоненту шума по данным
    (chi^2/dof = 1). Веса тени относительно яркой зоны у чисто аддитивной и у
    смешанной модели различаются на два порядка, и от этого зависит theta0,
    поэтому назначать модель шума вместо калибровки нельзя.
    """
    if eps_mult is None:
        eps_mult, _ = fl.calibrate_noise(theta, Y, tnoise, kappa=kappa)
    s2 = fl.noise_variance(Y, tnoise, kappa=kappa, eps_mult=eps_mult)

    free = fl.solve_bins(theta, Y, s2, order=2)
    wbin = 1.0 / np.maximum(np.median(s2, axis=0), 1e-300)
    r1 = fl.theta0_rank1(free.z, w=wbin)
    prof = fl.theta0_profile(theta, Y, s2)

    theta0 = prof["theta0_deg"]
    bound = fl.solve_bins(theta, Y, s2, theta0_deg=theta0)
    Tp, Tl, Cov2 = fl.channels_from_z(bound)
    r, var_lnr = fl.ratio_and_var(Tp, Tl, Cov2)

    snr_par = np.abs(Tl) ** 2 / np.maximum(Cov2[:, 1, 1], 1e-300)
    good = snr_par > 9.0                                    # SNR >= 3 по тени
    lr = fl.logratio_regression(freqs, r, var_lnr, valid=good)

    o4 = fl.order4_field_test(theta, Y, s2)
    geo = fl.extinction_geometry(Tp, Tl)
    cs = np.abs(np.cos(np.angle(Tp) - np.angle(Tl)))

    out = {
        "eps_mult": float(eps_mult),
        "theta0_deg": float(theta0),
        "sigma_theta0_deg": float(prof["sigma_theta0_deg"]),
        "theta0_rank1_deg": float(r1["theta0_deg"]),
        "lambda2_over_lambda1": float(r1["lambda2_over_lambda1"]),
        "imag_leak_median": float(np.median(np.abs(r1["imag_leak"]))),
        "chi2_red": float(np.mean(bound.chi2) / bound.dof),
        "dof_per_bin": int(bound.dof),
        "T_perp": Tp, "T_par": Tl, "Cov2": Cov2, "ratio": r,
        "eta_nu": np.abs(r),
        "eta_broadband": float(np.sqrt(np.sum(np.abs(Tl) ** 2)
                                       / np.sum(np.abs(Tp) ** 2))),
        "dphi_nu": -np.angle(r),
        "n_cs_violations": int(np.sum(cs > 1.0)),
        "n_bins": int(len(freqs)),
        "n_bins_snr": int(good.sum()),
        "order4_per_dof": float(o4["mean_per_dof"]),
        "rss_perp_total": float(np.sum(bound.chi2)),
        "inner_minimum_frac": float(np.mean(geo["inner"])),
        "dtheta_min_median_deg": float(np.median(geo["dtheta_min_deg"])),
        "leakage": {k: v for k, v in lr.items()
                    if k in ("eta0", "eta_exp", "sigma_eta_exp", "delta0_rad",
                             "tau_dphi_ps", "sigma_tau_dphi_ps", "tau_leak_m5_ps", "n_used",
                             "chi2_red_amp", "chi2_red_phase")},
        "unwrap": lr.get("unwrap", {}),
        "resid_stats": residual_stats(bound.resid, s2),
        "chi2_profile": {"grid_deg": prof["grid_deg"], "chi2": prof["chi2_curve"]},
    }

    if drift:
        Yc, g = fit_drift(theta, Y, s2, theta0)
        fb2 = fl.solve_bins(theta, Yc, s2, theta0_deg=theta0)
        o4c = fl.order4_field_test(theta, Yc, s2)
        r1c = fl.theta0_rank1(fl.solve_bins(theta, Yc, s2, order=2).z, w=wbin)
        out["drift"] = {
            "abs_rms": float(np.std(np.abs(g))),
            "arg_rms_deg": float(np.std(np.degrees(np.angle(g)))),
            "chi2_red_after": float(np.mean(fb2.chi2) / fb2.dof),
            "order4_per_dof_after": float(o4c["mean_per_dof"]),
            "theta0_after_deg": float(r1c["theta0_deg"]),
            "resid_stats_after": residual_stats(fb2.resid, s2),
        }
    return out


def run_power(theta, freqs, Y):
    """Мощностная ветка (аппарат A0) — на тех же данных, рецептом A."""
    P = np.abs(Y) ** 2
    per_nu = []
    for j in range(P.shape[1]):
        col = P[:, j]
        p, _, info = ml.fit_harmonics(theta, col,
                                      sigma=ml._weights_for(col, "relative"),
                                      order=4)
        ch = ml.channel_params(p)
        per_nu.append(ch)

    U = P.sum(axis=1)
    p_s, cov_s, info_s = ml.fit_harmonics(theta, U,
                                          sigma=ml._weights_for(U, "relative"),
                                          order=4)
    ch_s = ml.channel_params(p_s, cov_s)
    return {
        "theta0_deg": float(ch_s["theta0_deg_from_h2"]),
        "sigma_theta0_deg": float(ml.theta0_error(p_s, cov_s, info_s["chi2_red"])),
        "A4_over_A2": float(ch_s["harmonic4_over_harmonic2"]),
        "eta_broadband": float(ch_s["eta_amplitude"]),
        "eta_nu": np.array([c["eta_amplitude"] for c in per_nu]),
        "cos_dphi_nu": np.array([c["cos_dphi"] for c in per_nu]),
        "n_cs_violations": int(sum(0 if c["cauchy_schwarz_ok"] else 1
                                   for c in per_nu)),
        "n_neg_I_par": int(sum(1 for c in per_nu if c["I_par"] <= 0)),
        "n_bins": int(P.shape[1]),
    }


def compare_bg_policies(name, kappa=fl.KAPPA_ABS_TO_VAR):
    """Три политики выбора фона на одних данных. -> список словарей.

    Разница между политиками — прямая мера дрейфа между фоновыми трассами, а
    заодно проверка приёма владельца («усреднять два ближайших фона по
    амплитуде»). Мощностной метод к выбору политики почти безразличен, полевой
    нет: он использует фазу фона.
    """
    out = []
    for pol in (scan_loader.BG_PRECEDING, scan_loader.BG_NEAREST,
                scan_loader.BG_TWO_ABS):
        ang, freqs, exp, valid, meta, tn = load_dataset(name, policy=pol)
        wm = meta["water_mask"]
        Y, f, tnv = exp[:, wm], freqs[wm], tn[wm]
        s2_0 = fl.noise_variance(Y, tnv, kappa=kappa)
        t0_0 = fl.theta0_rank1(fl.solve_bins(ang, Y, s2_0, order=2).z)["theta0_deg"]
        fb0 = fl.solve_bins(ang, Y, s2_0, theta0_deg=t0_0)
        eps, _ = fl.calibrate_noise(ang, Y, tnv, kappa=kappa)
        s2 = fl.noise_variance(Y, tnv, kappa=kappa, eps_mult=eps)
        prof = fl.theta0_profile(ang, Y, s2)
        fb = fl.solve_bins(ang, Y, s2, theta0_deg=prof["theta0_deg"])
        Tp, Tl, _ = fl.channels_from_z(fb)
        out.append({
            "policy": pol, "eps_mult": float(eps),
            "chi2_red_additive_only": float(np.mean(fb0.chi2) / fb0.dof),
            "theta0_deg": float(prof["theta0_deg"]),
            "sigma_theta0_deg": float(prof["sigma_theta0_deg"]),
            "eta_broadband": float(np.sqrt(np.sum(np.abs(Tl) ** 2)
                                           / np.sum(np.abs(Tp) ** 2))),
        })
    return out


# ------------------------------------------------------------------ прогон
def analyse(name, kappa=fl.KAPPA_ABS_TO_VAR, verbose=True):
    ang, freqs, exp, valid, meta, tnoise = load_dataset(name)
    wm = meta["water_mask"]
    Y, f, tn = exp[:, wm], freqs[wm], tnoise[wm]

    field = run_field(ang, f, Y, tn, kappa=kappa)
    power = run_power(ang, f, Y)

    # Сравнивать повторы надо с ОТКАЛИБРОВАННОЙ моделью шума, иначе отношение
    # меряет лишь то, что аддитивной части мало (и так известно).
    s2_add = fl.noise_variance(Y, tn, kappa=kappa)
    s2_cal = fl.noise_variance(Y, tn, kappa=kappa, eps_mult=field["eps_mult"])
    reps = noise_from_reps(ang, Y, s2_cal)
    reps_add = noise_from_reps(ang, Y, s2_add)
    if reps is not None and reps_add is not None:
        reps["ratio_median_additive_only"] = reps_add["ratio_median"]

    res = {
        "dataset": name,
        "source": meta["source"],
        "n_traces": int(len(ang)),
        "n_distinct_angles": int(len(np.unique(ang))),
        "angle_min": float(np.min(ang)), "angle_max": float(np.max(ang)),
        "n_bins": int(wm.sum()),
        "f_min_thz": float(f[0]), "f_max_thz": float(f[-1]),
        "kappa": float(kappa),
        "field": field, "power": power,
        "reference": REFERENCE.get(name, {}),
    }
    if reps is not None:
        res["noise_from_reps"] = {
            "n_angles_with_reps": reps["n_angles_with_reps"],
            "n_traces": reps["n_traces"],
            "ratio_median": reps["ratio_median"],
            "ratio_median_additive_only": reps.get("ratio_median_additive_only"),
        }
    if name in DIR_DATASETS:
        res["bg_policy_scan"] = compare_bg_policies(name, kappa=kappa)

    if verbose:
        _report(res, field, power, reps)
    return res, (ang, f, Y, field, power)


def _report(res, field, power, reps):
    name = res["dataset"]
    print("=" * 78)
    print("A22 %s — трасс %d, различных углов %d (%.0f..%.0f), бинов %d (%.2f-%.2f ТГц)"
          % (name, res["n_traces"], res["n_distinct_angles"], res["angle_min"],
             res["angle_max"], res["n_bins"], res["f_min_thz"], res["f_max_thz"]))
    print("-" * 78)
    ref = res["reference"]
    t0ref = ("   [A13: %+.4f +- %.4f]" % (ref["theta0_deg"], ref["sigma_theta0_deg"])
             if "theta0_deg" in ref else "")
    print("theta0:  поле %+.4f +- %.4f   мощность %+.4f +- %.4f%s"
          % (field["theta0_deg"], field["sigma_theta0_deg"],
             power["theta0_deg"], power["sigma_theta0_deg"], t0ref))
    print("         rank-1 %+.4f, lambda2/lambda1 = %.2e (0 = единый theta0, тест HN13)"
          % (field["theta0_rank1_deg"], field["lambda2_over_lambda1"]))
    print("eta:     поле %.5f   мощность %.5f%s"
          % (field["eta_broadband"], power["eta_broadband"],
             ("   [Парсеваль A0: %.5f]" % ref["eta_parseval"])
             if "eta_parseval" in ref else ""))
    print("нарушений |cos dphi|<=1:  поле %d/%d   мощность %d/%d%s"
          % (field["n_cs_violations"], field["n_bins"],
             power["n_cs_violations"], power["n_bins"],
             ("   [A0 ранее: %s]" % ref["cs_violations"])
             if "cs_violations" in ref else ""))
    print("         |t_par|^2 < 0 в мощности: %d бинов" % power["n_neg_I_par"])
    lk = field["leakage"]
    print("утечка без модели: eta0 = %.5g, показатель a = %.3f +- %.3f"
          % (lk.get("eta0", np.nan), lk.get("eta_exp", np.nan),
             lk.get("sigma_eta_exp", np.nan)))
    print("tau_leak = %+.4f +- %.4f пс (бинов с SNR>=3: %d)%s"
          % (lk.get("tau_leak_m5_ps", np.nan), lk.get("sigma_tau_dphi_ps", np.nan),
             field["n_bins_snr"],
             ("   [ранее: %.3f пс]" % ref["tau_leak_ps"])
             if "tau_leak_ps" in ref else ""))
    uw = field.get("unwrap", {})
    if uw:
        print("         развёртка: кусков %d, сшивок %d, min dchi2 = %.1f%s"
              % (uw.get("n_segments", 0), uw.get("n_joins", 0),
                 uw.get("min_delta_chi2", np.nan),
                 "  [НЕОДНОЗНАЧНО]" if uw.get("ambiguous") else ""))
    print("-" * 78)
    print("chi2/dof = %.2f (dof %d на бин); 4-я гармоника в поле = %.1f на dof"
          % (field["chi2_red"], field["dof_per_bin"], field["order4_per_dof"]))
    d = field.get("drift")
    if d:
        print("дрейф уровня угла: |g| rms %.4f, arg rms %.3f град"
              % (d["abs_rms"], d["arg_rms_deg"]))
        print("  после снятия: chi2/dof %.2f (было %.2f), 4-я гармоника %.1f (было %.1f)"
              % (d["chi2_red_after"], field["chi2_red"],
                 d["order4_per_dof_after"], field["order4_per_dof"]))
    if reps is not None:
        print("шум по ПОВТОРАМ (углов с повторами %d, трасс %d): реальная "
              "дисперсия / модельная = %.2f"
              % (reps["n_angles_with_reps"], reps["n_traces"],
                 reps["ratio_median"]))
        print("         та же величина для ЧИСТО АДДИТИВНОЙ модели = %.1f "
              "-> мультипликативная часть обязательна"
              % reps.get("ratio_median_additive_only", float("nan")))
    bp = res.get("bg_policy_scan")
    if bp:
        print("политика фона (мера дрейфа между фонами):")
        for row in bp:
            print("   %-10s eps_mult %.4f, chi2/dof без мультипл. %6.2f, "
                  "theta0 %+.4f, eta %.5f"
                  % (row["policy"], row["eps_mult"],
                     row["chi2_red_additive_only"], row["theta0_deg"],
                     row["eta_broadband"]))
    rs = field["resid_stats"]
    print("остатки: mean %+.3f, std %.3f, асимметрия %+.3f, эксцесс %+.3f; "
          "Шапиро p = %.2e, Харке-Бера p = %.2e"
          % (rs["mean"], rs["std"], rs["skew"], rs["kurtosis"],
             rs["shapiro_p"], rs["jarque_bera_p"]))
    print("минимум |E| внутренний в %.0f %% бинов, медианный сдвиг %.2f град от theta0+90"
          % (100 * field["inner_minimum_frac"],
             90.0 - field["dtheta_min_median_deg"]))


# ------------------------------------------------------------------ графики
def make_plots(name, ang, f, Y, field, power):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    RESULTS.mkdir(parents=True, exist_ok=True)
    Tp, Tl = field["T_perp"], field["T_par"]
    err = np.sqrt(np.maximum(field["Cov2"][:, 0, 0], 0.0))
    errl = np.sqrt(np.maximum(field["Cov2"][:, 1, 1], 0.0))

    fig, ax = plt.subplots(2, 2, figsize=(12, 8))
    a = ax[0, 0]
    a.errorbar(f, 20 * np.log10(np.abs(Tp)), yerr=8.686 * err / np.abs(Tp),
               fmt="-", lw=1.2, label=r"$|T_\perp|$")
    a.errorbar(f, 20 * np.log10(np.abs(Tl)), yerr=8.686 * errl / np.abs(Tl),
               fmt="-", lw=1.2, label=r"$|T_\parallel|$")
    a.set_xlabel("частота, ТГц"); a.set_ylabel("дБ")
    a.set_title("спектры каналов (поле, модельно-независимо)")
    a.legend(); a.grid(alpha=.3)

    a = ax[0, 1]
    a.plot(f, field["eta_nu"], lw=1.2, label="поле")
    a.plot(f, power["eta_nu"], lw=1.0, alpha=.75, label="мощность (A0)")
    a.axhline(field["eta_broadband"], ls="--", c="k", lw=.8,
              label="поле, широкополосно")
    a.set_xlabel("частота, ТГц"); a.set_ylabel(r"$\eta = |T_\parallel/T_\perp|$")
    a.set_yscale("log"); a.set_title("утечка (HN2) без модели")
    a.legend(); a.grid(alpha=.3)

    a = ax[1, 0]
    a.plot(f, np.degrees(field["dphi_nu"]), lw=1.2, label=r"поле: $\Delta\varphi$")
    cs = np.clip(power["cos_dphi_nu"], -1, 1)
    a.plot(f, np.degrees(np.arccos(cs)), lw=1.0, alpha=.7,
           label=r"мощность: $\arccos$ (ветвь неизвестна)")
    bad = np.abs(power["cos_dphi_nu"]) > 1
    if bad.any():
        a.plot(f[bad], np.zeros(bad.sum()), "rx", ms=6,
               label="мощность: |cos| > 1")
    a.set_xlabel("частота, ТГц"); a.set_ylabel("градусы")
    a.set_title("относительная фаза каналов")
    a.legend(fontsize=8); a.grid(alpha=.3)

    a = ax[1, 1]
    g = field["chi2_profile"]
    a.plot(g["grid_deg"], g["chi2"] - np.min(g["chi2"]), lw=1.2)
    a.axvline(field["theta0_deg"], c="r", ls="--", lw=.8,
              label=r"$\theta_0$ = %+.4f" % field["theta0_deg"])
    a.set_ylim(0, 50); a.set_xlim(field["theta0_deg"] - 8, field["theta0_deg"] + 8)
    a.set_xlabel(r"$\theta_0$, град"); a.set_ylabel(r"$\Delta\chi^2$")
    a.set_title("профиль правдоподобия по угловому офсету")
    a.legend(); a.grid(alpha=.3)

    fig.suptitle("A22 %s — восстановление каналов линеаризацией в поле" % name)
    fig.tight_layout()
    out = RESULTS / ("a22_spectra_%s.png" % name)
    fig.savefig(out, dpi=110); plt.close(fig)
    return out


# ------------------------------------------------------------------ main
def _jsonable(o):
    if isinstance(o, dict):
        return {k: _jsonable(v) for k, v in o.items()}
    if isinstance(o, (list, tuple)):
        return [_jsonable(v) for v in o]
    if isinstance(o, np.ndarray):
        if np.iscomplexobj(o):
            return {"re": o.real.tolist(), "im": o.imag.tolist()}
        return o.tolist()
    if isinstance(o, (np.floating, np.integer)):
        return o.item()
    if isinstance(o, (np.bool_, bool)):
        return bool(o)
    return o


def main():
    ap = argparse.ArgumentParser(description="A22: поле против мощности")
    ap.add_argument("--dataset", type=str, default=None)
    ap.add_argument("--all", action="store_true",
                    help="test_grid_33_11 и specac-2")
    ap.add_argument("--kappa", type=float, default=fl.KAPPA_ABS_TO_VAR)
    ap.add_argument("--no-plots", action="store_true")
    args = ap.parse_args()

    names = (["test_grid_33_11", "specac-2"] if args.all
             else ([args.dataset] if args.dataset else []))
    if not names:
        ap.print_help()
        return 0

    RESULTS.mkdir(parents=True, exist_ok=True)
    summary = []
    for name in names:
        res, (ang, f, Y, field, power) = analyse(name, kappa=args.kappa)
        out = RESULTS / ("a22_field_%s.json" % name)
        out.write_text(json.dumps(_jsonable(res), ensure_ascii=False, indent=1),
                       encoding="utf-8")
        print("артефакт: %s" % out)
        if not args.no_plots:
            print("график:   %s" % make_plots(name, ang, f, Y, field, power))
        summary.append({k: res[k] for k in ("dataset", "n_traces", "n_bins")}
                       | {"theta0_field": field["theta0_deg"],
                          "theta0_power": power["theta0_deg"],
                          "eta_field": field["eta_broadband"],
                          "eta_power": power["eta_broadband"],
                          "cs_viol_field": field["n_cs_violations"],
                          "cs_viol_power": power["n_cs_violations"],
                          "chi2_red": field["chi2_red"],
                          "tau_leak_m5_ps": field["leakage"].get("tau_leak_m5_ps")})
        print()

    (RESULTS / "a22_field_summary.json").write_text(
        json.dumps(_jsonable(summary), ensure_ascii=False, indent=1),
        encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
