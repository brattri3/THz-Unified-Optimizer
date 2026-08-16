"""A24 — почему на приборе аттенюатора проваливается тест ранга 1: это ВТОРОЙ элемент.

ВОПРОС. Этап 1 (A22) раскладывает поле по углу в базисе `[1, cos2t, sin2t]` и
проверяет, что матрица `R = [[|z_c|^2, Re(conj(z_c) z_s)], [., |z_s|^2]]` имеет
ранг 1 — это условие существования ЕДИНОГО вещественного `theta0`. На чистых
образцах `lambda2/lambda1` = 1.9e-05 и 2.7e-05, на приборе ATT-11-16 — 2.11e-03
и 7.62e-04, на два порядка выше. Гипотеза: дело в том, что элементов ДВА, а
разложение выведено для одного вращающегося.

ОТВЕТ КОРОТКО. Да, и это считается в замкнутом виде. Линейный базис 2-го порядка
для двух элементов остаётся ТОЧНЫМ — ломается не он, а условие ранга 1.

ВЫВОД. Пусть вращается элемент 2, элемент 1 стоит на `phi1`; `d` — ось детектора
(угол `psi_d`), `E_in` — поляризация излучателя (угол `psi_e`). Матрица Джонса
поляризатора раскладывается на скаляр и бесследовую часть:

    J(phi) = t_ * I + dt * M(phi),   t_ = (a+b)/2,  dt = (a-b)/2,
    M(phi) = [[cos2phi, sin2phi], [sin2phi, -cos2phi]] = cos2phi*s3 + sin2phi*s1

где `a = t_perp`, `b = t_par`. Тогда

    E(theta) = d^T [t_2 I + dt_2 M(theta)] v,      v = J_1(phi1) E_in
    z0 = t_2 (d^T v),   z_c = dt_2 (d^T s3 v),   z_s = dt_2 (d^T s1 v)

то есть поле ЛИНЕЙНО по трём комплексным коэффициентам — базис order=2 точен и
здесь. Дальше, обозначив `g3 = d^T s3 E_in`, `h3 = d^T s3 M(phi1) E_in` и т.д.
(все четыре — ВЕЩЕСТВЕННЫЕ):

    z_c/dt_2 = t_1 g3 + dt_1 h3 ,   z_s/dt_2 = t_1 g1 + dt_1 h1
    Im(conj(z_c) z_s) = |dt_2|^2 * Im(conj(t_1) dt_1) * (g3 h1 - g1 h3)

Оба множителя считаются до конца. Спектральный:

    Im(conj(t_) dt) = (1/2) Im(conj(b) a) = (1/2) |a| |b| sin(dphi),
    dphi = arg(t_perp) - arg(t_par)

Геометрический — единичные векторы: `(g3, g1) = (cos S, sin S)` при
`S = psi_d + psi_e`, `(h3, h1) = (cos L, sin L)` при `L = psi_d - psi_e + 2 phi1`,
откуда `g3 h1 - g1 h3 = sin(L - S) = sin(2 phi1 - 2 psi_e)`. Итого

    Im(conj(z_c) z_s) = (1/2) |dt_2|^2 * |a_1| |b_1| * sin(dphi_1)
                        * sin(2*phi1 - 2*psi_e)                          (*)

Знаменатель тоже замкнут: `|z_c|^2 + |z_s|^2 = |dt_2|^2 * I_1`, где
`I_1 = |a_1|^2 cos^2(phi1 - psi_e) + |b_1|^2 sin^2(phi1 - psi_e)` — это в
точности закон Малюса для ПЕРВОГО элемента.

И, наконец, для вещественной симметричной 2x2 матрицы `R` выполняется ТОЖДЕСТВО

    det R = |z_c|^2 |z_s|^2 - Re(conj(z_c) z_s)^2 = Im(conj(z_c) z_s)^2

поэтому при почти-ранге-1 (`lambda1 ~ trace`)

    lambda2/lambda1 ~ det R / trace^2
                    = [ (1/2) |a_1| |b_1| sin(dphi_1) sin(2 phi1 - 2 psi_e) / I_1 ]^2  (**)

ЧТО ЭТО ДАЁТ, КРОМЕ ОБЪЯСНЕНИЯ.

1. Нарушение ранга 1 обнуляется при `phi1 = psi_e` и `phi1 = psi_e + 90` — то
   есть когда ось пропускания НЕПОДВИЖНОГО элемента параллельна или
   перпендикулярна поляризации излучателя. Это проверяемое предсказание и
   рецепт съёмки: так снятый двухэлементный прибор разложится как одиночный.
2. Величина управляется утечкой и фазой ЭЛЕМЕНТА, а не системы. У прибора
   измеренная широкополосная `eta` = 0.0129 — это системная величина (два
   скрещенных элемента); у элемента она примерно её корень, ~0.11, и именно она
   стоит в (**).
3. `lambda2/lambda1` обязана РАСТИ с частотой, потому что `dphi(nu)` растёт.
   Это проверяется прямо на данных прибора, без всякой симуляции.

ЕСЛИ ВРАЩАЕТСЯ ПЕРВЫЙ элемент (`att-11-16-s1`), формула та же по симметрии:
`E = d^T J_2(phi2) J_1(theta) E_in = (J_2 d)^T J_1(theta) E_in`, то есть роли `d`
и `E_in` меняются местами, и геометрический множитель становится
`sin(2*phi2 - 2*psi_d)`, а спектральный берётся у ВТОРОГО элемента.

СЛЕДСТВИЕ ДЛЯ МЕТОДА. Разложение `cos^2(t-t0) T_perp + sin^2(t-t0) T_par` для
двух элементов НЕПРИМЕНИМО: вещественного `theta0` не существует. Но достаточной
статистикой остаётся тройка `(z0, z_c, z_s)` — базис-то точен. Значит этап 2 для
прибора должен подгонять модель к ТРЁМ комплексным коэффициентам, а не к двум
каналам. Это ровно то, что умеет `model_2wgp.forward_field`.

Зона A. `model_2wgp`, `field_linearization` — свои; `model_core` только читается.

Запуск ИЗ КОРНЯ репозитория:
    .venv/Scripts/python.exe research/two_wgp/a24_two_element_rank.py --selftest
    .venv/Scripts/python.exe research/two_wgp/a24_two_element_rank.py --device
    .venv/Scripts/python.exe research/two_wgp/a24_two_element_rank.py --data
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
for p in (str(REPO), str(HERE), str(REPO / "research" / "experiments")):
    if p not in sys.path:
        sys.path.insert(0, p)

import field_linearization as fl                                # noqa: E402
import model_2wgp as m2                                         # noqa: E402
import a22_field_vs_power as a22                                # noqa: E402

OUT = REPO / "research" / "results" / "two_wgp"
OUT.mkdir(parents=True, exist_ok=True)

#: Паспорт прибора (`attenuator_app/passports/ATT-11-16-CA85_02721.json`, ЧТЕНИЕ).
PASSPORT = {"P_um": 15.5, "D_eff_um": 5.67, "D_phys_um": 11.0,
            "loss_factor": 0.255, "gamma": 1.58}

#: Измерено этапом 1 на реальных данных прибора (A22).
MEASURED = {
    "att-11-16-356": {"lambda2_over_lambda1": 2.11e-03, "eta_broadband": 0.01287,
                      "rotating": 2, "note": "вращается ВТОРОЙ WGP (у детектора)"},
    "att-11-16-s1": {"lambda2_over_lambda1": 7.62e-04, "eta_broadband": 0.01400,
                     "rotating": 1, "note": "вращается ПЕРВЫЙ WGP"},
    "test_grid_33_11": {"lambda2_over_lambda1": 1.85e-05, "eta_broadband": 0.05314,
                        "rotating": None, "note": "один элемент — контроль"},
    "specac-2": {"lambda2_over_lambda1": 2.71e-05, "eta_broadband": 0.03826,
                 "rotating": None, "note": "один элемент — контроль"},
}


# ------------------------------------------------------------------ аналитика
def imag_leak_analytic(a1, b1, dt2_abs2, phi_fixed_deg, psi_ref_deg):
    """Формула (*): Im(conj(z_c) z_s) для двух элементов.

    `a1`, `b1` — комплексные t_perp/t_par НЕПОДВИЖНОГО элемента (массивы по ν);
    `dt2_abs2` — |dt|^2 вращающегося (общий множитель, в отношение не входит);
    `phi_fixed_deg` — угол оси пропускания неподвижного элемента;
    `psi_ref_deg` — поляризация излучателя (если вращается второй элемент) либо
    ось детектора (если вращается первый).
    """
    geo = np.sin(2 * np.deg2rad(float(phi_fixed_deg) - float(psi_ref_deg)))
    spec = 0.5 * np.imag(np.conj(b1) * a1)
    return dt2_abs2 * spec * geo


def malus_intensity(a1, b1, phi_fixed_deg, psi_ref_deg):
    """`I_1` из знаменателя: закон Малюса для неподвижного элемента."""
    d = np.deg2rad(float(phi_fixed_deg) - float(psi_ref_deg))
    return np.abs(a1) ** 2 * np.cos(d) ** 2 + np.abs(b1) ** 2 * np.sin(d) ** 2


def ratio_analytic(a1, b1, phi_fixed_deg, psi_ref_deg):
    """Формула (**): lambda2/lambda1, побинно. Общий множитель dt2 сокращается."""
    num = imag_leak_analytic(a1, b1, 1.0, phi_fixed_deg, psi_ref_deg)
    den = malus_intensity(a1, b1, phi_fixed_deg, psi_ref_deg)
    return (num / np.maximum(den, 1e-300)) ** 2


# ------------------------------------------------------------------ симуляция
def simulate(theta_deg, freqs, phi_fixed_deg, *, rotating=2, psi_e=0.0, psi_d=0.0,
             wgp=None, wgp_fixed=None):
    """Поле двухэлементного прибора на угловой сетке. -> (n_ang, n_bin) комплексное."""
    w = wgp or m2.WGP(p_um=PASSPORT["P_um"], d_um=PASSPORT["D_eff_um"])
    wf = wgp_fixed or w
    th = np.asarray(theta_deg, dtype=float)
    if rotating == 2:
        E = m2.forward_field(float(phi_fixed_deg), th, freqs, wf, w,
                             e_in_deg=psi_e, det_deg=psi_d)
    else:
        E = m2.forward_field(th, float(phi_fixed_deg), freqs, w, wf,
                             e_in_deg=psi_e, det_deg=psi_d)
    return np.asarray(E)


def rank_of(theta_deg, Y, w=None):
    """Прогон этапа 1 без модели шума. -> словарь theta0_rank1."""
    return fl.theta0_rank1(fl.solve_bins(theta_deg, Y, order=2).z, w=w)


# ------------------------------------------------------------------- прогоны
def device_scan(n_phi=181, verbose=True):
    """Скан по углу НЕПОДВИЖНОГО элемента: где ранг 1 восстанавливается."""
    freqs = np.linspace(0.2, 1.5, 60)
    theta = np.arange(-90.0, 90.0, 5.0)
    w = m2.WGP(p_um=PASSPORT["P_um"], d_um=PASSPORT["D_eff_um"])
    a1, b1 = w.ab(freqs)

    phis = np.linspace(0.0, 180.0, n_phi)
    num, ana = [], []
    for phi in phis:
        Y = simulate(theta, freqs, phi)
        num.append(rank_of(theta, Y)["lambda2_over_lambda1"])
        ana.append(float(np.mean(ratio_analytic(a1, b1, phi, 0.0))))
    num, ana = np.array(num), np.array(ana)

    best = int(np.argmax(num))
    zeros = phis[np.r_[True, num[1:] < num[:-1]] & np.r_[num[:-1] < num[1:], True]]
    out = {"phi_deg": phis, "lambda_ratio_numeric": num,
           "lambda_ratio_analytic_binmean": ana,
           "max_at_phi_deg": float(phis[best]), "max_value": float(num[best]),
           "min_value": float(num.min()), "min_at_phi_deg": float(phis[int(np.argmin(num))]),
           "zeros_phi_deg": [float(z) for z in zeros]}
    if verbose:
        print("Скан по углу неподвижного элемента (psi_e = 0):")
        for phi in (0.0, 15.0, 30.0, 45.0, 60.0, 75.0, 90.0):
            i = int(np.argmin(np.abs(phis - phi)))
            print("  phi1 = %5.1f град -> lambda2/lambda1 = %.3e" % (phi, num[i]))
        print("  максимум %.3e при phi1 = %.1f, минимум %.3e при phi1 = %.1f"
              % (out["max_value"], out["max_at_phi_deg"], out["min_value"],
                 out["min_at_phi_deg"]))
    return out


def match_measured(dataset, verbose=True):
    """Какой угол неподвижного элемента воспроизводит ИЗМЕРЕННОЕ lambda2/lambda1."""
    meas = MEASURED[dataset]["lambda2_over_lambda1"]
    freqs = np.linspace(0.2, 1.5, 60)
    theta = np.arange(-90.0, 90.0, 5.0)
    rot = MEASURED[dataset]["rotating"]

    phis = np.linspace(0.0, 90.0, 181)
    vals = np.array([rank_of(theta, simulate(theta, freqs, p, rotating=rot))
                     ["lambda2_over_lambda1"] for p in phis])
    reach = bool(vals.max() >= meas)
    i = int(np.argmin(np.abs(vals - meas)))
    out = {"dataset": dataset, "measured": meas, "reachable": reach,
           "model_max": float(vals.max()),
           "phi_matching_deg": float(phis[i]) if reach else None,
           "value_at_match": float(vals[i])}
    if verbose:
        print("%s: измерено %.2e, модель прибора достигает максимум %.2e%s"
              % (dataset, meas, vals.max(),
                 ("; совпадает при phi = %.1f град" % phis[i]) if reach
                 else " — НЕ дотягивает"))
    return out


def joint_match(dataset="att-11-16-356", d_fix_grid=(5.67, 4.0, 3.0, 2.0, 1.5),
                phi_grid=(2.0, 5.0, 10.0, 15.0, 20.0, 30.0, 45.0), verbose=True):
    """РЕШАЮЩАЯ проверка: один угол обязан объяснить И lambda2/lambda1, И eta.

    Подогнать одно число под одну свободную ручку можно всегда. Механизм двух
    элементов имеет право называться объяснением только если ОДНА конфигурация
    (угол неподвижного элемента + его геометрия) воспроизводит обе измеренные
    величины сразу: нарушение ранга 1 и широкополосную утечку.

    Свободных ручек две (`phi1`, `D_eff` неподвижного элемента), измеренных
    величин две — то есть проверка не тавтологична только потому, что решения
    может не существовать.
    """
    meas = MEASURED[dataset]
    freqs = np.linspace(0.2, 1.5, 60)
    theta = np.arange(-90.0, 90.0, 5.0)
    rot = m2.WGP(p_um=PASSPORT["P_um"], d_um=PASSPORT["D_eff_um"])

    rows = []
    for dfix in d_fix_grid:
        fix = m2.WGP(p_um=PASSPORT["P_um"], d_um=dfix)
        af, bf = fix.ab(freqs)
        eta_fix = float(np.median(np.abs(bf) / np.abs(af)))
        for phi in phi_grid:
            Y = simulate(theta, freqs, phi, rotating=meas["rotating"] or 2,
                         wgp=rot, wgp_fixed=fix)
            r = rank_of(theta, Y)
            b = fl.solve_bins(theta, Y, theta0_deg=r["theta0_deg"])
            Tp, Tl, _ = fl.channels_from_z(b)
            eta = float(np.sqrt(np.sum(np.abs(Tl) ** 2) / np.sum(np.abs(Tp) ** 2)))
            rows.append({"D_fix_um": dfix, "eta_fix": eta_fix, "phi_deg": phi,
                         "lambda2_over_lambda1": r["lambda2_over_lambda1"],
                         "eta_broadband": eta,
                         "rank_ratio_to_measured": r["lambda2_over_lambda1"]
                         / meas["lambda2_over_lambda1"],
                         "eta_ratio_to_measured": eta / meas["eta_broadband"]})

    # где выполняется каждое условие по отдельности (в пределах x2)
    ok_rank = [r for r in rows if 0.5 < r["rank_ratio_to_measured"] < 2.0]
    ok_eta = [r for r in rows if 0.5 < r["eta_ratio_to_measured"] < 2.0]
    both = [r for r in rows if r in ok_rank and r in ok_eta]
    # насколько НЕ сходится: лучший компромисс по сумме логарифмов
    def _miss(r):
        return (np.log(r["rank_ratio_to_measured"]) ** 2
                + np.log(r["eta_ratio_to_measured"]) ** 2)
    best = min(rows, key=_miss)

    out = {"dataset": dataset, "measured": meas, "n_configs": len(rows),
           "n_matching_rank_only": len(ok_rank), "n_matching_eta_only": len(ok_eta),
           "n_matching_both": len(both),
           "best_compromise": best, "rows": rows}
    if verbose:
        print("Совместная проверка на %s (цель: lambda2/lambda1 = %.2e, eta = %.5f):"
              % (dataset, meas["lambda2_over_lambda1"], meas["eta_broadband"]))
        print("  конфигураций %d; воспроизводят ранг %d, eta %d, ОБА — %d"
              % (len(rows), len(ok_rank), len(ok_eta), len(both)))
        b = best
        print("  лучший компромисс: D_fix = %.2f (eta_fix = %.4f), phi = %.1f град"
              % (b["D_fix_um"], b["eta_fix"], b["phi_deg"]))
        print("    -> lambda2/lambda1 = %.3e (в %.2f раза от измеренного), "
              "eta = %.5f (в %.2f раза)"
              % (b["lambda2_over_lambda1"], b["rank_ratio_to_measured"],
                 b["eta_broadband"], b["eta_ratio_to_measured"]))
    return out


def data_frequency_test(dataset, verbose=True):
    """Проверка на РЕАЛЬНЫХ данных: растёт ли нарушение ранга 1 с частотой.

    Формула (*) даёт `Im(conj(z_c) z_s) ~ |a||b| sin(dphi)`, а `dphi` растёт с
    частотой; значит и нарушение обязано расти. Это проверка предсказания на
    данных, без всякой симуляции.
    """
    ang, freqs, exp, valid, meta, tnoise = a22.load_dataset(dataset)
    wm = meta["water_mask"]
    Y, f, tn = exp[:, wm], freqs[wm], tnoise[wm]
    eps, _ = fl.calibrate_noise(ang, Y, tn)
    s2 = fl.noise_variance(Y, tn, eps_mult=eps)
    fit = fl.solve_bins(ang, Y, s2, order=2)
    r1 = fl.theta0_rank1(fit.z)

    # побинное нарушение, нормированное на масштаб бина
    zc, zs = fit.z[1], fit.z[2]
    leak = np.imag(np.conj(zc) * zs)
    scale = np.abs(zc) ** 2 + np.abs(zs) ** 2
    rel = np.abs(leak) / np.maximum(scale, 1e-300)

    lo = f < np.median(f)
    sl = float(np.polyfit(f, rel, 1)[0])
    out = {"dataset": dataset, "n_bins": int(len(f)),
           "lambda2_over_lambda1": float(r1["lambda2_over_lambda1"]),
           "rel_leak_median": float(np.median(rel)),
           "rel_leak_low_half": float(np.median(rel[lo])),
           "rel_leak_high_half": float(np.median(rel[~lo])),
           "growth_factor_high_over_low": float(np.median(rel[~lo])
                                                / max(np.median(rel[lo]), 1e-300)),
           "slope_per_thz": sl,
           "f_min_thz": float(f[0]), "f_max_thz": float(f[-1])}
    if verbose:
        print("%-16s |Im|/масштаб: низ %.2e -> верх %.2e (в %.1f раза), наклон %+.2e 1/ТГц"
              % (dataset, out["rel_leak_low_half"], out["rel_leak_high_half"],
                 out["growth_factor_high_over_low"], sl))
    return out


# ------------------------------------------------------------------ selftest
def selftest(verbose=True):
    checks = []

    def chk(name, ok, detail):
        checks.append({"check": name, "ok": bool(ok), "detail": detail})
        if verbose:
            print(f"  [{'OK ' if ok else 'FAIL'}] {name}: {detail}")

    freqs = np.linspace(0.2, 1.5, 40)
    theta = np.arange(-90.0, 90.0, 5.0)          # строго периодично, без дубля
    w = m2.WGP(p_um=PASSPORT["P_um"], d_um=PASSPORT["D_eff_um"])
    a1, b1 = w.ab(freqs)

    # --- 1. базис order=2 ТОЧЕН и для двух элементов
    Y = simulate(theta, freqs, 33.0)
    fit = fl.solve_bins(theta, Y, order=2)
    rms = float(np.sqrt(np.mean(np.abs(fit.resid) ** 2)) / np.sqrt(np.mean(np.abs(Y) ** 2)))
    chk("базис order=2 описывает ДВА элемента точно", rms < 1e-14,
        "относит. rms невязки = %.2e" % rms)

    # --- 2. тождество det R = Im(conj(z_c) z_s)^2
    zc, zs = fit.z[1], fit.z[2]
    R11, R22 = np.abs(zc) ** 2, np.abs(zs) ** 2
    R12 = np.real(np.conj(zc) * zs)
    det = R11 * R22 - R12 ** 2
    im2 = np.imag(np.conj(zc) * zs) ** 2
    rel = float(np.max(np.abs(det - im2)) / max(np.max(im2), 1e-300))
    chk("тождество det R = Im(conj(z_c) z_s)^2", rel < 1e-9,
        "макс. относит. расхождение = %.2e" % rel)

    # --- 3. формула (*) воспроизводит побинное нарушение
    a2, b2 = w.ab(freqs)
    dt2 = np.abs((a2 - b2) / 2.0) ** 2
    ana = imag_leak_analytic(a1, b1, dt2, 33.0, 0.0)
    num = np.imag(np.conj(zc) * zs)
    rel = float(np.max(np.abs(ana - num)) / max(np.max(np.abs(num)), 1e-300))
    chk("формула (*) воспроизводит Im(conj(z_c) z_s) побинно", rel < 1e-9,
        "макс. относит. расхождение = %.2e" % rel)

    # --- 4. геометрические нули: phi1 = psi_e и psi_e + 90
    for phi in (0.0, 90.0, 180.0):
        Yz = simulate(theta, freqs, phi)
        rz = rank_of(theta, Yz)["lambda2_over_lambda1"]
        chk("ранг 1 ВОССТАНАВЛИВАЕТСЯ при phi1 = %.0f (ось вдоль/поперёк E_in)" % phi,
            rz < 1e-20, "lambda2/lambda1 = %.2e" % rz)

    # --- 5. один элемент — контроль: нарушения нет
    Y1 = np.stack([_one_element(theta, a, b, 0.0) for a, b in zip(a1, b1)], axis=1)
    r1 = rank_of(theta, Y1)["lambda2_over_lambda1"]
    chk("один элемент: ранг 1 точный (контроль)", r1 < 1e-20,
        "lambda2/lambda1 = %.2e" % r1)

    # --- 6. ГДЕ максимум. Числитель ~ sin(2*delta) гасится к delta = 90, но
    #        знаменатель (Малюс неподвижного элемента) гаснет БЫСТРЕЕ, как
    #        |a|^2 eps^2 + |b|^2 при eps = 90 - delta. Максимум отношения — при
    #        eps = |b|/|a| = eta, то есть при phi1 ПОЧТИ СКРЕЩЕННОМ с E_in, а
    #        вовсе не на 45 град. Практический смысл прямой: нарушение ранга 1
    #        сильнее всего в рабочем режиме аттенюатора — при большом затухании.
    eta_el = float(np.median(np.abs(b1) / np.abs(a1)))
    phi_pred = 90.0 - np.degrees(eta_el)
    sc = device_scan(n_phi=901, verbose=False)
    imax = int(np.argmax(sc["lambda_ratio_numeric"]))
    phi_num = sc["phi_deg"][imax]
    # Максимум симметричен относительно 90: отклонение +-eta в обе стороны.
    chk("максимум при |phi1 - 90| = eta (почти скрещенный неподвижный элемент)",
        abs(abs(phi_num - 90.0) - (90.0 - phi_pred)) < 0.3,
        "численно |phi1 - 90| = %.2f град, формула eta = %.2f град"
        % (abs(phi_num - 90.0), 90.0 - phi_pred))

    # --- 6b. и величина максимума ~ (sin dphi / 2)^2
    dphi = np.angle(a1) - np.angle(b1)
    pred_max = float(np.max((np.sin(dphi) / 2.0) ** 2))
    got = float(sc["max_value"])
    chk("величина максимума согласуется с (sin dphi / 2)^2 по порядку",
        0.2 < got / max(pred_max, 1e-300) < 5.0,
        "численно %.3f, формула %.3f (отношение %.2f)"
        % (got, pred_max, got / max(pred_max, 1e-300)))

    # --- 7. масштаб: нарушение ~ eta_ЭЛЕМЕНТА^2, а не eta_системы^2
    Y45 = simulate(theta, freqs, 45.0)
    r45 = rank_of(theta, Y45)["lambda2_over_lambda1"]
    eta_el = float(np.median(np.abs(b1) / np.abs(a1)))
    chk("масштаб нарушения задаётся eta ЭЛЕМЕНТА",
        0.02 < np.sqrt(r45) / eta_el < 2.0,
        "sqrt(lambda2/lambda1) = %.4f при eta элемента = %.4f (отношение %.2f)"
        % (np.sqrt(r45), eta_el, np.sqrt(r45) / eta_el))

    # --- 8. измеренное на приборе ДОСТИЖИМО моделью двух элементов
    m = match_measured("att-11-16-356", verbose=False)
    chk("измеренное 2.11e-03 достижимо моделью двух элементов", m["reachable"],
        "максимум модели %.2e, совпадение при phi = %s град"
        % (m["model_max"], m["phi_matching_deg"]))

    n_ok = sum(c["ok"] for c in checks)
    if verbose:
        print("\n  Итог: %d/%d проверок пройдено" % (n_ok, len(checks)))
    return {"n_ok": n_ok, "n_total": len(checks), "checks": checks}


def _one_element(theta_deg, a, b, t0):
    t = np.deg2rad(np.asarray(theta_deg, dtype=float) - t0)
    return np.cos(t) ** 2 * a + np.sin(t) ** 2 * b


# ----------------------------------------------------------------------- CLI
def _jsonable(o):
    if isinstance(o, np.ndarray):
        return [_jsonable(x) for x in o.tolist()]
    if isinstance(o, (np.floating, np.integer)):
        return o.item()
    if isinstance(o, dict):
        return {k: _jsonable(v) for k, v in o.items()}
    if isinstance(o, (list, tuple)):
        return [_jsonable(x) for x in o]
    return o


def main():
    ap = argparse.ArgumentParser(
        description="A24 — провал ранга 1 на приборе: два элемента")
    ap.add_argument("--selftest", action="store_true")
    ap.add_argument("--device", action="store_true", help="скан по углу неподвижного элемента")
    ap.add_argument("--data", action="store_true", help="проверка роста с частотой на данных")
    a = ap.parse_args()

    res = {"task": "A24", "passport": PASSPORT, "measured": MEASURED}
    if a.selftest:
        print("A24 selftest — алгебра двух элементов")
        res["selftest"] = selftest()
    if a.device:
        print("\n--- модель прибора ---")
        res["device_scan"] = device_scan()
        res["match"] = [match_measured(d) for d in ("att-11-16-356", "att-11-16-s1")]
        print()
        res["joint"] = [joint_match(d) for d in ("att-11-16-356", "att-11-16-s1")]
    if a.data:
        print("\n--- реальные данные: растёт ли нарушение с частотой ---")
        res["data"] = [data_frequency_test(d) for d in
                       ("att-11-16-356", "att-11-16-s1", "test_grid_33_11", "specac-2")]
    if not (a.selftest or a.device or a.data):
        ap.error("нужен --selftest, --device или --data")

    (OUT / "a24_two_element_rank.json").write_text(
        json.dumps(_jsonable(res), ensure_ascii=False, indent=2), encoding="utf-8")
    st = res.get("selftest")
    return 0 if (st is None or st["n_ok"] == st["n_total"]) else 1


if __name__ == "__main__":
    raise SystemExit(main())
