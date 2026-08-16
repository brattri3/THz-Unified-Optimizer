#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
A22_field_linearization — линеаризация закона Малюса в КОМПЛЕКСНОМ ПОЛЕ.

Идея одной фразой
-----------------
Аппарат A0 (`malus_linearization.py`) линеаризует МОЩНОСТЬ U = |E|^2 и потому
работает в базисе 4-го порядка. Но данные проекта комплексные до самого конца
(`exp_builder.py:94`: exp = FFT(sig)/FFT(bg)), а глобальный фит фазу уже
использует (`model_core.py:329-337`). Если линеаризовать САМО ПОЛЕ, базис
становится 2-го порядка, а на выходе получаются КОМПЛЕКСНЫЕ спектры каналов.

    E(theta) = T_perp*cos^2(theta-theta0) + T_par*sin^2(theta-theta0)
             = z0 + z_c*cos(2theta) + z_s*sin(2theta)

Три КОМПЛЕКСНЫХ коэффициента, базис {1, cos2t, sin2t} — вещественный и от
данных не зависящий. Отсюда:

    tan(2*theta0) = z_s/z_c    <- обязано быть ВЕЩЕСТВЕННЫМ и одним на все бины
    B = z_c*cos(2theta0) + z_s*sin(2theta0)
    T_perp = z0 + B,   T_par = z0 - B          (комплексные, с фазами)

Что это даёт сверх мощностной версии
------------------------------------
1. arg(T_perp) и arg(T_par) ПО ОТДЕЛЬНОСТИ, а не только cos(dphi). Ветвь знака
   `arccos` (открытый вопрос NOTE_malus_linearization.md §8) исчезает: фаза
   берётся как аргумент комплексного числа.
2. tau_leak — наклоном развёрнутой фазы, БЕЗ закрепления D_eff (в M5 он
   идентифицируем только при пине, HN8).
3. |cos dphi| <= 1 и |T_par|^2 >= 0 выполняются ТОЖДЕСТВЕННО, при любом шуме.
   В мощности они нарушаются при низком SNR (измерено: SNR 3.6 -> 19.8 %
   нарушений), что и объясняет 12/151 у 356att и 43/94 у series1.
4. Четвёртая гармоника в ПОЛЕ обязана быть нулём. Любой значимый (z_c4, z_s4) —
   выход за схему «один вращающийся WGP x скаляр». Мощностная формулировка
   такого теста сделать не может: там A4/A2 = 0.25 есть ожидаемая физика.

Чего этот метод НЕ даёт
-----------------------
Вырождение D_eff <-> eta он НЕ снимает: оно структурное, на уровне отображения
параметров модели в пару (T_perp, T_par), а корреляции -0.9999993 получены
фитом, который фазу уже использовал. Это зафиксировано отрицательным контролем
в selftest.

Веса: модель шума калибруется, а не назначается
-----------------------------------------------
sigma^2(theta,nu) = kappa*tnoise(nu)*(1+|t|^2) + (eps_mult*|t|)^2 — аддитивный
шум детектора плюс мультипликативный дрейф уровня трассы.

Общий множитель kappa*tnoise(nu) от УГЛА не зависит и в решении ВМНК
сокращается, поэтому САМ ПО СЕБЕ шумовой пол точечные оценки не двигает. Но
угловую зависимость дают оба остальных члена, и она существенна: вес тени
относительно яркой зоны у чисто аддитивной и у смешанной модели различается на
два порядка. Измерено на test_grid_33_11: чисто аддитивные веса дают
theta0 = -1.3837, веса ~|t|^2 (рецепт A) -1.2700, а откалиброванная смесь
-1.2723 +- 0.0063 против эталона A13 -1.2770 +- 0.0354.

Поэтому eps_mult подбирается из chi^2/dof = 1 (`calibrate_noise`), а не
берётся с потолка. Проверка независимая: на specac/exp_data_2 есть 2-4 повтора
на угол, и разброс МЕЖДУ ПОВТОРАМИ сравнивается с моделью напрямую — чисто
аддитивная занижена в 42 раза, откалиброванная согласуется (0.58).

Это уточняет, а не отменяет вердикт A по вопросу C-5: там речь шла о
мощностном фите, где sigma по шумовому полу двигает саму оценку eta. В поле
опасность иная — не сам шумовой пол, а ПОЛНОТА модели шума.

Обозначения
-----------
theta   угол поворота WGP, град (конвенция репозитория: ось ПРОПУСКАНИЯ)
T_perp  комплексное пропускание «яркого» канала (E перпендикулярно проводам)
T_par   то же для «тёмного» канала (утечка)
eta     амплитудная утечка |T_par|/|T_perp|
dphi    относительная фаза arg(T_perp) - arg(T_par) = delta0 + 2*pi*nu*tau_leak
theta0  угловой офсет (нуль шкалы поворотного столика)
r       отношение T_par/T_perp; инвариантно к общему множителю (нормировка на
        фон, скаляр t_perp2(nu) второго WGP по теореме выравнивания T8)

Зона A. Ядро зоны B (fit_lib/model_core) не правится; `malus_linearization.py`
не правится тоже — у него порт-копия в track_viewer (вопрос C-1).

Запуск:  .venv\\Scripts\\python.exe research/two_wgp/field_linearization.py --selftest
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

#: E|n| = sigma*sqrt(pi/4) для кругового комплексного гауссова шума, а
#: `noise_floor` берёт nfa = среднее МОДУЛЯ FFT (a3_eps_stability.py:133).
#: Значит tnoise занижен по ДИСПЕРСИИ ровно в 4/pi раз.
KAPPA_ABS_TO_VAR = 4.0 / np.pi


# --------------------------------------------------------------- базис
def field_design(theta_deg, order=2):
    """Матрица плана X: столбцы [1, cos2t, sin2t] (+ [cos4t, sin4t] при order=4).

    Поле ЛИНЕЙНО по коэффициентам при order=2 — это точная структура схемы
    «один вращающийся WGP». Столбцы 4-й гармоники добавляются только как
    ТЕСТОВЫЕ: в поле их коэффициенты обязаны быть нулём.
    """
    t = np.asarray(theta_deg, dtype=np.float64) * D2R
    cols = [np.ones_like(t), np.cos(2 * t), np.sin(2 * t)]
    if order >= 4:
        cols += [np.cos(4 * t), np.sin(4 * t)]
    return np.column_stack(cols)


def bound_design(theta_deg, theta0_deg):
    """Связанный базис [1, cos2(t-t0)] — 2 столбца при ИЗВЕСТНОМ офсете.

    Даёт T_perp/T_par с меньшей дисперсией и честным dof = 2n-4.
    """
    t = (np.asarray(theta_deg, dtype=np.float64) - float(theta0_deg)) * D2R
    return np.column_stack([np.ones_like(t), np.cos(2 * t)])


# --------------------------------------------------- веса и модель шума
def noise_variance(exp, tnoise, kappa=KAPPA_ABS_TO_VAR, amp=None, eps_mult=0.0):
    """sigma^2[theta,nu] = kappa*tnoise[nu]*(1+|t|^2) + (eps_mult*|t|)^2.

    ДВЕ компоненты, и обе измерены на данных:

    * АДДИТИВНАЯ — шум детектора. Шумит и сигнал, и фон: t = E_s/E_b, откуда
      в первом порядке var(t) = (var_s + |t|^2 var_b)/|E_b|^2.
    * МУЛЬТИПЛИКАТИВНАЯ — дрейф уровня трассы (мощность лазера,
      перепозиционирование). Одной аддитивной моделью chi^2/dof выходит 8.7
      на test_grid_33_11, а снятие комплексного множителя на трассу роняет
      его до 2.6 при разбросе |g| ~ 1 % — то есть шум существенно
      мультипликативный (это ровно HN6).

    Практическое следствие крупное: вес тени относительно яркой зоны
    отличается у двух моделей на ДВА порядка, и от него заметно зависит
    theta0 (измерено: -1.3837 при чисто аддитивных весах против -1.2700 при
    весах ~|t|^2, то есть рецепт A). Поэтому eps_mult не назначается, а
    калибруется по данным — `calibrate_noise`.

    `amp` — модуль, по которому берётся |t|. По умолчанию из самих данных;
    правильнее подавать ПОДОГНАННЫЙ |E| (одна итерация IRLS), иначе вес точки
    коррелирует с её собственным шумом.
    """
    exp = np.asarray(exp)
    tn = np.asarray(tnoise, dtype=np.float64)
    a = np.abs(exp) if amp is None else np.abs(np.asarray(amp))
    return float(kappa) * tn[None, :] * (1.0 + a ** 2) + (float(eps_mult) * a) ** 2


def calibrate_noise(theta_deg, Y, tnoise, kappa=KAPPA_ABS_TO_VAR,
                    theta0_deg=None, lo=0.0, hi=0.2, n_iter=40):
    """Подбирает eps_mult из условия chi^2/dof = 1. -> (eps_mult, chi2_red).

    Модель шума калибруется ПО ДАННЫМ, а не назначается: половинным делением
    по монотонной функции chi^2_red(eps). Если даже при eps = hi невязка выше
    единицы, возвращается hi и фактический chi^2_red — это сигнал, что часть
    невязки структурная, а не шумовая, и её надо смотреть в RSS_perp.
    """
    def _chi(eps):
        s2 = noise_variance(Y, tnoise, kappa=kappa, eps_mult=eps)
        t0 = theta0_deg
        if t0 is None:
            t0 = theta0_rank1(solve_bins(theta_deg, Y, s2, order=2).z)["theta0_deg"]
        fb = solve_bins(theta_deg, Y, s2, theta0_deg=t0)
        return float(np.mean(fb.chi2) / fb.dof)

    c_lo, c_hi = _chi(lo), _chi(hi)
    if c_hi > 1.0:
        return float(hi), c_hi
    if c_lo < 1.0:
        return float(lo), c_lo
    a, b = lo, hi
    for _ in range(n_iter):
        m = 0.5 * (a + b)
        if _chi(m) > 1.0:
            a = m
        else:
            b = m
    m = 0.5 * (a + b)
    return float(m), _chi(m)


# --------------------------------------------------- побинный комплексный ВМНК
class BinFit(object):
    """Результат побинного полевого ВМНК."""

    __slots__ = ("z", "C", "chi2", "dof", "resid", "N", "order", "theta0_deg",
                 "n_ang", "n_bin")

    def __init__(self):
        self.z = None            # (k, n_bin) комплексные коэффициенты
        self.C = None            # (n_bin, k, k) вещественная = inv(X^T W X)
        self.chi2 = None         # (n_bin,)
        self.dof = 0             # 2*n_ang - 2*k  (комплексная невязка)
        self.resid = None        # (n_ang, n_bin) комплексные
        self.N = None            # (n_bin, k, k) информационная матрица
        self.order = 2
        self.theta0_deg = None   # не None => связанный фит


def solve_bins(theta_deg, Y, sigma2=None, theta0_deg=None, order=2):
    """Комплексный взвешенный МНК, независимо по каждому частотному бину.

    Y — (n_ang, n_bin) комплексная матрица t(theta,nu) = FFT(sig)/FFT(bg).
    sigma2 — (n_ang, n_bin) либо None (равные веса).

    Матрица плана ВЕЩЕСТВЕННАЯ, поэтому Re и Im решаются одной и той же
    системой; «комплексный МНК» в смысле эрмитовых наименьших квадратов не
    нужен. Для кругового шума Cov(z) = inv(X^T W X) — тоже вещественная, а
    псевдоковариация равна нулю (круговость сохраняется).
    """
    theta_deg = np.asarray(theta_deg, dtype=np.float64)
    Y = np.asarray(Y)
    if Y.ndim == 1:
        Y = Y[:, None]
    n_ang, n_bin = Y.shape

    if theta0_deg is None:
        X = field_design(theta_deg, order)
    else:
        X = bound_design(theta_deg, theta0_deg)
    k = X.shape[1]
    if n_ang < k:
        raise ValueError("углов %d < параметров %d: система недоопределена"
                         % (n_ang, k))

    W = (np.ones_like(Y, dtype=np.float64) if sigma2 is None
         else 1.0 / np.asarray(sigma2, dtype=np.float64))

    # N[b] = sum_i w[i,b] * X[i,:] X[i,:]^T ;  rhs[b] = sum_i w[i,b] X[i,:] Y[i,b]
    N = np.einsum("ib,ik,il->bkl", W, X, X)
    rhs = np.einsum("ib,ik,ib->bk", W, X, Y)
    z = np.linalg.solve(N, rhs[..., None])[..., 0].T          # (k, n_bin)

    resid = Y - X @ z
    # Множитель 2: sigma^2 — ПОЛНАЯ дисперсия комплексной величины (E|n|^2),
    # поэтому Re и Im имеют по sigma^2/2. Без него E[chi^2] вышло бы вдвое
    # меньше числа ВЕЩЕСТВЕННЫХ степеней свободы (измерено в selftest: 0.4984).
    chi2 = 2.0 * np.einsum("ib,ib->b", W, np.abs(resid) ** 2)

    res = BinFit()
    res.z, res.N = z, N
    res.C = np.linalg.inv(N)
    res.chi2 = chi2
    res.dof = 2 * n_ang - 2 * k
    res.resid = resid
    res.order = int(order) if theta0_deg is None else 2
    res.theta0_deg = theta0_deg
    res.n_ang, res.n_bin = n_ang, n_bin
    return res


# ------------------------------------------------------- угловой офсет
def theta0_rank1(z, w=None):
    """theta0 в ЗАМКНУТОЙ форме + мера нарушения гипотезы единого офсета.

    Вектор (z_c, z_s) равен B_j*u, где u = (cos2t0, sin2t0) — ВЕЩЕСТВЕННЫЙ
    единичный вектор, общий для всех бинов. Значит

        R_j = [[|z_c|^2, Re(conj(z_c) z_s)], [.., |z_s|^2]] = |B_j|^2 * u u^T

    точно ранга 1, и M = sum_j w_j R_j — тоже. u — старший собственный вектор.

    lambda2/lambda1 — прямая скалярная мера нарушения единого theta0, то есть
    модельно-независимый тест HN13. Побинный Im(conj(z_c) z_s) — та же величина
    с разрешением по частоте.
    """
    z = np.asarray(z)
    zc, zs = z[1], z[2]
    if w is None:
        w = np.ones(zc.shape, dtype=np.float64)
    w = np.asarray(w, dtype=np.float64)

    cross = np.conj(zc) * zs
    m11 = float(np.sum(w * np.abs(zc) ** 2))
    m22 = float(np.sum(w * np.abs(zs) ** 2))
    m12 = float(np.sum(w * cross.real))
    M = np.array([[m11, m12], [m12, m22]], dtype=np.float64)

    lam, vec = np.linalg.eigh(M)                     # по возрастанию
    u = vec[:, -1]
    lam1, lam2 = float(lam[-1]), float(lam[0])

    # Ветвь: u определён с точностью до ЗНАКА, а смена знака меняет каналы
    # местами и сдвигает theta0 на 90 град. Фиксируем требованием, что ЯРКИЙ
    # канал ярче: sum|T_perp|^2 > sum|T_par|^2 по всей полосе.
    B = z[1] * u[0] + z[2] * u[1]
    if float(np.sum(np.abs(z[0] + B) ** 2)) < float(np.sum(np.abs(z[0] - B) ** 2)):
        u = -u
    theta0 = 0.5 * np.arctan2(u[1], u[0]) * R2D

    denom = np.maximum(np.abs(zc) * np.abs(zs), 1e-300)
    return {
        "theta0_deg": float(theta0),
        "lambda2_over_lambda1": float(lam2 / lam1) if lam1 > 0 else np.nan,
        "imag_leak": cross.imag / denom,             # (n_bin,) безразмерная
        "u": u,
    }


def fix_branch(theta0_deg, theta_deg, Y, sigma2=None):
    """Приводит theta0 к ветви, в которой ЯРКИЙ канал действительно ярче.

    chi^2(theta0) имеет период 90 град: сдвиг на 90 просто меняет T_perp и
    T_par местами, невязка та же. Поэтому профиль сам ветвь не выбирает.
    """
    fb = solve_bins(theta_deg, Y, sigma2, theta0_deg=theta0_deg)
    Tp, Tl, _ = channels_from_z(fb)
    if float(np.sum(np.abs(Tp) ** 2)) < float(np.sum(np.abs(Tl) ** 2)):
        theta0_deg = theta0_deg + 90.0
    return float((theta0_deg + 90.0) % 180.0 - 90.0)


def theta0_profile(theta_deg, Y, sigma2=None, span_deg=(-45.0, 45.0), n=901,
                   birge=True):
    """ML-профиль chi^2(theta0): одномерный скан ЦЕЛИКОМ + уточнение параболой.

    Сканируется вся область, а не ищется локальный минимум: мультимодальность
    видна глазом, и мультистарт не нужен (урок A2/A3 — угловые добавочные
    параметры мультимодальны).

    sigma(theta0) — из Delta chi^2 = 1. При `birge=True` профиль сначала
    делится на chi^2/dof: это поправка Бирджа, честная когда известен не
    абсолютный уровень шума, а лишь его ФОРМА. Здесь именно такой случай —
    tnoise видит только ВЧ-полку спектра и не знает про шум УРОВНЯ УГЛА
    (дрейф лазера, перепозиционирование), из-за чего chi^2/dof выходит
    заметно больше единицы. Без поправки sigma(theta0) занижалась бы в
    sqrt(chi^2/dof) раз. Тот же приём — в `malus_linearization.theta0_error`.
    """
    grid = np.linspace(span_deg[0], span_deg[1], int(n))
    chi = np.empty_like(grid)
    for i, t0 in enumerate(grid):
        chi[i] = float(np.sum(solve_bins(theta_deg, Y, sigma2, theta0_deg=t0).chi2))

    i0 = int(np.argmin(chi))
    best = float(grid[i0])
    h = float(grid[1] - grid[0])
    den = np.nan
    if 0 < i0 < len(grid) - 1:                        # параболическое уточнение
        y0, y1, y2 = chi[i0 - 1], chi[i0], chi[i0 + 1]
        den = float(y0 - 2 * y1 + y2)
        if den > 0:
            best = float(grid[i0] + 0.5 * (y0 - y2) / den * h)

    chi_min = float(chi[i0])
    nb = Y.shape[1] if np.asarray(Y).ndim > 1 else 1
    dof_tot = max((2 * len(theta_deg) - 4) * nb, 1)
    scale = max(chi_min / dof_tot, 1e-300) if birge else 1.0

    # sigma — из КРИВИЗНЫ минимума, а не перебором узлов: минимум обычно острее
    # шага сетки (здесь sigma ~ 0.03 град при шаге 0.1), и перебор дал бы nan.
    sigma_t0 = (float(h * np.sqrt(2.0 * scale / den)) if np.isfinite(den) and den > 0
                else float("nan"))
    best = fix_branch(best, theta_deg, Y, sigma2)
    return {"theta0_deg": best, "sigma_theta0_deg": sigma_t0,
            "chi2_min": chi_min, "birge_scale": float(scale),
            "grid_deg": grid, "chi2_curve": chi}


# ------------------------------------------------------- каналы и отношение
def channels_from_z(fit):
    """(T_perp, T_par, Cov2) из СВЯЗАННОГО фита (2 столбца, theta0 известен).

    T_perp = z0 + zB,  T_par = z0 - zB;  M = [[1,1],[1,-1]];  Cov2 = M C M^T.
    Cov2 вещественная 2x2 на бин: Cov(dT_a, conj(dT_b)).
    """
    if fit.theta0_deg is None:
        raise ValueError("нужен связанный фит (theta0_deg задан)")
    z0, zB = fit.z[0], fit.z[1]
    M = np.array([[1.0, 1.0], [1.0, -1.0]])
    Cov2 = np.einsum("ak,bkl,cl->bac", M, fit.C, M)          # (n_bin, 2, 2)
    return z0 + zB, z0 - zB, Cov2


def ratio_and_var(T_perp, T_par, Cov2):
    """r = T_par/T_perp и дисперсия ln r. -> (r, var_lnr).

    d(ln r) = dT_par/T_par - dT_perp/T_perp, откуда

        E|d ln r|^2 = S_par/|T_par|^2 + S_perp/|T_perp|^2
                      - 2*Re( S_cross/(T_par*conj(T_perp)) )

    Круговость сохраняется, поэтому Var(ln|r|) = Var(arg r) = E|d ln r|^2 / 2 —
    амплитудная и фазовая части имеют ОДНУ И ТУ ЖЕ дисперсию.
    """
    Sp = Cov2[:, 0, 0]                # perp
    Sl = Cov2[:, 1, 1]                # par
    Sx = Cov2[:, 0, 1]                # cross (вещественная)
    r = T_par / T_perp
    with np.errstate(divide="ignore", invalid="ignore"):
        term = (Sl / np.abs(T_par) ** 2 + Sp / np.abs(T_perp) ** 2
                - 2.0 * np.real(Sx / (T_par * np.conj(T_perp))))
    return r, 0.5 * term


# ------------------------------------------------------- развёртка фазы
def unwrap_segmented(freqs, phase, valid, sigma=None):
    """np.unwrap ПО КУСКАМ непрерывных валидных бинов + сшивка целым 2*pi*k.

    Водяная маска рвёт частотную ось, а через разрыв развёртка неоднозначна.
    Куски сшиваются целым k, минимизирующим отклонение от линейного тренда
    предыдущего куска; в отчёт идёт Delta chi^2 между лучшим и вторым k —
    при Delta chi^2 < 4 сшивка помечается НЕОДНОЗНАЧНОЙ.
    """
    freqs = np.asarray(freqs, dtype=np.float64)
    phase = np.asarray(phase, dtype=np.float64)
    valid = np.asarray(valid, dtype=bool)
    out = np.full_like(phase, np.nan)

    idx = np.where(valid)[0]
    if not len(idx):
        return out, {"n_segments": 0, "n_joins": 0, "min_delta_chi2": np.nan,
                     "ambiguous": False}

    cuts = np.where(np.diff(idx) > 1)[0]
    segs = np.split(idx, cuts + 1)

    joins, dchi = 0, []
    prev_f = prev_p = None
    for seg in segs:
        p = np.unwrap(phase[seg])
        if prev_f is not None and len(prev_f) >= 2:
            a, b = np.polyfit(prev_f, prev_p, 1)          # тренд предыдущего куска
            pred = a * freqs[seg[0]] + b
            ks = np.arange(-6, 7)
            err = np.abs(p[0] + 2 * np.pi * ks - pred)
            order = np.argsort(err)
            p = p + 2 * np.pi * ks[order[0]]
            joins += 1
            dchi.append(float(err[order[1]] ** 2 - err[order[0]] ** 2))
        out[seg] = p
        prev_f, prev_p = freqs[seg], p

    md = float(np.min(dchi)) if dchi else np.inf
    return out, {"n_segments": len(segs), "n_joins": joins,
                 "min_delta_chi2": md, "ambiguous": bool(md < 4.0)}


def logratio_regression(freqs, r, var_lnr, valid=None):
    """Обе части ln r сразу: Re -> степенной закон утечки, Im -> задержка.

        ln|r| = ln(eta0) + a*ln(nu)                   (HN2, модельно-независимо)
        dphi(nu) = arg(T_perp) - arg(T_par) = -arg r  (развёртывается по частоте)

    Измеряется НАКЛОН dphi по частоте — величина `tau_dphi_ps`, групповая
    задержка МЕЖДУ КАНАЛАМИ. Это то, что даёт эксперимент; связь с параметром
    модели требует модели, и её надо выписывать явно. В конвенции M5
    (`model_core.py:275,293-294`) утечка входит как exp(+i*2*pi*nu*tau_leak), а
    общая задержка канала — как exp(-i*2*pi*nu*tau), откуда для режима, где
    утечка доминирует над блансовским t_par:

        dphi(nu) ~ const + 2*pi*nu*(tau_par - tau_leak)
        =>  tau_dphi = tau_par - tau_leak,  и при tau_par ~ 0  tau_leak = -tau_dphi

    Поэтому рядом кладётся `tau_leak_m5_ps = -tau_dphi_ps` — пересчёт в
    конвенцию проекта, а не переименование: на att-11-16-356 он даёт +0.187 пс
    против 0.195 пс у A0 и 0.19 пс у M5. Знак здесь не косметика: у A0 фаза
    восстанавливалась через arccos и знаковой ветви не имела вовсе (открытый
    вопрос NOTE §8), а в поле она измеряется прямо.

    Возвращает eta0, eta_exp, delta0, задержки с ошибками, chi2_red и остатки —
    вывод о качестве делается по статистике остатков, а не по одному числу.
    """
    freqs = np.asarray(freqs, dtype=np.float64)
    r = np.asarray(r)
    var = np.asarray(var_lnr, dtype=np.float64)
    if valid is None:
        valid = np.isfinite(var) & (var > 0) & np.isfinite(np.abs(r)) & (np.abs(r) > 0)
    valid = np.asarray(valid, dtype=bool)

    out = {"n_used": int(valid.sum())}
    if valid.sum() < 3:
        out.update({k: np.nan for k in ("eta0", "eta_exp", "delta0_rad",
                                        "tau_dphi_ps", "sigma_tau_dphi_ps",
                                        "tau_leak_m5_ps")})
        return out

    ph, uinfo = unwrap_segmented(freqs, np.angle(r), valid)
    out["unwrap"] = uinfo

    f = freqs[valid]
    w = 1.0 / var[valid]

    def _wls(x, y, w):
        X = np.column_stack([np.ones_like(x), x])
        N = X.T @ (w[:, None] * X)
        p = np.linalg.solve(N, X.T @ (w * y))
        C = np.linalg.inv(N)
        rss = float(np.sum(w * (y - X @ p) ** 2))
        return p, C, rss, y - X @ p

    # --- амплитуда: степенной закон
    pa, Ca, rss_a, res_a = _wls(np.log(f), np.log(np.abs(r[valid])), w)
    out["eta0"] = float(np.exp(pa[0]))
    out["eta_exp"] = float(pa[1])
    out["sigma_eta_exp"] = float(np.sqrt(max(Ca[1, 1], 0.0)))
    out["chi2_red_amp"] = rss_a / max(len(f) - 2, 1)
    # Полная 2x2 ковариация в координатах (ln eta0, a) — она нужна для тестов
    # однородности между датасетами: ln eta0 и a скоррелированы сильно
    # (общий рычаг ln nu), и по одним диагоналям тест построить нельзя.
    out["cov_amp"] = Ca
    out["p_amp"] = pa

    # --- фаза: dphi = -arg r; наклон = 2*pi*tau_dphi (см. докстринг)
    pp, Cp, rss_p, res_p = _wls(f, -ph[valid], w)
    out["delta0_rad"] = float(pp[0])
    out["tau_dphi_ps"] = float(pp[1] / (2 * np.pi))
    out["sigma_tau_dphi_ps"] = float(np.sqrt(max(Cp[1, 1], 0.0)) / (2 * np.pi))
    out["tau_leak_m5_ps"] = -out["tau_dphi_ps"]        # конвенция model_core
    out["chi2_red_phase"] = rss_p / max(len(f) - 2, 1)
    out["cov_phase"] = Cp
    out["resid_amp"], out["resid_phase"] = res_a, res_p
    out["weights_used"] = w
    out["freqs_used"] = f
    return out


# ------------------------------------------------------- структурные тесты
def order4_field_test(theta_deg, Y, sigma2=None):
    """Четвёртая гармоника В ПОЛЕ обязана быть нулём.

    В мощности A4/A2 = 0.25 — ожидаемая физика когерентной композиции, и тест
    невозможен. В поле любой значимый (z_c4, z_s4) означает выход за схему
    «один вращающийся WGP x скаляр»: Фабри-Перо, нелинейность детектора,
    угол-зависимая фаза (джиттер).

    -> chi2 разности (4 вещественные степени свободы на бин) и p-уровень.
    """
    f2 = solve_bins(theta_deg, Y, sigma2, order=2)
    f4 = solve_bins(theta_deg, Y, sigma2, order=4)
    d = f2.chi2 - f4.chi2                       # >= 0
    n_bin = len(d)
    total = float(np.sum(d))
    dof = 4 * n_bin
    try:
        from scipy import stats
        p = float(stats.chi2.sf(total, dof)) if dof > 0 else np.nan
    except Exception:                            # noqa: BLE001
        p = np.nan
    return {"delta_chi2_per_bin": d, "delta_chi2_total": total,
            "dof": dof, "p_value": p,
            "mean_per_dof": total / dof if dof else np.nan}


def out_of_span(fit):
    """RSS_perp(nu) — часть невязки, недостижимая НИКАКОЙ моделью вида
    cos^2*T_perp + sin^2*T_par.

    При фиксированном theta0 модель целиком лежит в span(X), поэтому невязка
    ортогональной части не зависит от параметров вообще. Это разделяет
    «модель могла бы объяснить» и «не могла бы никогда».
    """
    return {"rss_perp_per_bin": fit.chi2, "rss_perp_total": float(np.sum(fit.chi2)),
            "dof": fit.dof, "chi2_red": float(np.mean(fit.chi2) / fit.dof)
            if fit.dof > 0 else np.nan}


def extinction_geometry(T_perp, T_par):
    """Где на самом деле минимум |E(theta)| — и почему не обязательно на 90 град.

        x* = -(|T_perp|^2 - |T_par|^2) / |T_perp - T_par|^2

    При |x*| < 1 минимум ВНУТРЕННИЙ: theta_min = theta0 + 0.5*arccos(x*), то
    есть угол наилучшего гашения отличается от theta0+90 БЕЗ всякой асимметрии
    образца. Критерий внутренности сводится к cos(dphi) < eta.

    Кандидат-объяснение «необъяснённых ~3 град» у test_grid_33_11 и части
    асимметрии HN13 — без единого нового параметра.
    """
    T_perp = np.asarray(T_perp)
    T_par = np.asarray(T_par)
    num = np.abs(T_perp) ** 2 - np.abs(T_par) ** 2
    den = np.abs(T_perp - T_par) ** 2
    with np.errstate(divide="ignore", invalid="ignore"):
        xs = -num / den
    inner = np.abs(xs) < 1.0
    dtheta = np.where(inner, 0.5 * np.arccos(np.clip(xs, -1, 1)) * R2D, 90.0)
    return {"x_star": xs, "inner": inner, "dtheta_min_deg": dtheta}


# ------------------------------------------------------- мощностная сверка
def power_bias_prediction(T_perp, T_par, tnoise, kappa=KAPPA_ABS_TO_VAR):
    """Что ДОЛЖЕН был показать мощностный аппарат A0 при данном шумовом поле.

    E|E+n|^2 = |E|^2 + sigma^2: шум добавляется к оценке мощности аддитивно и
    одинаково на всех углах, то есть садится в c0. Для НЕвзвешенного фита это
    даёт завышение eta; при рецепте A (sigma ~ U) знак переворачивается, и
    смещение мало. Функция возвращает оценку «наивного» варианта — она нужна
    как опорная точка, а не как предсказание рабочего режима.
    """
    s2 = float(kappa) * np.asarray(tnoise, dtype=np.float64)
    ip = np.abs(T_perp) ** 2 + s2
    il = np.abs(T_par) ** 2 + s2
    return np.sqrt(np.maximum(il, 0.0) / np.maximum(ip, 1e-300))


# =========================================================== самопроверка
def _synth_field(theta_deg, T_perp, T_par, theta0_deg):
    """Прямая задача в поле: E(theta) = T_perp cos^2 + T_par sin^2."""
    t = (np.asarray(theta_deg, dtype=np.float64) - theta0_deg) * D2R
    return T_perp * np.cos(t) ** 2 + T_par * np.sin(t) ** 2


def _power_channels(theta_deg, U):
    """Мощностная ветка (аппарат A0) — для сверки бок о бок.

    Импортируем `malus_linearization`, а не копируем: копия потребовала бы
    собственной машинной сверки (урок C-1).
    """
    import sys
    sys.path.insert(0, str(Path(__file__).resolve().parent))
    import malus_linearization as ml               # noqa: E402

    sigma = ml._weights_for(U, "relative")
    p, cov, info = ml.fit_harmonics(theta_deg, U, sigma=sigma, order=4)
    return ml.channel_params(p, cov)


def selftest(verbose=True):
    checks = []

    def chk(name, ok, detail):
        checks.append({"check": name, "ok": bool(ok), "detail": detail})
        if verbose:
            print(f"  [{'OK ' if ok else 'FAIL'}] {name}: {detail}")

    rng = np.random.default_rng(20260816)
    th = np.arange(-90.0, 91.0, 10.0)                 # 19 углов, как 356att
    T_perp, T_par, t0 = 0.93 + 0.11j, 0.031 * np.exp(1j * 0.7), 1.35
    eta_true = abs(T_par) / abs(T_perp)
    dphi_true = np.angle(T_perp) - np.angle(T_par)

    # --- 1. базис порядка 2 описывает поле ТОЧНО
    E = _synth_field(th, T_perp, T_par, t0)
    f2 = solve_bins(th, E[:, None], order=2)
    rms = float(np.sqrt(np.mean(np.abs(f2.resid) ** 2)))
    chk("полевой базис order=2 описывает схему точно", rms < 1e-15,
        f"rms невязки = {rms:.2e}")

    # --- 2. восстановление комплексных каналов и офсета
    r1 = theta0_rank1(f2.z)
    chk("theta0 замкнутой формулой (rank-1)", abs(r1["theta0_deg"] - t0) < 1e-10,
        f"{r1['theta0_deg']:.12f} против {t0}")
    chk("матрица (z_c,z_s) точно ранга 1", r1["lambda2_over_lambda1"] < 1e-25,
        f"lambda2/lambda1 = {r1['lambda2_over_lambda1']:.2e}")
    chk("Im(conj(z_c)*z_s) = 0 (единый theta0)",
        abs(float(r1["imag_leak"][0])) < 1e-12,
        f"Im-leak = {float(r1['imag_leak'][0]):.2e}")

    fb = solve_bins(th, E[:, None], theta0_deg=r1["theta0_deg"])
    Tp, Tl, Cov2 = channels_from_z(fb)
    chk("восстановление T_perp (комплексное)", abs(Tp[0] - T_perp) < 1e-12,
        f"|dT_perp| = {abs(Tp[0] - T_perp):.2e}")
    chk("восстановление T_par (комплексное)", abs(Tl[0] - T_par) < 1e-12,
        f"|dT_par| = {abs(Tl[0] - T_par):.2e}")
    dphi_rec = float(np.angle(Tp[0]) - np.angle(Tl[0]))
    chk("относительная фаза БЕЗ ветви arccos",
        abs(dphi_rec - dphi_true) < 1e-10,
        f"dphi = {dphi_rec:.10f} против {dphi_true:.10f} рад")

    # --- 3. стационарность по theta0 (эффект ВТОРОГО порядка).
    #        Свойство ОРТОГОНАЛЬНОГО плана, поэтому берём строго периодичный
    #        набор: у -90..+90 шагом 10 углы -90 и +90 физически совпадают
    #        (период 180), базис неортогонален и появляется линейная примесь.
    th_per = np.arange(-90.0, 90.0, 10.0)             # 18 углов, без дубля
    E_per = _synth_field(th_per, T_perp, T_par, t0)

    def _err(t0_shift, angles, data):
        Tp_, _, _ = channels_from_z(solve_bins(angles, data[:, None],
                                               theta0_deg=t0 + t0_shift))
        return abs(Tp_[0] - T_perp)

    e1, e2 = _err(0.2, th_per, E_per), _err(0.4, th_per, E_per)
    ratio = e2 / e1 if e1 > 0 else np.nan
    chk("T_perp стационарна по theta0 (ошибка ~ delta^2)",
        abs(ratio - 4.0) < 0.01, f"e(2d)/e(d) = {ratio:.4f} (ждём 4.00)")

    # то же на наборе с дублем -90/+90 — не проверка, а ИЗМЕРЕНИЕ примеси
    gram = field_design(th).T @ field_design(th) / len(th)
    chk("измерено: дубль -90/+90 делает базис неортогональным",
        True, f"Gram[1,const] = {gram[0, 1]:+.4f} (у строго периодичного 0), "
              f"e(2d)/e(d) = {_err(0.4, th, E) / max(_err(0.2, th, E), 1e-300):.2f}")

    # --- 4. четвёртая гармоника в поле = ноль
    o4 = order4_field_test(th, E[:, None])
    chk("4-я гармоника в ПОЛЕ равна нулю", o4["delta_chi2_total"] < 1e-20,
        f"delta chi^2 = {o4['delta_chi2_total']:.2e}")

    # --- 5. инвариантность к eps_cross (угол-независимая добавка)
    eps = 0.004 - 0.002j
    fe = solve_bins(th, (E + eps)[:, None], order=2)
    re = theta0_rank1(fe.z)
    fbe = solve_bins(th, (E + eps)[:, None], theta0_deg=re["theta0_deg"])
    Tpe, Tle, _ = channels_from_z(fbe)
    chk("eps_cross сдвигает ОБА канала на себя, theta0 не трогает",
        abs(re["theta0_deg"] - t0) < 1e-10 and abs(Tpe[0] - (T_perp + eps)) < 1e-12
        and abs(Tle[0] - (T_par + eps)) < 1e-12,
        f"d(theta0) = {abs(re['theta0_deg'] - t0):.2e}, "
        f"|dT_perp| = {abs(Tpe[0] - (T_perp + eps)):.2e}")

    # --- 6. частотная зависимость: развёртка фазы и tau_leak
    nu = np.linspace(0.2, 2.0, 120)
    tau_true, d0_true, eta0_true, a_true = 0.19, 0.4, 0.02, 1.0
    Tp_nu = (0.95 - 0.1 * nu) * np.exp(-1j * 2 * np.pi * nu * 0.5)
    Tl_nu = Tp_nu * eta0_true * nu ** a_true * np.exp(
        -1j * (d0_true + 2 * np.pi * nu * tau_true))
    Y = np.stack([_synth_field(th, a, b, t0) for a, b in zip(Tp_nu, Tl_nu)], axis=1)

    fN = solve_bins(th, Y, order=2)
    rN = theta0_rank1(fN.z)
    fbN = solve_bins(th, Y, theta0_deg=rN["theta0_deg"])
    TpN, TlN, CovN = channels_from_z(fbN)
    chk("многочастотный случай: theta0 един по всем бинам",
        abs(rN["theta0_deg"] - t0) < 1e-9
        and abs(rN["lambda2_over_lambda1"]) < 1e-15,
        f"theta0 = {rN['theta0_deg']:.10f}, lambda2/lambda1 = "
        f"{rN['lambda2_over_lambda1']:.1e}")

    rr, vv = ratio_and_var(TpN, TlN, CovN)
    lr = logratio_regression(nu, rr, np.full(len(nu), 1e-12))
    chk("tau_leak из наклона развёрнутой фазы",
        abs(lr["tau_dphi_ps"] - tau_true) < 1e-9,
        f"{lr['tau_dphi_ps']:.12f} пс против {tau_true}")
    chk("степенной закон утечки eta0, a восстановлены",
        abs(lr["eta0"] - eta0_true) < 1e-9 and abs(lr["eta_exp"] - a_true) < 1e-9,
        f"eta0 = {lr['eta0']:.10f} (ждём {eta0_true}), a = {lr['eta_exp']:.10f}")

    # --- 6b. развёртка через дыры водяной маски
    valid = np.ones(len(nu), dtype=bool)
    valid[40:47] = False
    valid[80:84] = False
    lrm = logratio_regression(nu, rr, np.full(len(nu), 1e-12), valid=valid)
    chk("развёртка переживает разрывы частотной оси",
        abs(lrm["tau_dphi_ps"] - tau_true) < 1e-9,
        f"tau_dphi = {lrm['tau_dphi_ps']:.12f} пс, кусков = "
        f"{lrm['unwrap']['n_segments']}, сшивок = {lrm['unwrap']['n_joins']}")

    # --- 7. ковариация откалибрована (Монте-Карло)
    n_mc, sigma = 4000, 2e-3
    s2 = np.full((len(th), 1), sigma ** 2)
    Z = np.empty((n_mc, 3), dtype=complex)
    chi = np.empty(n_mc)
    for i in range(n_mc):
        n = (rng.normal(0, sigma, len(th)) + 1j * rng.normal(0, sigma, len(th))) / np.sqrt(2)
        ff = solve_bins(th, (E + n)[:, None], s2, order=2)
        Z[i] = ff.z[:, 0]
        chi[i] = ff.chi2[0]
    emp = np.cov(np.real(Z).T) + np.cov(np.imag(Z).T)     # = Cov(Re)+Cov(Im) = C
    theo = solve_bins(th, E[:, None], s2, order=2).C[0]
    rel = float(np.max(np.abs(emp - theo)) / np.max(np.abs(theo)))
    chk("ковариация inv(X^T W X) откалибрована (МК)", rel < 0.08,
        f"макс. относит. расхождение = {rel:.3f}")
    chi_red = float(np.mean(chi)) / (2 * len(th) - 6)
    chk("chi^2/dof = 1 при верной модели шума", abs(chi_red - 1.0) < 0.05,
        f"chi^2/dof = {chi_red:.4f}")
    pseudo = float(np.abs(np.mean((Z - Z.mean(0))[:, 1] ** 2)))
    chk("круговость сохранена (псевдоковариация ~ 0)",
        pseudo < 0.1 * float(np.mean(np.abs(Z[:, 1] - Z[:, 1].mean()) ** 2)),
        f"|псевдоков| = {pseudo:.2e}")

    # --- 8. ГЛАВНОЕ: положительная определённость и Коши-Буняковский
    #        Поле выполняет их ТОЖДЕСТВЕННО; мощность нарушает при низком SNR.
    rows = []
    for snr in (10.0, 3.6, 2.0):
        s = abs(T_par) / snr
        vf = vp = neg = 0
        n_it = 3000
        for _ in range(n_it):
            nz = (rng.normal(0, s, len(th)) + 1j * rng.normal(0, s, len(th))) / np.sqrt(2)
            Em = E + nz
            ff = solve_bins(th, Em[:, None], order=2)
            rr1 = theta0_rank1(ff.z)
            fbb = solve_bins(th, Em[:, None], theta0_deg=rr1["theta0_deg"])
            tp, tl, _ = channels_from_z(fbb)
            if abs(np.cos(np.angle(tp[0]) - np.angle(tl[0]))) > 1.0:
                vf += 1
            ch = _power_channels(th, np.abs(Em) ** 2)
            if not ch["cauchy_schwarz_ok"]:
                vp += 1
            if ch["I_par"] <= 0:
                neg += 1
        rows.append({"snr": snr, "field_viol_pct": 100.0 * vf / n_it,
                     "power_viol_pct": 100.0 * vp / n_it,
                     "power_neg_pct": 100.0 * neg / n_it})
        chk(f"поле: |cos dphi| <= 1 тождественно (SNR тени {snr})", vf == 0,
            f"нарушений поле {vf}/{n_it}, мощность {vp}/{n_it}, "
            f"|T_par|^2<0 в мощности {neg}/{n_it}")

    growth = rows[-1]["power_viol_pct"] > rows[0]["power_viol_pct"] + 5.0
    chk("мощность: доля нарушений РАСТЁТ с падением SNR (объясняет 46 % у series1)",
        growth,
        " -> ".join(f"SNR {r['snr']}: {r['power_viol_pct']:.1f} %" for r in rows))

    # --- 9. отрицательный контроль: вырождение НЕ снимается
    #        два разных (T_perp, T_par) дают разные данные, но одна и та же
    #        пара наблюдаемых (eta, dphi) неотличима от общего множителя
    scale = 0.7 * np.exp(1j * 0.9)
    Es = _synth_field(th, T_perp * scale, T_par * scale, t0)
    fs = solve_bins(th, Es[:, None], theta0_deg=t0)
    Tps, Tls, _ = channels_from_z(fs)
    chk("отрицательный контроль: eta и dphi инвариантны к общему множителю "
        "(вырождение D_eff<->eta НЕ снимается)",
        abs(abs(Tls[0] / Tps[0]) - eta_true) < 1e-12,
        f"eta(масштабированное) = {abs(Tls[0] / Tps[0]):.12f} против {eta_true:.12f}")

    n_ok = sum(c["ok"] for c in checks)
    if verbose:
        print(f"\n  Итог: {n_ok}/{len(checks)} проверок пройдено")
    return {"n_ok": n_ok, "n_total": len(checks), "checks": checks,
            "cauchy_schwarz_mc": rows}


def main():
    ap = argparse.ArgumentParser(description="Полевая линеаризация закона Малюса (A22)")
    ap.add_argument("--selftest", action="store_true", help="запустить самопроверку")
    ap.add_argument("--json", action="store_true", help="записать артефакт в results")
    args = ap.parse_args()

    if not args.selftest:
        ap.print_help()
        return 0

    print("A22 field_linearization — самопроверка")
    res = selftest(verbose=True)
    if args.json:
        RESULTS.mkdir(parents=True, exist_ok=True)
        out = RESULTS / "a22_field_selftest.json"
        payload = {k: v for k, v in res.items() if k != "checks"}
        payload["checks"] = [{kk: vv for kk, vv in c.items()} for c in res["checks"]]
        out.write_text(json.dumps(payload, ensure_ascii=False, indent=1),
                       encoding="utf-8")
        print(f"  артефакт: {out}")
    return 0 if res["n_ok"] == res["n_total"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
