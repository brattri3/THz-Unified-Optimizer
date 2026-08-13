# -*- coding: utf-8 -*-
u"""Линеаризованный фит закона Малюса: гармонический базис и замкнутые формулы.

Наблюдаемая нелинейна по УГЛУ, но линейна по КОЭФФИЦИЕНТАМ, если разложить её
по гармоникам угла. Тогда подгонка — одно нормальное уравнение: ни итераций, ни
стартовой точки, ни мультистарта, ни бутстрапа. Ковариация ``(X^T W X)^-1``
точна, потому что модель линейна. Это и позволяет считать фит на лету на старом
ПК без scipy.

Метод выведен и проверен сессией A (``research/two_wgp/NOTE_malus_linearization.md``,
задача ``A0_malus_linearization``, самопроверка 12/12). Здесь он **скопирован**, а
не импортирован: комплект для Windows 7 каталога ``research/`` не содержит.
Бит-в-бит сверка порта с оригиналом — ``tools/verify_fit_malus.py``.

Два порядка
-----------
``order=4`` — ТОЧНЫЙ для схемы этого проекта. Поле складывается когерентно,
``E = t_perp·cos²t + t_par·sin²t``, и только потом берётся ``U = |E|²``. Квадрат
модуля величины, линейной по ``cos2t``, есть тригонометрический многочлен ровно
4-го порядка:

    U(t) = c0 + a2·cos2t + b2·sin2t + a4·cos4t + b4·sin4t

``order=2`` — учебная 3-параметрическая форма. Она верна для ИНТЕНСИВНОСТНОЙ
(некогерентной) композиции и для нашей установки систематически неверна: теряет
четверть сигнала и полностью слепа к фазе утечки, потому что фаза входит ТОЛЬКО
через 4-ю гармонику. На диапазоне 0…90° это смещает угловой офсет до −12.6°.
Оставлена как демонстрация ловушки, а не как рабочий режим.

Частотная зависимость модели НЕ требует
---------------------------------------
На каждой частоте ``U(nu, t)`` — многочлен 4-го порядка по углу; сумма по
частотам сохраняет порядок. Поэтому ``T(t)`` фитируется тем же линейным МНК при
любом поведении ``t_perp(nu)`` и ``t_par(nu)``. Широкополосный фит даёт при этом
не «значение на неизвестной частоте», а среднее по спектру с весом мощности
референса — так его и подписывать. Сами спектры каналов даёт побинный фит.

Модуль не импортирует tkinter и matplotlib (проверяется приёмкой).
"""
from __future__ import annotations

import numpy as np

D2R = np.pi / 180.0
R2D = 180.0 / np.pi

WEIGHT_RELATIVE = "relative"
WEIGHT_NOISE = "noise"
WEIGHT_UNIFORM = "uniform"

WHICH_TIME = "time"
WHICH_BAND = "band"

#: минимальная доля от максимума, ниже которой вес обрезается — иначе точка с
#: околонулевым U получила бы бесконечный вес и одна определила бы весь фит
WEIGHT_FLOOR_REL = 1e-6


# --------------------------------------------------------------- базис и МНК
def harmonic_design(theta_deg, order=4):
    u"""Матрица плана X: столбцы [1, cos2t, sin2t, (cos4t, sin4t)].

    Модель U = X @ p строго ЛИНЕЙНА по p — отсюда весь замкнутый аппарат МНК.
    """
    t = np.asarray(theta_deg, dtype=np.float64) * D2R
    cols = [np.ones_like(t), np.cos(2 * t), np.sin(2 * t)]
    if order >= 4:
        cols += [np.cos(4 * t), np.sin(4 * t)]
    return np.column_stack(cols)


def fit_harmonics(theta_deg, U, sigma=None, order=4):
    u"""Взвешенный линейный МНК. -> (p, cov, info).

    cov = (X^T W X)^-1 — ТОЧНАЯ ковариация (модель линейна), бутстрап не нужен.
    Если sigma не задана, дисперсия оценивается из невязки: sigma² = RSS/(n−k).
    """
    X = harmonic_design(theta_deg, order)
    y = np.asarray(U, dtype=np.float64)
    n, k = X.shape
    if n < k:
        raise ValueError(u"углов %d < параметров %d: система недоопределена" % (n, k))

    if sigma is None:
        w = np.ones(n)
    else:
        s = np.asarray(sigma, dtype=np.float64)
        if np.any(s <= 0):
            raise ValueError(u"sigma должна быть > 0")
        w = 1.0 / s ** 2

    XtWX = np.dot(X.T, w[:, None] * X)
    XtWy = np.dot(X.T, w * y)
    p = np.linalg.solve(XtWX, XtWy)
    resid = y - np.dot(X, p)
    dof = n - k
    rss = float(np.sum(w * resid ** 2))

    if sigma is None:
        s2 = rss / dof if dof > 0 else np.nan
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
    u"""Учебный случай U = a0 + a1·cos2t + b1·sin2t -> контраст, утечка, офсет."""
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
        var_t0 = (b1 ** 2 * v_a1 + a1 ** 2 * v_b1 - 2 * a1 * b1 * c_ab) / (4 * A ** 4)
        out["theta0_err_deg"] = float(np.sqrt(max(var_t0, 0.0)) * R2D)
    return out


# --------------------------------------- замкнутые формулы (order=4, когерентно)
def theta0_error(p, cov, chi2_red=None):
    u"""sigma(theta0) из 2-й гармоники — аналитически, дельта-методом.

    При заданном chi2_red ковариация масштабируется на него (поправка Бирджа):
    это честно, когда известен не абсолютный уровень шума, а лишь его ФОРМА —
    ровно наш случай с весами sigma_i ~ U_i и неизвестным множителем.
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
    u"""Пять гармонических коэффициентов -> физика каналов, в замкнутом виде.

    В системе с офсетом: U = c0 + c2·cos2(t−t0) + c4·cos4(t−t0). Тождества:

        |t_perp|²              = c0 + c2 + c4
        |t_par|²               = c0 − c2 + c4
        Re(t_perp·conj(t_par)) = c0 − 3·c4

    Отсюда утечка и КОСИНУС относительной фазы — модельно-независимо, из одного
    линейного МНК. Фаза сидит только в 4-й гармонике, поэтому order=2 её теряет.

    ``c4`` берётся как ПРОЕКЦИЯ на систему офсета, а не как модуль ``hypot(a4,b4)``:
    у неё есть знак, и он входит в Re-часть.
    """
    c0, a2, b2, a4, b4 = (float(v) for v in p[:5])
    A2, A4 = float(np.hypot(a2, b2)), float(np.hypot(a4, b4))
    t0_h2 = 0.5 * np.arctan2(b2, a2) * R2D
    t0_h4 = 0.25 * np.arctan2(b4, a4) * R2D
    c2 = A2
    ph = 4 * t0_h2 * D2R
    c4 = a4 * np.cos(ph) + b4 * np.sin(ph)

    I_perp = c0 + c2 + c4
    I_par = c0 - c2 + c4
    cross = c0 - 3 * c4
    denom = np.sqrt(I_perp * I_par) if I_perp > 0 and I_par > 0 else np.nan
    cos_dphi = cross / denom if denom and np.isfinite(denom) else np.nan

    # два независимых офсета обязаны совпасть; сверка по модулю 45 градусов
    d = (t0_h4 - t0_h2 + 22.5) % 45.0 - 22.5

    ok = I_perp > 0 and I_par > 0
    return {
        "theta0_deg_from_h2": t0_h2,
        "theta0_deg_from_h4_mod45": t0_h4,
        "theta0_h4_minus_h2_deg": float(d),
        "A2": A2, "A4": A4, "c4_signed": float(c4),
        "harmonic4_over_harmonic2": A4 / A2 if A2 else np.nan,
        "I_perp": float(I_perp), "I_par": float(I_par),
        "eta_amplitude": float(np.sqrt(I_par / I_perp)) if ok else np.nan,
        "extinction_db": float(10 * np.log10(I_perp / I_par)) if ok else np.inf,
        "Re_cross": float(cross),
        "cos_dphi": float(cos_dphi),
        "cauchy_schwarz_ok": bool(np.isfinite(cos_dphi) and abs(cos_dphi) <= 1.0 + 1e-9),
    }


def channel_params_array(coef):
    u"""То же, но сразу для стопки коэффициентов ``(5, nбин)`` — для побинного фита.

    Возвращает словарь массивов. Отдельная реализация нужна ради скорости: цикл
    ``channel_params`` по двум сотням бинов свёл бы на нет смысл векторизации.
    """
    c0, a2, b2, a4, b4 = [np.asarray(x, dtype=np.float64) for x in coef[:5]]
    A2 = np.hypot(a2, b2)
    A4 = np.hypot(a4, b4)
    t0_h2 = 0.5 * np.arctan2(b2, a2) * R2D
    ph = 4 * t0_h2 * D2R
    c4 = a4 * np.cos(ph) + b4 * np.sin(ph)

    I_perp = c0 + A2 + c4
    I_par = c0 - A2 + c4
    cross = c0 - 3 * c4
    ok = (I_perp > 0) & (I_par > 0)
    with np.errstate(invalid="ignore", divide="ignore"):
        denom = np.where(ok, np.sqrt(np.abs(I_perp * I_par)), np.nan)
        cos_dphi = np.where(ok, cross / denom, np.nan)
        eta = np.where(ok, np.sqrt(np.abs(I_par / I_perp)), np.nan)
    return {
        "I_perp": I_perp, "I_par": I_par, "Re_cross": cross,
        "eta_amplitude": eta, "cos_dphi": cos_dphi,
        "theta0_deg_from_h2": t0_h2,
        "A2": A2, "A4": A4, "ok": ok,
    }


# ------------------------------------------------------ диагностика плана углов
def design_diagnostics(theta_deg, order=4, sigma_rel=0.01, U_scale=1.0, A2_rel=1.0):
    u"""Что угловой НАБОР (а не данные) делает с точностью офсета.

    Чистая геометрия плана эксперимента: число обусловленности и предсказанная
    sigma(theta0) при заданном относительном шуме. Именно здесь видно, что
    неполный трек в начале съёмки физику определить не может.
    """
    theta_deg = np.asarray(theta_deg, dtype=np.float64)
    X = harmonic_design(theta_deg, order)
    XtX = np.dot(X.T, X)
    cov = np.linalg.inv(XtX) * (sigma_rel * U_scale) ** 2
    a1, b1 = A2_rel * U_scale, 0.0            # без потери общности t0 = 0
    A = np.hypot(a1, b1)
    var_t0 = (b1 ** 2 * cov[1, 1] + a1 ** 2 * cov[2, 2]
              - 2 * a1 * b1 * cov[1, 2]) / (4 * A ** 4)
    n = len(theta_deg)
    G = XtX / n
    off24 = float(np.max(np.abs(G[1:3, 3:5]))) if order >= 4 else 0.0
    return {
        "n_angles": int(n),
        "n_distinct": int(len(np.unique(np.round(theta_deg, 6)))),
        "theta_min": float(np.min(theta_deg)), "theta_max": float(np.max(theta_deg)),
        "cond_XtX": float(np.linalg.cond(XtX)),
        "sigma_theta0_deg": float(np.sqrt(max(var_t0, 0.0)) * R2D),
        "max_overlap_h2_h4": off24,
    }


def high_harmonic_share(theta_deg, resid, sigma=None, order=4):
    u"""Доля мощности невязки, объяснимая гармониками 6 и 8. -> (доля, случайный уровень).

    Тест на выход за пределы модели: Фабри-Перо и нелинейность детектора рождают
    гармоники выше 4-й, которых у ``|E|²`` быть не может.

    Считается **в той же метрике, в которой минимизировался сам фит**, и по
    базису, ортогонализованному к базису фита. Оба условия обязательны, иначе
    статистика бессмысленна:

    * без весов доля считается там, где невязка не минимизировалась, и яркая
      зона (где абсолютная невязка велика) накачивает число. На эталонном
      прогоне невзвешенный счёт даёт 24.2 % против 6.1 % в метрике весов —
      разница между «есть выход за модель» и «нет»;
    * без ортогонализации гармоники 6 и 8 подбирают то, что уже описано
      базисом порядка 4, и доля завышается ещё раз.

    Вторым числом возвращается уровень, ожидаемый по чистой случайности:
    четыре параметра на ``dof`` степеней свободы. Без него доля не читается —
    сама по себе она не отвечает на вопрос «много это или мало».
    """
    t = np.asarray(theta_deg, dtype=np.float64) * D2R
    r = np.asarray(resid, dtype=np.float64)
    n = len(t)
    k = 5 if order >= 4 else 3
    X = np.column_stack([np.cos(6 * t), np.sin(6 * t), np.cos(8 * t), np.sin(8 * t)])
    dof = n - k
    if dof <= X.shape[1]:
        return float("nan"), float("nan")

    sw = np.ones(n) if sigma is None else 1.0 / np.asarray(sigma, dtype=np.float64)
    Xw = X * sw[:, None]
    X0 = harmonic_design(theta_deg, order) * sw[:, None]
    rw = r * sw
    tot = float(np.sum(rw ** 2))
    chance = float(X.shape[1]) / float(dof)
    if tot <= 0:
        return 0.0, chance
    try:
        Xw = Xw - np.dot(X0, np.linalg.lstsq(X0, Xw, rcond=None)[0])
        c = np.linalg.lstsq(Xw, rw, rcond=None)[0]
    except np.linalg.LinAlgError:
        return float("nan"), chance
    return float(np.sum(np.dot(Xw, c) ** 2) / tot), chance


# ------------------------------------------------------------------ веса
def weights_for(U, mode, noise=None):
    u"""Вектор sigma по выбранному режиму. -> (sigma, пометка или None).

    ``relative`` — постоянная ОТНОСИТЕЛЬНАЯ ошибка, sigma_i ~ U_i. Рецепт сессии A
    и рабочий режим. Без него невзвешенный МНК взвешивает точку как U², яркая
    зона полностью заглушает тень, а вся информация об утечке сидит именно в
    тени: у A на series1 такая подгонка дала |t_par|² < 0.

    ``noise`` — sigma из измеренного шумового пола. Физически честнее в тени, но
    в глубоком гашении это оценка, а не измерение.

    ``uniform`` — равные веса. Только для демонстрации ловушки.
    """
    U = np.asarray(U, dtype=np.float64)
    top = float(np.max(np.abs(U))) if U.size else 1.0
    floor = max(top * WEIGHT_FLOOR_REL, 1e-300)
    if mode == WEIGHT_UNIFORM:
        return np.ones_like(U), None
    if mode == WEIGHT_NOISE:
        if noise is None:
            return np.maximum(np.abs(U), floor), u"шумовой пол недоступен, веса взяты относительные"
        s = np.asarray(noise, dtype=np.float64)
        bad = ~np.isfinite(s) | (s <= 0)
        if np.all(bad):
            return np.maximum(np.abs(U), floor), u"шумовой пол недоступен, веса взяты относительные"
        if np.any(bad):
            s = s.copy()
            s[bad] = np.nanmax(s[~bad])
        return np.maximum(s, floor), None
    return np.maximum(np.abs(U), floor), None


# ------------------------------------------------------------ угловой фит
class AngularFit(object):
    u"""Результат подгонки одной угловой зависимости по всему треку."""

    __slots__ = ("which", "order", "weight_mode", "theta", "U", "sigma",
                 "p", "cov", "info", "params", "theta0_err", "diag",
                 "used", "n_excluded", "hh_share", "hh_chance", "notes")

    def __init__(self):
        self.which = WHICH_TIME
        self.order = 4
        self.weight_mode = WEIGHT_RELATIVE
        self.theta = None
        self.U = None
        self.sigma = None
        self.p = None
        self.cov = None
        self.info = None
        self.params = None       # channel_params (order=4) либо malus_params (order=2)
        self.theta0_err = float("nan")
        self.diag = None
        self.used = None         # индексы точек трека, попавших в фит
        self.n_excluded = 0
        self.hh_share = float("nan")
        self.hh_chance = float("nan")
        self.notes = []

    def curve(self, theta_grid):
        u"""Непрерывная подогнанная кривая на произвольной сетке углов."""
        return np.dot(harmonic_design(theta_grid, self.order), self.p)

    def residuals(self):
        u"""Невязки в порядке ``used``."""
        return self.info["resid"]


def fit_angular(theta_deg, U, weight_mode=WEIGHT_RELATIVE, order=4, noise=None,
                which=WHICH_TIME):
    u"""Подгонка ``U(theta)`` гармоническим базисом. -> AngularFit.

    Бросает ValueError, если различных углов меньше числа параметров: это не
    сбой, а физика — угловой скан на одной частоте ограничивает не более пяти
    вещественных чисел, и меньшим числом углов их не определить.
    """
    theta_deg = np.asarray(theta_deg, dtype=np.float64)
    U = np.asarray(U, dtype=np.float64)
    k = 5 if order >= 4 else 3
    n_distinct = len(np.unique(np.round(theta_deg, 6)))
    if n_distinct < k:
        raise ValueError(u"различных углов %d, а параметров %d — фит невозможен"
                         % (n_distinct, k))

    res = AngularFit()
    res.which = which
    res.order = int(order)
    res.weight_mode = weight_mode
    res.theta = theta_deg
    res.U = U
    sigma, note = weights_for(U, weight_mode, noise)
    res.sigma = sigma
    if note:
        res.notes.append(note)

    p, cov, info = fit_harmonics(theta_deg, U, sigma=sigma, order=order)
    res.p, res.cov, res.info = p, cov, info
    res.params = channel_params(p, cov) if order >= 4 else malus_params(p, cov)
    res.theta0_err = theta0_error(p, cov, info["chi2_red"])
    res.diag = design_diagnostics(theta_deg, order)
    res.hh_share, res.hh_chance = high_harmonic_share(
        theta_deg, info["resid"], sigma, order)
    return res


# --------------------------------------------------------- частотный фит
class SpectralFit(object):
    u"""Побинный фит: тот же базис, применённый к каждой частоте отдельно.

    Даёт спектры каналов |t_perp|²(nu) и |t_par|²(nu) модельно-независимо —
    бины подгоняются независимо друг от друга, никакой модели частотной
    зависимости для этого не нужно.
    """

    __slots__ = ("freqs", "coef", "order", "weight_mode", "ch",
                 "n_bins", "n_neg", "n_cs_violations", "used", "notes")

    def __init__(self):
        self.freqs = None
        self.coef = None         # (5, nбин)
        self.order = 4
        self.weight_mode = WEIGHT_RELATIVE
        self.ch = None           # словарь массивов из channel_params_array
        self.n_bins = 0
        self.n_neg = 0
        self.n_cs_violations = 0
        self.used = None
        self.notes = []

    def transmission(self, theta_deg):
        u"""Восстановленная фитом T(nu) для заданного угла.

        Она сглажена по построению: каждая точка получена по ВСЕМ трассам трека,
        а не по одной паре образец/референс.
        """
        x = harmonic_design(np.atleast_1d(np.asarray(theta_deg, dtype=np.float64)),
                            self.order)[0]
        return np.dot(x, self.coef)


def fit_spectral(theta_deg, Y, freqs, weight_mode=WEIGHT_RELATIVE, order=4,
                 noise=None, common_sigma=None):
    u"""Фит по каждому частотному бину сразу. -> SpectralFit.

    ``Y`` — матрица ``(углов, бинов)`` со значениями ``T(nu, theta)``.

    Все бины решаются одним ``np.linalg.solve`` над стопкой матриц ``(бин,5,5)``:
    двести бинов обходятся в единицы миллисекунд, поэтому фит можно считать на
    каждое обновление.

    ``common_sigma`` задаёт ОДИН вектор угловых весов на все бины. Тогда, по
    линейности МНК, спектрально-взвешенное среднее побинных коэффициентов точно
    равно коэффициентам широкополосного фита — это тождество проверяется приёмкой.
    """
    theta_deg = np.asarray(theta_deg, dtype=np.float64)
    Y = np.asarray(Y, dtype=np.float64)
    k = 5 if order >= 4 else 3
    n_distinct = len(np.unique(np.round(theta_deg, 6)))
    if n_distinct < k:
        raise ValueError(u"различных углов %d, а параметров %d — фит невозможен"
                         % (n_distinct, k))

    res = SpectralFit()
    res.order = int(order)
    res.weight_mode = weight_mode
    res.freqs = np.asarray(freqs, dtype=np.float64)
    res.n_bins = int(Y.shape[1])

    X = harmonic_design(theta_deg, order)                      # (nang, k)
    if common_sigma is not None:
        s = np.asarray(common_sigma, dtype=np.float64)[:, None]
        W = np.repeat(1.0 / s, Y.shape[1], axis=1)
    else:
        top = np.max(np.abs(Y), axis=0)
        floor = np.maximum(top * WEIGHT_FLOOR_REL, 1e-300)
        if weight_mode == WEIGHT_UNIFORM:
            W = np.ones_like(Y)
        elif weight_mode == WEIGHT_NOISE and noise is not None:
            s = np.asarray(noise, dtype=np.float64)
            s = np.where(np.isfinite(s) & (s > 0), s, np.nan)
            if np.all(np.isnan(s)):
                W = 1.0 / np.maximum(np.abs(Y), floor)
                res.notes.append(u"шумовой пол недоступен, веса взяты относительные")
            else:
                s = np.where(np.isnan(s), np.nanmax(s), s)
                W = 1.0 / np.maximum(s, floor)
        else:
            W = 1.0 / np.maximum(np.abs(Y), floor)

    Xw = X[:, :, None] * W[:, None, :]                         # (nang, k, nbin)
    A = np.einsum("akb,ajb->bkj", Xw, Xw)
    r = np.einsum("akb,ab->bk", Xw, Y * W)
    # правая часть подаётся столбцом: без явного [..., None] трактовка формы
    # (n, k) разъезжается между numpy 1.21 и 2.x
    res.coef = np.linalg.solve(A, r[..., None])[..., 0].T      # (k, nbin)

    if order >= 4:
        ch = channel_params_array(res.coef)
        res.ch = ch
        res.n_neg = int(np.sum(~ch["ok"]))
        cs = ch["cos_dphi"]
        res.n_cs_violations = int(np.sum(np.isfinite(cs) & (np.abs(cs) > 1.0 + 1e-9)))
    return res


# ---------------------------------------------------- связка с треком прогона
#: ниже этого числа различных углов фит не выполняется вовсе
MIN_DISTINCT_ANGLES = 5
#: выше этой обусловленности плана углов результат помечается ненадёжным
COND_WARN = 50.0


class TrackFit(object):
    u"""Всё, что подгонка говорит о треке: две угловые кривые и спектр каналов."""

    __slots__ = ("time", "band", "spectral", "order", "weight_mode",
                 "n_points", "n_used", "n_excluded", "notes")

    def __init__(self):
        self.time = None        # AngularFit по T во времени
        self.band = None        # AngularFit по T в полосе (при включённом Парсевале)
        self.spectral = None    # SpectralFit
        self.order = 4
        self.weight_mode = WEIGHT_RELATIVE
        self.n_points = 0
        self.n_used = 0
        self.n_excluded = 0
        self.notes = []

    @property
    def ok(self):
        return self.time is not None

    def delta_curve(self, theta_grid):
        u"""Разность двух подогнанных функций: фит(время) − фит(полоса).

        Непрерывный аналог точек Δ на среднем ряду панелей. Показывает, сколько
        пропускания набирается за пределами полосы интегрирования.
        """
        if self.time is None or self.band is None:
            return None
        return self.time.curve(theta_grid) - self.band.curve(theta_grid)

    def summary_lines(self):
        u"""Блок фита для панели «измерение» и для консоли.

        Символы, которых нет в кодировке консоли Windows, здесь не используются:
        строка обязана читаться и в окне, и в `.bat`-е у прибора.
        """
        out = []
        if not self.ok:
            for n in self.notes:
                out.append(u"фит не выполнен: %s" % n)
            return out

        f = self.time
        pr = f.params
        order = u"порядок %d" % f.order
        wm = {WEIGHT_RELATIVE: u"веса относительные",
              WEIGHT_NOISE: u"веса по шумовому полу",
              WEIGHT_UNIFORM: u"веса равные"}.get(f.weight_mode, f.weight_mode)
        out.append(u"ФИТ Малюса: %s, %s, точек %d из %d, cond %.1f"
                   % (order, wm, self.n_used, self.n_points, f.diag["cond_XtX"]))
        if f.order >= 4:
            out.append(u"  |t_перп|^2 = %.6g   |t_прод|^2 = %.4g   экстинкция %.2f дБ"
                       % (pr["I_perp"], pr["I_par"], pr["extinction_db"]))
            out.append(u"  утечка eta = %.5g   cos(разн. фаз) = %+.4f%s"
                       % (pr["eta_amplitude"], pr["cos_dphi"],
                          u"" if pr["cauchy_schwarz_ok"] else u"  [!] |cos| > 1"))
            out.append(u"  угол столика theta0 = %+.4f +- %.4f град   "
                       u"(из 4-й гармоники %+.4f, расхождение %+.4f)"
                       % (pr["theta0_deg_from_h2"], f.theta0_err,
                          pr["theta0_deg_from_h4_mod45"], pr["theta0_h4_minus_h2_deg"]))
            out.append(u"  A4/A2 = %.4f (идеал cos^4 даёт 0.25); невязка %.3g; "
                       u"гармоники 6-8 несут %.1f %% невязки при случайных %.1f %%"
                       % (pr["harmonic4_over_harmonic2"], f.info["resid_rms"],
                          100.0 * f.hh_share, 100.0 * f.hh_chance))
        else:
            out.append(u"  экстинкция %.2f дБ   утечка eta = %.5g   theta0 = %+.4f град"
                       % (pr["extinction_db"], pr["eta_amplitude"], pr["theta0_deg"]))
        if self.band is not None and self.band.order >= 4:
            b = self.band.params
            out.append(u"  по полосе: |t_перп|^2 = %.6g   |t_прод|^2 = %.4g   eta = %.5g"
                       % (b["I_perp"], b["I_par"], b["eta_amplitude"]))
        sp = self.spectral
        if sp is not None and sp.ch is not None:
            out.append(u"  по частотам: бинов %d, |t_прод|^2 < 0 в %d, "
                       u"нарушений Коши-Буняковского %d"
                       % (sp.n_bins, sp.n_neg, sp.n_cs_violations))
        for n in self.notes:
            out.append(u"  [!] %s" % n)
        for n in f.notes:
            out.append(u"  [!] %s" % n)
        return out

    def point_line(self, index):
        u"""Строка об отклонении конкретной точки от подогнанной кривой.

        Прямой детектор сбойной трассы во время съёмки: точка, ушедшая от кривой
        на несколько своих сигм, обычно означает не физику, а сдвиг образца.
        """
        f = self.time
        if f is None or f.used is None:
            return None
        where = np.nonzero(np.asarray(f.used) == int(index))[0]
        if not len(where):
            return u"эта точка в фит не вошла"
        i = int(where[0])
        r = float(f.info["resid"][i])
        s = float(f.sigma[i])
        rel = r / s if s > 0 else float("nan")
        return (u"отклонение от фита %+.4g п.п. (%+.2f сигмы своего веса)"
                % (100.0 * r, rel))


def _mean_noise(point):
    u"""Скалярная оценка шумового пола точки в единицах пропускания."""
    tn = getattr(point, "T_noise", None)
    if tn is None:
        return float("nan")
    tn = np.asarray(tn, dtype=np.float64)
    if not tn.size:
        return float("nan")
    good = np.isfinite(tn) & (tn > 0)
    if not np.any(good):
        return float("nan")
    return float(np.mean(tn[good]))


def _common_grid(points, idx):
    u"""Общая частотная сетка для побинного фита.

    У точек она может различаться: маска динамического диапазона строится по
    СВОЁМУ референсу, а трассы бывают разной длины. Берём пересечение — молча
    подставить чужую сетку было бы подменой данных.
    """
    common = None
    for i in idx:
        f = getattr(points[i], "freqs", None)
        if f is None or not len(f):
            return None
        common = f if common is None else np.intersect1d(common, f)
        if not len(common):
            return None
    return common


def fit_track(points, settings):
    u"""Подгонка всего трека по текущим настройкам. -> TrackFit."""
    res = TrackFit()
    res.n_points = len(points)
    order = int(getattr(settings, "fit_order", 4))
    weight_mode = getattr(settings, "fit_weights", WEIGHT_RELATIVE)
    drop_warned = bool(getattr(settings, "fit_drop_warned", False))
    res.order = order
    res.weight_mode = weight_mode

    idx = []
    for i, p in enumerate(points):
        if not np.isfinite(p.T) or p.T <= 0:
            continue
        if drop_warned and p.warnings:
            continue
        idx.append(i)
    res.n_used = len(idx)
    res.n_excluded = res.n_points - res.n_used
    if drop_warned and res.n_excluded:
        res.notes.append(u"исключено точек с предупреждениями: %d" % res.n_excluded)

    if not idx:
        res.notes.append(u"нет пригодных точек")
        return res

    theta = np.array([points[i].trace.angle for i in idx], dtype=np.float64)
    n_distinct = len(np.unique(theta))
    need = 5 if order >= 4 else 3
    if n_distinct < max(need, MIN_DISTINCT_ANGLES if order >= 4 else need):
        res.notes.append(u"снято %d различных углов, для порядка %d нужно не менее %d"
                         % (n_distinct, order, max(need, MIN_DISTINCT_ANGLES
                                                   if order >= 4 else need)))
        return res

    noise = np.array([_mean_noise(points[i]) for i in idx], dtype=np.float64)

    U_time = np.array([points[i].T for i in idx], dtype=np.float64)
    try:
        res.time = fit_angular(theta, U_time, weight_mode, order, noise, WHICH_TIME)
    except (ValueError, np.linalg.LinAlgError) as exc:
        res.notes.append(u"угловой фит не сошёлся: %s" % exc)
        return res
    res.time.used = idx
    res.time.n_excluded = res.n_excluded
    if res.time.diag["cond_XtX"] > COND_WARN:
        res.notes.append(u"план углов неполон: cond = %.0f, оценки ненадёжны"
                         % res.time.diag["cond_XtX"])
    if order >= 4 and res.time.params["I_par"] <= 0:
        res.notes.append(u"|t_прод|^2 вышло отрицательным — смените веса "
                         u"или доснимите зону гашения")

    if getattr(settings, "parseval_on", False):
        U_band = np.array([points[i].T_freq for i in idx], dtype=np.float64)
        if np.all(np.isfinite(U_band)) and np.all(U_band > 0):
            try:
                res.band = fit_angular(theta, U_band, weight_mode, order, noise,
                                       WHICH_BAND)
                res.band.used = idx
            except (ValueError, np.linalg.LinAlgError) as exc:
                res.notes.append(u"фит по полосе не сошёлся: %s" % exc)

    if getattr(settings, "fit_spectral_on", False):
        grid = _common_grid(points, idx)
        if grid is None or len(grid) < 8:
            res.notes.append(u"общих частотных бинов у трасс меньше восьми — "
                             u"побинный фит не выполнен")
        else:
            Y = np.empty((len(idx), len(grid)), dtype=np.float64)
            N = np.empty_like(Y)
            for row, i in enumerate(idx):
                p = points[i]
                pos = np.searchsorted(p.freqs, grid)
                Y[row] = np.asarray(p.T_nu, dtype=np.float64)[pos]
                tn = getattr(p, "T_noise", None)
                N[row] = (np.asarray(tn, dtype=np.float64)[pos]
                          if tn is not None else np.nan)
            try:
                res.spectral = fit_spectral(theta, Y, grid, weight_mode, order, N)
                res.spectral.used = idx
            except (ValueError, np.linalg.LinAlgError) as exc:
                res.notes.append(u"побинный фит не сошёлся: %s" % exc)
    return res
