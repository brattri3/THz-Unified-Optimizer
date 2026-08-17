"""Быстрый прототип: угловая и спектральная кривая аттенюатора по измерениям.

Читает серию `att-11-16-s2`/`s3` (`data_pool/two_wgp_attenuator`, первый WGP
вращается, второй статичен вертикально — см. README этого каталога), считает:

  1. пропускание во времени T(theta) на каждый угол -> линеаризованный фит
     закона Малюса (`track_viewer.core.fit_malus`, порядок 4, замкнутый МНК);
  2. модельно-независимые спектральные каналы I_perp(nu), I_par(nu)
     (`fit_spectral` того же модуля) -> подгонка ЧИСТОЙ теории Бланко
     (`attenuator_app.core.blanco.blanco_t`, только период P и диаметр D, без
     феноменологии потерь/утечки) в энергетическом представлении.

Прототип для отчёта — минимум проверок, без сверки с зоной A (санкция
владельца 2026-08-17): не тянуть время до появления её результатов на
`s2`/`s3`.

Запуск из корня репозитория:
    .venv\\Scripts\\python.exe -m attenuator_app.tools.measured_curve
    .venv\\Scripts\\python.exe -m attenuator_app.tools.measured_curve --dataset att-11-16-s3
"""
from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))

from track_viewer.core.compat import rfft_freqs, trapezoid          # noqa: E402
from track_viewer.core.fit_malus import (WEIGHT_RELATIVE, fit_angular,  # noqa: E402
                                          fit_spectral)
from track_viewer.core import fit_malus as fm                       # noqa: E402
from track_viewer.core import physics as ph                         # noqa: E402
from track_viewer.core.scan import Scan                              # noqa: E402
from attenuator_app.core.blanco import blanco_t, dressed_t           # noqa: E402

#: НАЧАЛЬНОЕ приближение для показателя степенного закона потерь (паспортный
#: номинал, attenuator_app/tools/make_passport.py:SPEC["gamma"]) и дефолт для
#: мест, где gamma не передаётся явно. Раньше был буквально зафиксирован, но
#: владелец 2026-08-17 указал на несовпадение наклона модели с данными даже
#: после подгонки P/D/потерь -- теперь тоже свободный параметр `fit_pure_blanco`.
GAMMA_FIXED = 1.58

DATA_ROOT = REPO / "data_pool" / "two_wgp_attenuator"
FNAME_RE = re.compile(r"^(?P<dataset>.+)_(?P<angle>[+-]?\d+)deg_rep\d+_(?P<kind>sig|bg)\.txt$")


def normalize_angle(a: float) -> float:
    """0..360 -> -180..180, чтобы 270 град читался как -90."""
    a = a % 360.0
    return a - 360.0 if a > 180.0 else a


def load_txt(path: Path):
    """Время (пс) + поле; лишние колонки (счётчик подвижки) игнорируются."""
    raw = np.loadtxt(path)
    return raw[:, 0].astype(float), raw[:, 1].astype(float)


def scan_series(root: Path, dataset: str) -> dict:
    """angle_deg -> (t_sig, E_sig, t_bg, E_bg). Sig/bg уже связаны по имени файла.

    Старая плоская конвенция `{dataset}_{angle}deg_rep{N}_{sig|bg}.txt` — один
    bg жёстко на угол, история выбора референса потеряна при переименовании.
    """
    found: dict = {}
    for f in root.glob(f"{dataset}_*deg_rep*_*.txt"):
        m = FNAME_RE.match(f.name)
        if not m or m.group("dataset") != dataset:
            continue
        angle = normalize_angle(float(m.group("angle")))
        found.setdefault(angle, {})[m.group("kind")] = f
    out = {}
    for angle, kv in found.items():
        if "sig" not in kv or "bg" not in kv:
            continue
        t_s, E_s = load_txt(kv["sig"])
        t_b, E_b = load_txt(kv["bg"])
        out[angle] = (t_s, E_s, t_b, E_b)
    return out




def baseline_correct(t: np.ndarray, E: np.ndarray, frac: float = 0.15) -> np.ndarray:
    """Линейный (не константный) тренд базовой линии по пред-/постимпульсу.

    Обычное вычитание среднего по предымпульсу ловит только DC-смещение;
    дрейф от движения подвижки линии задержки и неидеальной юстировки луча
    вдоль оптической оси даёт ЛИНЕЙНЫЙ наклон в пределах трассы, который
    константа не убирает.
    """
    n = len(E)
    k = max(4, int(round(frac * n)))
    t_edge = np.concatenate([t[:k], t[-k:]])
    E_edge = np.concatenate([E[:k], E[-k:]])
    slope, intercept = np.polyfit(t_edge, E_edge, 1)
    return E - (intercept + slope * t)


def time_transmission(t_s, E_s, t_b, E_b) -> float:
    Es = baseline_correct(t_s, E_s)
    Eb = baseline_correct(t_b, E_b)
    num = trapezoid(Es ** 2, t_s)
    den = trapezoid(Eb ** 2, t_b)
    return num / den if den > 0 else float("nan")


def spectrum_ratio(t_s, E_s, t_b, E_b):
    """(freqs_thz, T(nu)) = |rFFT(sig)|^2 / |rFFT(bg)|^2, после коррекции базовой линии."""
    Es = baseline_correct(t_s, E_s)
    Eb = baseline_correct(t_b, E_b)
    n = min(len(Es), len(Eb))
    dt = float(t_s[1] - t_s[0])
    freqs = rfft_freqs(n, dt)
    Ps = np.abs(np.fft.rfft(Es[:n])) ** 2
    Pb = np.abs(np.fft.rfft(Eb[:n])) ** 2
    with np.errstate(divide="ignore", invalid="ignore"):
        T_nu = Ps / Pb
    return freqs, T_nu


def reference_power_spectrum(bg_traces):
    """Средний спектр мощности референсных (bg) импульсов -- |Ê_ref(nu)|^2.

    Это и есть аппаратная функция всей системы (источник+оптика+детектор) --
    не теория, а прямое измерение. Используется как ВЕС при частотном
    усреднении модельной кривой (`blanco_angular_curve`), чтобы оно было
    согласовано с теоремой Парсеваля: `T_time = ∫T(nu)|Ê_ref(nu)|^2 dnu /
    ∫|Ê_ref(nu)|^2 dnu` (полная полоса, `physics.py`), а не плоское среднее по
    узкой полосе фита каналов. Вес сам обнуляется там, где у импульса нет
    энергии -- отдельно урезать диапазон не нужно (владелец 2026-08-19).

    `bg_traces` -- список пар `(t, E)`. Частотные сетки усредняемых трасс
    могут чуть отличаться по длине (разные дни/повторы) -- берём минимальную
    длину, обрезая остальные (низкие частоты у всех совпадают, отличается
    только хвост).
    """
    specs = []
    for t, E in bg_traces:
        Ec = baseline_correct(np.asarray(t, dtype=float), np.asarray(E, dtype=float))
        dt = float(t[1] - t[0])
        f = rfft_freqs(len(Ec), dt)
        P = np.abs(np.fft.rfft(Ec)) ** 2
        specs.append((f, P))
    m = min(len(f) for f, _ in specs)
    freqs_ref = specs[0][0][:m]
    power_ref = np.mean([P[:m] for _, P in specs], axis=0)
    # DC-бин (nu=0) исключаем: модель Бланко там не определена (деление на
    # частоту в переходе к длине волны, core/blanco.py) -- физически это не
    # ТГц-сигнал, а остаток смещения после baseline_correct.
    keep = freqs_ref > 0
    return freqs_ref[keep], power_ref[keep]


def db(x):
    with np.errstate(divide="ignore", invalid="ignore"):
        return 10.0 * np.log10(np.asarray(x, dtype=float))


def channels_at_theta0(spec_fit, theta0_deg):
    """I_perp(nu), I_par(nu) при ЗАДАННОМ офсете ротатора theta0 (не авто по бину).

    `SpectralFit.ch["I_perp"/"I_par"]` (channel_params_array в fit_malus.py)
    берёт офсет theta0 из 2-й гармоники КАЖДОГО бина отдельно -- разумный
    дефолт, но не даёт проверить, насколько разложение t_perp/t_par
    чувствительно к неточности установки поляризатора в ротаторе (вопрос
    владельца 2026-08-17). Здесь -- один и тот же ручной theta0 для всех
    бинов, пересчитано из сырых гармонических коэффициентов
    (`spec_fit.coef`, форма (5, n_бинов)) по тем же тождествам, что и
    `fit_malus.channel_params_array`, но без переопределения фазы по данным.
    """
    c0, a2, b2, a4, b4 = (spec_fit.coef[0], spec_fit.coef[1], spec_fit.coef[2],
                          spec_fit.coef[3], spec_fit.coef[4])
    A2 = np.hypot(a2, b2)
    ph = 4.0 * np.deg2rad(theta0_deg)
    c4 = a4 * np.cos(ph) + b4 * np.sin(ph)
    return c0 + A2 + c4, c0 - A2 + c4


def blanco_channel_model(freqs, P, D, loss_factor, gamma=GAMMA_FIXED):
    """Модель каналов I_perp(nu), I_par(nu) -- ПОЛНАЯ Джонс-матричная схема S1,
    оба WGP считаются ИДЕНТИЧНЫМИ (одна и та же геометрия Бланко P,D,потери,
    gamma; владелец 2026-08-17 подтвердил -- по паспорту так и есть).

    Вывод (лаб. система отсчёта: ось x = ось пропускания статичного WGP2 =
    ось входной поляризации = ось детектора, theorem of alignment; WGP1 на
    угле theta = theta_reading - theta0):

        J(phi) = [[tp*c^2+ta*s^2, (tp-ta)*c*s], [(tp-ta)*c*s, tp*s^2+ta*c^2]]
        M = J(0) . J(theta)                         # WGP2 статичен при phi=0
        E = M . (1, 0)^T = (tp*(tp*c^2+ta*s^2), ta*(tp-ta)*c*s)
        I(theta,nu) = |E_x|^2 = |tp|^2 * |tp*cos^2(theta) + ta*sin^2(theta)|^2

    Каналы `I_perp`/`I_par`, которые извлекает `fit_spectral`/`channels_at_theta0`
    из ИЗМЕРЕННОГО T(theta,nu), -- это как раз значения этой формулы при
    theta=0 (тета=theta0) и theta=90 град:

        I_perp = |tp|^2 * |tp|^2 = |tp|^4
        I_par  = |tp|^2 * |ta|^2

    А НЕ `|tp|^2`, `|ta|^2` одного WGP -- та более простая формула стояла
    здесь раньше и подгоняла (P,D,потери) под ПРОИЗВЕДЕНИЕ двух одинаковых
    WGP как будто это один WGP; из-за этого угловая кривая в абсолютном
    режиме проседала на лишний множитель (найдено владельцем 2026-08-17:
    пик ~0.8 против ~0.9 у данных при отличном совпадении самих каналов --
    ровно то расхождение, которого и следовало ждать от `|tp|^2` вместо
    `|tp|^4`). Множитель потерь `exp(-(loss_factor/4.343)*nu^gamma)` --
    общий для tp/ta, без eta0/tau_leak (спорная феноменология зоны A, сюда
    не берём, `core/blanco.py:dressed_t`).
    """
    tp, ta, _ = dressed_t(freqs, P, D, loss_factor=loss_factor, gamma=gamma)
    tp_e, ta_e = np.abs(tp) ** 2, np.abs(ta) ** 2
    return tp_e ** 2, tp_e * ta_e


def blanco_angular_curve(theta_deg, freqs, P, D, loss_factor, theta0_deg=0.0,
                         gamma=GAMMA_FIXED, relative=True, weight=None):
    """T(theta), ВОССТАНОВЛЕННОЕ из ПОЛНОЙ Джонс-матричной схемы S1, два
    ИДЕНТИЧНЫХ WGP (вывод и обозначения -- `blanco_channel_model`, тот же P,D,
    потери,gamma для обоих):

        I(theta, nu) = |tp(nu)|^2 * |tp(nu)*cos^2(theta-theta0) + ta(nu)*sin^2(theta-theta0)|^2

    `tp`, `ta` -- ИСТИННЫЕ параметры ОДНОГО WGP (подогнаны в `fit_pure_blanco`
    к `blanco_channel_model`, а не к сырым `|tp|^2`/`|ta|^2` -- см. её докстрок
    про задвоение/недобор, найденный владельцем 2026-08-17 по проседанию пика
    в абсолютном режиме).

    `weight=None` (по умолчанию) -- ПЛОСКОЕ среднее по `freqs`, как раньше.
    `weight=power_ref(freqs)` -- среднее ВЗВЕШЕНО спектральной плотностью
    референсного импульса (`reference_power_spectrum`) -- согласовано с
    теоремой Парсеваля, которой подчиняется T_time (`physics.py`:
    `T_time = ∫T(nu)|Ê_ref(nu)|^2 dnu / ∫|Ê_ref(nu)|^2 dnu`, полная полоса).
    Плоское среднее по узкой полосе фита (0.3-1.5 ТГц) и взвешенное по всей
    записанной полосе -- РАЗНЫЕ операции, они не обязаны давать одно число;
    расхождение пика ~5-6% в абсолютном режиме, найденное владельцем
    2026-08-19, было объяснено именно этим, не качеством фита (проверено:
    увеличение веса I_perp в целевой функции пик не двигает, только портит
    I_par). `freqs` тогда должен быть ШИРЕ, чем полоса фита -- сам вес
    обнулится там, где у импульса нет энергии, отдельно урезать не нужно.

    `relative=True` (по умолчанию) — нормировка на theta0 (T(theta0)=1 по
    построению). Множитель `|tp|^2` -- общий для всех theta (не зависит от
    угла), поэтому в этом режиме сокращается ПРИБЛИЖЁННО (точно -- только
    поточечно по частоте, до усреднения; после усреднения по-разному
    взвешенных `|tp|^2 E1^2` и `|tp|^4` строгого равенства ratio-of-means нет,
    но при плавном `loss_factor` разница мала).
    `relative=False` — АБСОЛЮТНОЕ пропускание всей сборки, без деления на U0;
    `bg` в этих сериях -- пучок БЕЗ аттенюатора, так что измеренные точки и
    фит Малюса и без того на абсолютной шкале.
    """
    tp, ta, _ = dressed_t(freqs, P, D, loss_factor=loss_factor, gamma=gamma)
    d = np.deg2rad(np.asarray(theta_deg, dtype=float) - theta0_deg)
    c2, s2 = np.cos(d) ** 2, np.sin(d) ** 2
    E1 = tp[None, :] * c2[:, None] + ta[None, :] * s2[:, None]

    def wavg(x):
        return np.mean(x, axis=-1) if weight is None else np.average(x, axis=-1, weights=weight)

    if relative:
        U = wavg(np.abs(E1) ** 2)
        U0 = wavg(np.abs(tp) ** 2)
        return U / U0 if U0 > 0 else U
    return wavg(np.abs(tp[None, :]) ** 2 * np.abs(E1) ** 2)


def fit_pure_blanco(freqs, I_perp, I_par, p0=(16.0, 11.0, 0.25, GAMMA_FIXED)):
    """Подгонка ПОЛНОЙ Джонс-матричной модели (два ИДЕНТИЧНЫХ WGP, S1) к
    измеренным энергетическим каналам, `blanco_channel_model` (не сырые
    `|tp|^2`/`|ta|^2` одного WGP -- задвоение/недобор, см. её докстрок).

    Четыре свободных параметра: период `P`, диаметр `D`, `loss_factor`
    (дБ/ТГц^gamma по мощности) и сам показатель `gamma` -- раньше был
    зафиксирован на паспортном номинале (`GAMMA_FIXED`), но владелец
    2026-08-17 заметил, что НАКЛОН модельной кривой не совпадает с данными
    даже после подгонки P/D/потерь — это симптом неверного показателя
    степени, не только амплитуды потерь. Проектный диапазон для gamma —
    0.2…4.2 (research/SYNTHESIS.md, разброс по конфигурациям Фазы 3).
    Без `eta0`/`tau_leak` — это отдельная утечка/фаза, спорная феноменология
    зоны A (`research/two_wgp/NOTE_stage2_blanco.md`), сюда не берём.

    Невязка нормирована ОТДЕЛЬНО по каждому каналу (на масштаб самого канала):
    без этого `I_perp ~ 1`, `I_par ~ 1e-3`, и абсолютная невязка была почти
    целиком про t_perp — оптимизатор игнорировал форму t_par (замечено
    владельцем 2026-08-17: "t_par совершенно не сходится по характеру
    поведения", верный диагноз).

    ДВА ПРОХОДА с 3σ-отбраковкой ПО `I_perp` (владелец 2026-08-17: несколько
    точек уходят далеко за 3σ и тянут фит вниз), И ТОТ ЖЕ ШАБЛОН ЧАСТОТ
    применяется к `I_par` (владелец 2026-08-19). Частоты выбросов на `s2`
    (1.097, 1.113, 1.159, 1.166, 1.412 ТГц) совпадают со стандартным набором
    линий поглощения атмосферного водяного пара в непродуваемом ТГц-тракте --
    это ослабление всего импульса на этих частотах, общее для обоих каналов
    (свойство тракта, не поляризации), поэтому маска, найденная по чистому
    `I_perp` (где линии видны на фоне гладкой модели), переносится на `I_par`
    буквально по индексу бина, а не считается там заново -- отдельный 3σ по
    `I_par` даёт вместо этого широкий систематический уход почти на трети
    верхней полосы (неполнота 4-параметрической модели, не выброс) и портит
    отбраковку. σ оценивается робастно (через MAD, `1.4826*MAD` ~ sigma для
    нормального распределения) ИМЕННО ПОТОМУ, что обычное std само раздуто
    выбросами. Второй проход стартует из результата первого.
    """
    from scipy.optimize import least_squares

    scale_perp = float(np.max(np.abs(I_perp))) or 1.0
    scale_par = float(np.max(np.abs(I_par))) or 1.0
    lo, hi = np.array([2.0, 0.5, 0.0, 0.2]), np.array([120.0, 60.0, 3.0, 4.2])

    def resid(x, mask_perp, mask_par):
        P, D, loss, gamma = x
        model_perp, model_par = blanco_channel_model(freqs, P, D, loss, gamma=gamma)
        r_perp = (model_perp - I_perp) / scale_perp
        r_par = (model_par - I_par) / scale_par
        return np.concatenate([r_perp[mask_perp], r_par[mask_par]])

    all_perp = np.ones(len(I_perp), dtype=bool)
    all_par = np.ones(len(I_par), dtype=bool)
    res1 = least_squares(resid, list(p0), bounds=(lo, hi), args=(all_perp, all_par))

    def robust_mask(model, data):
        r = model - data
        mad = float(np.median(np.abs(r - np.median(r))))
        sigma = 1.4826 * mad if mad > 0 else (float(np.std(r)) or 1.0)
        bad = np.abs(r) > 3.0 * sigma
        # линия поглощения имеет ширину -- захватываем и соседний бин с каждой
        # стороны от выброса, не только сам пик (владелец 2026-08-19).
        bad_wide = bad | np.r_[bad[1:], False] | np.r_[False, bad[:-1]]
        return ~bad_wide

    model_perp1, _ = blanco_channel_model(freqs, *res1.x[:3], gamma=res1.x[3])
    mask_perp = robust_mask(model_perp1, I_perp)
    # Тот же шаблон бинов -- на I_par (поглощение атмосферы общее для обоих
    # каналов); отдельный 3σ по I_par вместо этого ловит широкий систематический
    # уход на верхней части полосы (неполнота модели, не выброс) -- не годится.
    mask_par = mask_perp.copy()

    res = least_squares(resid, list(res1.x), bounds=(lo, hi), args=(mask_perp, mask_par))
    res.mask_perp, res.mask_par = mask_perp, mask_par
    res.n_excluded_perp = int((~mask_perp).sum())
    res.n_excluded_par = int((~mask_par).sum())
    d_over_p = res.x[1] / res.x[0]
    res.at_bound = bool(np.any(np.isclose(res.x, lo, rtol=1e-3))
                        or np.any(np.isclose(res.x, hi, rtol=1e-3))
                        or d_over_p > 0.9)
    return res


def _expand_outlier_masks(ok: np.ndarray, blanco_res) -> tuple:
    """`blanco_res.mask_perp/mask_par` -- относительно `freqs_band[ok]`, короче
    полного `freqs_band`. Разворачивает обратно на полную длину (True = точка
    ИСКЛЮЧЕНА из фита -- не прошла `isfinite` или отбракована по 3σ) для
    отрисовки на графике; не-`ok` бины тоже считаются исключёнными."""
    outlier_perp = np.ones(ok.shape, dtype=bool)
    outlier_par = np.ones(ok.shape, dtype=bool)
    outlier_perp[ok] = ~blanco_res.mask_perp
    outlier_par[ok] = ~blanco_res.mask_par
    return outlier_perp, outlier_par


def _analyze_points(label: str, points: list, band) -> dict:
    """Общее ядро расчёта по списку `(angle, t_s, E_s, t_b, E_b)`.

    Углы могут повторяться (несколько точек на угол за разные дни/повторы) --
    это только улучшает обусловленность гармонического МНК, не мешает.
    """
    if not points:
        raise ValueError(f"нет точек для {label!r}")
    points = sorted(points, key=lambda p: p[0])
    angles = np.array([p[0] for p in points])

    T_time = np.array([time_transmission(*p[1:]) for p in points])
    warn_bad_t = angles[T_time > 1.0]

    fit_time = fit_angular(angles, T_time, WEIGHT_RELATIVE, order=4)
    freqs_ref, power_ref = reference_power_spectrum([(p[3], p[4]) for p in points])

    lo, hi = band
    spectra = [spectrum_ratio(*p[1:]) for p in points]
    freqs = spectra[0][0]
    mask = (freqs >= lo) & (freqs <= hi)
    freqs_band = freqs[mask]
    Y = np.array([T_nu[mask] for _, T_nu in spectra])

    spec_fit = fit_spectral(angles, Y, freqs_band, WEIGHT_RELATIVE, order=4)
    I_perp = spec_fit.ch["I_perp"]
    I_par = spec_fit.ch["I_par"]

    ok = np.isfinite(I_perp) & np.isfinite(I_par)
    blanco_res = fit_pure_blanco(freqs_band[ok], I_perp[ok], I_par[ok])
    P_fit, D_fit, loss_fit, gamma_fit = blanco_res.x
    outlier_perp, outlier_par = _expand_outlier_masks(ok, blanco_res)

    return dict(dataset=label, band=(lo, hi), angles=angles, T_time=T_time,
               warn_bad_t=warn_bad_t, fit_time=fit_time, freqs_band=freqs_band,
               I_perp=I_perp, I_par=I_par, spec_fit=spec_fit,
               outlier_perp=outlier_perp, outlier_par=outlier_par,
               freqs_ref=freqs_ref, power_ref=power_ref,
               blanco_res=blanco_res, P_fit=P_fit, D_fit=D_fit, loss_fit=loss_fit,
               gamma_fit=gamma_fit, theta0_fit=fit_time.params["theta0_deg_from_h2"])


def analyze(dataset: str, root: Path = DATA_ROOT, band=(0.3, 1.5)) -> dict:
    """Плоская конвенция `{dataset}_{angle}deg_rep{N}_{sig|bg}.txt` -> словарь для CLI/GUI.

    Бросает ValueError, если по датасету нет данных -- вызывающий код сам
    решает, как это показать (консоль или messagebox).
    """
    data = scan_series(Path(root), dataset)
    if not data:
        raise ValueError(f"нет данных для {dataset!r} в {root}")
    points = [(a, *vals) for a, vals in data.items()]
    return _analyze_points(dataset, points, band)


def analyze_road(folder: Path, band=(0.3, 1.5)) -> dict:
    """Дорожная конвенция (`NNN_aXX.txt`/`NNN_bg.txt`) -- через настоящий движок track_viewer.

    Не свой код: `Scan` сам разбирает дорогу и сам выбирает референс
    (ближайший/ранний/поздний, с явной пометкой на краях прогона), а
    `core.physics.compute_track` считает DC (предымпульс), ВЧ-шумовой пол
    (метод Нафтали) и, что здесь главное, **маскирует по динамическому
    диапазону каждый частотный бин отдельно для каждой точки** -- этого не
    было в моей первой версии, и без этого спектральные каналы на тёмных
    углах превращаются в шум, выданный за сигнал. Настройки -- дефолты
    владельца (`core/physics.py:Settings`, решения 2026-08-12/13), не свои.
    """
    scan = Scan(str(folder))
    if not scan.signals:
        raise ValueError(f"в {folder} не нашлось трасс дорожной конвенции NNN_aXX.txt/NNN_bg.txt")

    lo, hi = band
    # ref=later, не дефолт CLI (earlier): на границах дней (эта дорога снята за
    # 4 дня с разрывами) "earlier" утаскивает референс из ПРЕДЫДУЩЕГО дня и даёт
    # T=179% на 017_a0.txt (проверено: earlier/later/average сравнены явно,
    # только later убирает выброс -- см. чат 2026-08-17). Совпадает с тем, что
    # реально использовалось при ручном экспорте CSV в GUI (референс 018, не 016).
    settings = ph.Settings(band_lo=lo, band_hi=hi, ref_mode=ph.REF_LATER,
                           fit_on=True, fit_spectral_on=True,
                           fit_order=4, fit_weights=ph.FIT_RELATIVE)
    points = ph.compute_track(scan, settings)
    fit = fm.fit_track(points, settings)
    if not fit.ok:
        raise ValueError(f"фит не выполнен: {'; '.join(fit.notes) or 'см. точки'}")
    if fit.spectral is None:
        raise ValueError(f"побинный фит не выполнен: {'; '.join(fit.notes) or '?'}")

    angles, T_time = fit.time.theta, fit.time.U
    warn_bad_t = angles[T_time > 1.0]
    freqs_ref, power_ref = reference_power_spectrum([(tr.t, tr.E) for tr in scan.backgrounds])

    freqs_band = fit.spectral.freqs
    ch = fit.spectral.ch
    I_perp, I_par = ch["I_perp"], ch["I_par"]

    ok = np.isfinite(I_perp) & np.isfinite(I_par)
    blanco_res = fit_pure_blanco(freqs_band[ok], I_perp[ok], I_par[ok])
    P_fit, D_fit, loss_fit, gamma_fit = blanco_res.x
    outlier_perp, outlier_par = _expand_outlier_masks(ok, blanco_res)

    return dict(dataset=folder.name, band=(lo, hi), angles=angles, T_time=T_time,
               warn_bad_t=warn_bad_t, fit_time=fit.time, freqs_band=freqs_band,
               I_perp=I_perp, I_par=I_par, spec_fit=fit.spectral,
               outlier_perp=outlier_perp, outlier_par=outlier_par,
               freqs_ref=freqs_ref, power_ref=power_ref,
               blanco_res=blanco_res, P_fit=P_fit, D_fit=D_fit, loss_fit=loss_fit,
               gamma_fit=gamma_fit, theta0_fit=fit.time.params["theta0_deg_from_h2"])


def summary_text(r: dict) -> str:
    """Тот же текст, что печатает CLI, -- один источник и для консоли, и для GUI."""
    pr = r["fit_time"].params
    lo, hi = r["band"]
    lines = [f"{r['dataset']}: {len(r['angles'])} углов, "
             f"{r['angles'].min():+.0f}..{r['angles'].max():+.0f} град"]
    if r["warn_bad_t"].size:
        lines.append(f"[!] T > 100% на углах {list(r['warn_bad_t'])} -- проверить sig/bg")
    ft = r["fit_time"]
    hh_flag = "" if ft.hh_share <= ft.hh_chance else "  [!] выше случайного уровня"
    lines += ["", "фит Малюса (по времени, порядок 4, веса относительные):",
             f"  theta0 (офсет установки в ротаторе) = "
             f"{pr['theta0_deg_from_h2']:+.3f} +- {ft.theta0_err:.3f} град",
             f"  eta (утечка) = {pr['eta_amplitude']:.4f}",
             f"  экстинкция = {pr['extinction_db']:.2f} дБ",
             f"  гармоники 6-8 (перекос/биение сверх офсета): "
             f"{100*ft.hh_share:.1f}% при случайных {100*ft.hh_chance:.1f}%{hh_flag}"]
    lines += ["", f"спектральные каналы (модельно-независимо, {lo:.2f}-{hi:.2f} ТГц, "
             f"{r['spec_fit'].n_bins} бинов): |t_par|^2 < 0 в "
             f"{int(np.sum(r['I_par'] < 0))} бинах"]

    res = r["blanco_res"]
    P_fit, D_fit, loss_fit = r["P_fit"], r["D_fit"], r["loss_fit"]
    gamma_fit = r.get("gamma_fit", GAMMA_FIXED)
    lines += ["", "фит теории Бланко + общий множитель потерь (P, D, потери, gamma):",
             f"  P = {P_fit:.3f} мкм, D = {D_fit:.3f} мкм  (номинал паспорта: P=16, D=11)",
             f"  loss_factor = {loss_fit:.4f} дБ/ТГц^gamma, gamma = {gamma_fit:.3f}  "
             f"(номинал паспорта: 0.255 дБ/ТГц^1.58)",
             f"  невязка (cost) = {res.cost:.4g}, оптимизатор: {res.message}",
             f"  отбраковано по 3σ (MAD): I_perp {res.n_excluded_perp}, "
             f"I_par {res.n_excluded_par} бинов из {len(r['freqs_band'])}"]
    if res.at_bound:
        lines += [f"  [!] упёрся в границу поиска или D/P = {D_fit / P_fit:.3f} > 0.9 "
                 "-- не считать сходимостью, значения P/D/loss ниже НЕ измерение."]
    return "\n".join(lines)


def build_figure(r: dict):
    """matplotlib Figure с угловой и спектральной панелью -- общая для CLI и GUI."""
    import matplotlib.pyplot as plt

    angles, T_time, fit_time = r["angles"], r["T_time"], r["fit_time"]
    freqs_band, I_perp, I_par = r["freqs_band"], r["I_perp"], r["I_par"]
    P_fit, D_fit = r["P_fit"], r["D_fit"]
    loss_fit = r.get("loss_fit", 0.0)
    gamma_fit = r.get("gamma_fit", GAMMA_FIXED)
    theta0_used = r.get("theta0_manual", r.get("theta0_fit", 0.0))
    angular_absolute = r.get("angular_absolute", False)

    theta_grid = np.linspace(angles.min(), angles.max(), 400)
    T_fit_grid = fit_time.curve(theta_grid)
    I_perp_fit, I_par_fit = blanco_channel_model(freqs_band, P_fit, D_fit, loss_fit,
                                                 gamma=gamma_fit)
    # Угловая кривая -- на ШИРОКОЙ записанной полосе с весом |Ê_ref(nu)|^2, не
    # на узкой полосе фита каналов: согласовано с теоремой Парсеваля, которой
    # подчиняется T_time (владелец 2026-08-19, см. докстрок blanco_angular_curve).
    freqs_wide = r.get("freqs_ref", freqs_band)
    power_ref = r.get("power_ref")
    T_blanco_grid = blanco_angular_curve(theta_grid, freqs_wide, P_fit, D_fit, loss_fit,
                                         theta0_deg=theta0_used, gamma=gamma_fit,
                                         relative=not angular_absolute, weight=power_ref)

    fig, axes = plt.subplots(2, 2, figsize=(11, 7.5))
    (ax_a_lin, ax_a_db), (ax_s_lin, ax_s_db) = axes

    if angular_absolute:
        T_time_plot, T_fit_plot = T_time, T_fit_grid
        pts_label, fit_label = "измерено (абс., sig/bg)", "фит Малюса (по данным)"
        blanco_label = "Бланко, АБСОЛЮТНОЕ (Джонс, 2 идент. WGP, P,D,потери)"
        note = ("данные и фит Малюса -- абсолютное T=sig/bg (bg = референс БЕЗ\n"
                "аттенюатора); кривая Бланко -- полная Джонс-схема S1, оба WGP\n"
                "идентичны (P/D/потери/gamma -- параметры ОДНОГО WGP, не суммы)")
    else:
        # По умолчанию данные/фит Малюса -- НЕ пик theta0 модели, а пик самого
        # фита Малюса (владелец 2026-08-17): те же данные, тот же коэффициент
        # для точек и кривой, никакой отдельной нормировки на theta0 не нужно --
        # гармонический фит порядка 4 и так пикует практически в theta0.
        norm = float(np.max(T_fit_grid)) if T_fit_grid.size else 1.0
        norm = norm if norm > 0 else 1.0
        T_time_plot, T_fit_plot = T_time / norm, T_fit_grid / norm
        pts_label = "измерено (норм. на пик фита Малюса)"
        fit_label = "фит Малюса (норм., пик=1)"
        blanco_label = "Бланко относительное, норм. на θ0 (P,D,потери)"
        note = ("данные и фит Малюса нормированы на пик фита Малюса (один и тот\n"
                "же коэффициент для обеих кривых); кривая Бланко нормирована\n"
                "отдельно, на своё значение в θ0 -- пики близки, но не тождественны")

    for ax, y_pts, y_fit, y_blanco, ylabel in (
            (ax_a_lin, T_time_plot, T_fit_plot, T_blanco_grid, "T (лин.)"),
            (ax_a_db, db(T_time_plot), db(T_fit_plot), db(T_blanco_grid), "затухание, дБ")):
        ax.scatter(angles, y_pts, color="#2a78d6", zorder=3, label=pts_label)
        ax.plot(theta_grid, y_fit, color="#eb6834", lw=1.6, zorder=2,
               label=fit_label)
        ax.plot(theta_grid, y_blanco, color="#7a3fd6", lw=1.6, ls=(0, (4, 2)), zorder=2,
               label=blanco_label)
        ax.set_xlabel("угол первого WGP, град")
        ax.set_ylabel(ylabel)
        ax.legend(fontsize=7.5, loc="center")  # пик и края у кривой купола/ямы -- в центре только средние точки
        ax.set_title(f"угловая кривая ({r['dataset']})", fontsize=9)
        ax.text(0.02, 0.02, note, transform=ax.transAxes, fontsize=6.5, color="#888",
               ha="left", va="bottom")

    # точки, отбракованные по 3σ в fit_pure_blanco (обычно линии поглощения
    # водяного пара -- отдельная физика, не WGP) -- показываем крестиками, не
    # молча выкидываем.
    out_p = r.get("outlier_perp", np.zeros(len(freqs_band), dtype=bool))
    out_a = r.get("outlier_par", np.zeros(len(freqs_band), dtype=bool))
    keep_p, keep_a = ~out_p, ~out_a

    for ax, y_perp, y_par, y_perp_f, y_par_f, ylabel in (
            (ax_s_lin, I_perp, I_par, I_perp_fit, I_par_fit, "канал (лин.)"),
            (ax_s_db, db(I_perp), db(I_par), db(I_perp_fit), db(I_par_fit), "канал, дБ")):
        ax.scatter(freqs_band[keep_p], y_perp[keep_p], s=10, color="#2a78d6",
                  label="I_perp = |tp|^4 (данные)")
        ax.scatter(freqs_band[keep_a], y_par[keep_a], s=10, color="#1baf7a",
                  label="I_par = |tp|^2|ta|^2 (данные)")
        if out_p.any():
            ax.scatter(freqs_band[out_p], y_perp[out_p], s=28, marker="x", color="#999999",
                      label=f"I_perp искл. >3σ ({int(out_p.sum())}, водяной пар?)")
        if out_a.any():
            ax.scatter(freqs_band[out_a], y_par[out_a], s=28, marker="x", color="#555555",
                      label=f"I_par искл. >3σ ({int(out_a.sum())}, водяной пар?)")
        ax.plot(freqs_band, y_perp_f, color="#2a78d6", lw=1.6, ls="--", label="Бланко, 2 WGP (P,D)")
        ax.plot(freqs_band, y_par_f, color="#1baf7a", lw=1.6, ls="--")
        ax.set_xlabel("частота, ТГц")
        ax.set_ylabel(ylabel)
        ax.legend(fontsize=6.5)
        ax.set_title(f"каналы (Джонс, 2 идент. WGP): P={P_fit:.2f} D={D_fit:.2f} мкм, "
                     f"потери={loss_fit:.3f}, gamma={gamma_fit:.2f}", fontsize=8.5)

    fig.suptitle(f"{r['dataset']} -- прототип, без сверки с зоной A", fontsize=10)
    fig.tight_layout()
    return fig


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--dataset", default="att-11-16-s2",
                    help="плоская конвенция (игнорируется, если задан --road)")
    ap.add_argument("--root", default=str(DATA_ROOT))
    ap.add_argument("--road", default=None,
                    help="каталог дорожной конвенции (NNN_aXX.txt/NNN_bg.txt + CSV track_viewer)")
    ap.add_argument("--band", nargs=2, type=float, default=(0.3, 1.5),
                    help="полоса для спектральных каналов и фита Бланко, ТГц")
    ap.add_argument("--out", default=None, help="каталог для PNG (по умолчанию — рядом со скриптом)")
    args = ap.parse_args()

    try:
        if args.road:
            r = analyze_road(Path(args.road), tuple(args.band))
        else:
            r = analyze(args.dataset, Path(args.root), tuple(args.band))
    except ValueError as e:
        print(e)
        return 1
    print(summary_text(r))

    try:
        import matplotlib
        matplotlib.use("Agg")
    except Exception as e:                                          # noqa: BLE001
        print(f"\nmatplotlib недоступен ({e}) -- графики пропущены, только числа выше")
        return 0

    fig = build_figure(r)
    out_dir = Path(args.out) if args.out else HERE
    out_dir.mkdir(parents=True, exist_ok=True)
    out_path = out_dir / f"measured_curve_{r['dataset']}.png"
    fig.savefig(out_path, dpi=140)
    print(f"\nграфик: {out_path}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
