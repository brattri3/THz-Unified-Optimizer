#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""A25 — утечка `eta` как свойство ЭЛЕМЕНТА, а не тракта: измерение с ошибками.

ЧТО ПРОВЕРЯЕТСЯ И ПОЧЕМУ ЭТО НЕ РЫБАЛКА
---------------------------------------
Этап 1 (A22) даёт два комплексных спектра каналов, а вместе с ними —
модельно-независимое отношение

    r(nu) = T_par(nu) / T_perp(nu),   eta(nu) = |r(nu)|

Ключевое свойство `r`: она ИНВАРИАНТНА к любому общему комплексному множителю
поля. Всё, что тракт делает с лучом ПОЛЯРИЗАЦИОННО-НЕЗАВИСИМО — мощность
лазера, положение перетяжки, коллимация, окна, дрейф между днями, нормировка на
фон, скаляр `t_perp2(nu)` статичного элемента по теореме выравнивания T8 — это
общий множитель, и в `r` оно сокращается ТОЖДЕСТВЕННО (проверка 9 в
`field_linearization.selftest`).

Отсюда ПРЕДСКАЗАНИЕ, а не наблюдение: `eta(nu)` обязана быть характеристикой
самого поляризатора. Тракт может её сдвинуть только там, где он сам
поляризационно-зависим — перекошенный второй поляризатор (механизм выведен в
A24), неидеальный детектор как анализатор, смена ОСВЕЩАЕМОГО УЧАСТКА решётки.

Проверка построена как три контраста на одних и тех же данных:

  * тракт меняем, элемент тот же   -> `eta` обязана СОВПАСТЬ
        specac (перетяжка) / specac-2 (коллимированный пучок), разные дни;
        att-11-16-s1 (без плёнок) / -s2 / -s3 (обе плёнки TYDEX)
  * элемент меняем                 -> `eta` обязана РАЗОЙТИСЬ
        шесть элементов проекта
  * УЧАСТОК элемента меняем        -> положительный контроль
        purewave — единственный образец с известной нерегулярностью решётки,
        два дня светили в разные зоны

ЧЕСТНЫЕ ОШИБКИ
--------------
Статистическая ошибка регрессии здесь заведомо не главная: остатки степенного
закона СИЛЬНО скоррелированы по частоте (rho_1 = 0.30...0.84 на реальных
данных), а сам степенной закон — приближение. Поэтому бюджет ошибки собирается
из четырёх частей, и каждая измеряется:

  1. `stat`  — inv(X^T W X) регрессии;
  2. `birge` — умножение на chi^2/dof регрессии (несовершенство степенного
     закона), chi^2_red = 0.57...4.96;
  3. `boot`  — БЛОЧНЫЙ бутстрап остатков, длина блока из интегрального времени
     автокорреляции. Именно он ловит скоррелированность: наивный ВМНК при
     rho_1 = 0.5 занижает ошибку наклона примерно втрое (измерено в selftest);
  4. `jack`  — складной нож ПО УГЛАМ (leave-one-angle-out): ловит и качество
     углового плана, и одиночную плохую трассу;
  5. `syst`  — разброс по вариантам анализа: kappa, калибровка мультипликативной
     компоненты шума, порог SNR, снятие дрейфа, края полосы, политика фона.

Итоговая ковариация — сумма (2 или 3, что больше) + 4 + 5. Тест однородности
идёт по ней.

ЧТО ЗДЕСЬ НЕ ДЕЛАЕТСЯ (гардрейлы)
---------------------------------
* `data_pool/**`, `unified_optimizer/**` — только чтение; ядро зоны B не
  правится и не импортируется вовсе (этап 1 — чистая линейная алгебра).
* **В-21** замораживает ПЕРЕФИТЫ `purewave` и `test_grid_40_20` и пересчёт
  закона H5 до ответа зоны L по В-19. Здесь фитов модели нет: `eta` берётся
  линейным МНК в угловом базисе, ни одного параметра Бланко, ни `D_phys`, ни
  `D_eff/D_phys` не считается. Это, на мой взгляд, вне запрета — но допущение
  явное, и если владелец читает В-21 шире, `--skip-frozen` исключает оба
  датасета, а все выводы пересчитываются без них.
* `eta` на `att-11-16-*` — ЭФФЕКТИВНАЯ, не физическая: A24 показал, что тест
  ранга 1 там провален (lambda2/lambda1 до 2.3e-3 против 2.8e-5 на чистых
  образцах), а перекос неподвижного элемента ИМИТИРУЕТ утечку. Пометка
  `rank1_ok = False` идёт в каждую строку отчёта.

Запуск ИЗ КОРНЯ репозитория:
    .venv/Scripts/python.exe research/two_wgp/a25_eta_invariance.py --selftest
    .venv/Scripts/python.exe research/two_wgp/a25_eta_invariance.py --run
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
from scipy import stats as sps

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
for p in (str(REPO), str(HERE), str(REPO / "research" / "experiments")):
    if p not in sys.path:
        sys.path.insert(0, p)

import field_linearization as fl                                  # noqa: E402
import a22_field_vs_power as a22                                  # noqa: E402
import scan_loader                                                # noqa: E402

OUT = REPO / "research" / "results" / "two_wgp"

#: Общая полоса сравнения. `eta` растёт как ~nu, поэтому широкополосные числа с
#: РАЗНЫХ полос несравнимы напрямую; пересечение полос всех девяти прогонов —
#: [0.399, 1.358] ТГц, берём его с запасом внутрь.
COMMON_BAND = (0.40, 1.36)

#: Опорная частота, на которой сравниваются элементы. eta0 степенного закона —
#: это ровно eta(1 ТГц), поэтому пересчёта не требуется.
NU_REF_THZ = 1.0


class DS(object):
    """Паспорт прогона: какой ЭЛЕМЕНТ и какой ТРАКТ."""

    __slots__ = ("name", "element", "tract", "rank1_ok", "angles_limit",
                 "policy", "note", "frozen")

    def __init__(self, name, element, tract, rank1_ok=True, angles_limit=None,
                 policy=None, note="", frozen=False):
        self.name, self.element, self.tract = name, element, tract
        self.rank1_ok, self.angles_limit, self.policy = rank1_ok, angles_limit, policy
        self.note, self.frozen = note, frozen


#: `element` — физический поляризатор. Совпадение элемента при разном тракте и
#: есть проверяемая гипотеза; источники — data_pool/DATASETS.md и
#: data_pool/two_wgp_attenuator/README.md (конфигурации подтверждены владельцем
#: 2026-07-27), для specac — владелец 2026-08-16 («поляризатор извлекался из
#: ротатора», то есть образец тот же, схема другая).
DATASETS = [
    DS("specac", "SPECAC", "перетяжка, 13 углов 0..100",
       note="exp_data_1"),
    DS("specac-2", "SPECAC", "коллимированный пучок, 21 угол -70..100, повторы",
       note="exp_data_2, снят 15.08"),

    DS("att-11-16-s1", "ATT-A", "прибор без плёнок вовсе", rank1_ok=False,
       note="вращается WGP_A"),
    DS("att-11-16-s2", "ATT-A", "прибор + ОБЕ плёнки TYDEX", rank1_ok=False,
       note="вращается WGP_A"),
    DS("att-11-16-s3", "ATT-A", "прибор + обе плёнки, повтор съёмки", rank1_ok=False,
       note="вращается WGP_A, 12 углов"),
    DS("att-11-16-356", "ATT-B", "прибор + входная плёнка", rank1_ok=False,
       note="вращается WGP_B того же прибора"),

    DS("test_grid_33_11", "TG33", "одиночный WGP между плёнками, 37 углов"),
    DS("test_grid_40_20", "TG40", "одиночный WGP между плёнками, 25 углов",
       frozen=True, note="перефит заморожен В-21; здесь фита нет"),

    DS("purewave", "PW-all", "одиночный WGP, ОБА дня вместе", frozen=True,
       note="продакшн-диапазон смешивает дни — служебная строка"),
    DS("purewave@d1", "PW-d1", "день 28.07, углы -95..80", frozen=True,
       angles_limit=(-95.0, 80.0), note="восстановленный образец, зона A"),
    DS("purewave@d2", "PW-d2", "день 03.08, углы 82..110", frozen=True,
       angles_limit=(82.0, 110.0), note="тот же образец, ДРУГОЙ участок"),
]

#: Группы «один элемент — разный тракт». Только они входят в тест однородности.
TRACT_GROUPS = {"SPECAC": ["specac", "specac-2"],
                "ATT-A": ["att-11-16-s1", "att-11-16-s2", "att-11-16-s3"]}

#: Что считать РАЗНЫМИ элементами при оценке межэлементного разброса.
#: `PW-all` служебный (смешивает два участка), в разложение дисперсии не идёт.
ELEMENTS_FOR_SPREAD = ["SPECAC", "ATT-A", "ATT-B", "TG33", "TG40", "PW-d1", "PW-d2"]

#: Варианты анализа для систематики. Каждый — одна ручка, сдвинутая в разумную
#: сторону; систематика берётся как эмпирический разброс по ним.
#:
#: Третий элемент — идёт ли вариант В БЮДЖЕТ. `eps_mult=0` в бюджет НЕ идёт:
#: чисто аддитивная модель шума не «другой разумный выбор», а измеренно
#: неверная (chi^2/dof = 8.7 против 1.0, и независимо — разброс по повторам
#: specac-2 занижен ею в 42 раза, A22 §5). Включать её в ошибку значило бы
#: раздувать интервал заведомо забракованным вариантом, а раздутый интервал в
#: тесте СОГЛАСИЯ работает в пользу гипотезы — то есть в самую опасную сторону.
#: Число всё равно считается и кладётся в артефакт отдельной строкой: на
#: att-11-16-s1 оно даёт x1.74, и это само по себе диагноз.
VARIANTS = [
    ("kappa=1", {"kappa": 1.0}, True),
    ("eps_mult=0", {"eps_mult": 0.0}, False),
    ("snr>25", {"snr_min": 25.0}, True),
    ("drift", {"drift": True}, True),
    ("band_lo_cut", {"band_frac": (0.30, 0.0)}, True),
    ("band_hi_cut", {"band_frac": (0.0, 0.30)}, True),
]

_POLICY_VARIANTS = [("bg_preceding", {"policy": scan_loader.BG_PRECEDING}, True),
                    ("bg_nearest", {"policy": scan_loader.BG_NEAREST}, True)]


# ====================================================================== этап 1
def channels(name, kappa=fl.KAPPA_ABS_TO_VAR, eps_mult=None, drift=False,
             policy=None, angles_limit=None):
    """Этап 1 с полным набором ручек. -> словарь каналов и весов.

    Это те же вызовы `field_linearization`, что в `a23.stage1`; расширены только
    ручки (угловой диапазон, назначенный `eps_mult`), которые нужны для
    складного ножа и систематики. Совпадение с `a23.stage1` на умолчаниях
    проверяется в `selftest`, чтобы расширение не разъехалось с A23.
    """
    kw = {}
    if policy is not None:
        kw["policy"] = policy
    ang, freqs, exp, valid, meta, tnoise = a22.load_dataset(
        name, angles_limit=angles_limit, **kw)
    wm = meta["water_mask"]
    Y, f, tn = exp[:, wm], freqs[wm], tnoise[wm]

    if eps_mult is None:
        eps_mult, _ = fl.calibrate_noise(ang, Y, tn, kappa=kappa)
    s2 = fl.noise_variance(Y, tn, kappa=kappa, eps_mult=eps_mult)

    if drift:
        t0_pre = fl.theta0_rank1(fl.solve_bins(ang, Y, s2, order=2).z)["theta0_deg"]
        Y, _ = a22.fit_drift(ang, Y, s2, t0_pre)
        s2 = fl.noise_variance(Y, tn, kappa=kappa, eps_mult=eps_mult)

    prof = fl.theta0_profile(ang, Y, s2)
    theta0 = float(prof["theta0_deg"])
    bound = fl.solve_bins(ang, Y, s2, theta0_deg=theta0)
    Tp, Tl, Cov2 = fl.channels_from_z(bound)
    r, var_lnr = fl.ratio_and_var(Tp, Tl, Cov2)
    rank1 = fl.theta0_rank1(fl.solve_bins(ang, Y, s2, order=2).z)

    return {"angles": ang, "freqs": f, "Y": Y, "sigma2": s2, "tnoise": tn,
            "eps_mult": float(eps_mult), "theta0_deg": theta0,
            "sigma_theta0_deg": float(prof["sigma_theta0_deg"]),
            "T_perp": Tp, "T_par": Tl, "Cov2": Cov2, "ratio": r,
            "var_lnr": var_lnr, "bound": bound,
            "lambda2_over_lambda1": float(rank1["lambda2_over_lambda1"]),
            "cond_design": float(np.linalg.cond(
                fl.field_design(ang, 2).T @ fl.field_design(ang, 2))),
            "chi2_red_bound": float(np.mean(bound.chi2) / bound.dof)}


# =================================================================== регрессия
def _wls(x, y, w):
    """Прямая ВМНК y = p0 + p1*x. -> (p, C, rss, resid)."""
    X = np.column_stack([np.ones_like(x), x])
    N = X.T @ (w[:, None] * X)
    p = np.linalg.solve(N, X.T @ (w * y))
    C = np.linalg.inv(N)
    res = y - X @ p
    return p, C, float(np.sum(w * res ** 2)), res


def autocorr_time(res):
    """Интегральное время автокорреляции остатков (начальная положит. серия).

    Блочный бутстрап без него бессмыслен: длина блока обязана покрывать
    масштаб корреляции, иначе он вырождается в обычный бутстрап и повторяет
    ошибку наивного ВМНК.
    """
    r = np.asarray(res, dtype=np.float64)
    r = r - r.mean()
    n = len(r)
    if n < 8:
        return 1.0, 0.0
    var = float(np.dot(r, r) / n)
    if var <= 0:
        return 1.0, 0.0
    tau, rho1 = 1.0, 0.0
    for k in range(1, min(n // 2, 40)):
        rho = float(np.dot(r[:-k], r[k:]) / (n * var))
        if k == 1:
            rho1 = rho
        if rho <= 0.0:
            break
        tau += 2.0 * rho
    return float(max(tau, 1.0)), rho1


def block_bootstrap_cov(x, y, w, res, n_boot=600, block=None, rng=None):
    """Ковариация (p0, p1) круговым БЛОЧНЫМ бутстрапом остатков.

    Остатки переставляются блоками длины `block` и возвращаются на подогнанную
    прямую при НЕИЗМЕННОМ плане — так сохраняются и рычаги, и локальная
    корреляция. При block = 1 сводится к обычному бутстрапу; разница между ними
    и есть цена скоррелированности.
    """
    rng = rng or np.random.default_rng(20260817)
    n = len(x)
    if n < 6 or n_boot < 2:
        return np.zeros((2, 2)), int(block or 1)
    if block is None:
        tau, _ = autocorr_time(res)
        block = int(min(max(int(np.ceil(tau)), 1), max(n // 4, 1)))
    fitted = y - res
    n_blk = int(np.ceil(n / block))
    P = np.empty((n_boot, 2))
    for b in range(n_boot):
        st = rng.integers(0, n, size=n_blk)
        idx = ((st[:, None] + np.arange(block)[None, :]).ravel() % n)[:n]
        P[b] = _wls(x, fitted + res[idx], w)[0]
    return np.cov(P.T), int(block)


def powerlaw(freqs, r, var_lnr, valid, n_boot=600, rng=None):
    """ln|r| = ln(eta0) + a*ln(nu) со всеми оценками ошибки сразу.

    Возвращает z = (ln eta0, a) и ЧЕТЫРЕ ковариации: статистическую, с
    поправкой Бирджа, бутстраповскую и итоговую `cov_stat` = наибольшая из
    двух последних по следу. Кривизна проверяется отдельно: если квадратичный
    член значим, степенной закон неполон, и это надо знать до сравнения.
    """
    v = np.asarray(valid, dtype=bool)
    v = v & np.isfinite(var_lnr) & (var_lnr > 0) & (np.abs(r) > 0)
    out = {"n_used": int(v.sum())}
    if v.sum() < 6:
        out["ok"] = False
        return out
    x = np.log(np.asarray(freqs)[v])
    y = np.log(np.abs(np.asarray(r)[v]))
    w = 1.0 / np.asarray(var_lnr)[v]

    p, C, rss, res = _wls(x, y, w)
    dof = len(x) - 2
    chi2_red = rss / max(dof, 1)
    C_birge = C * max(chi2_red, 1.0)
    C_boot, block = block_bootstrap_cov(x, y, w, res, n_boot=n_boot, rng=rng)
    C_stat = C_boot if np.trace(C_boot) > np.trace(C_birge) else C_birge

    # --- кривизна: значим ли квадратичный член в ln nu
    #     Ошибка квадратичного члена умножается на sqrt(tau_int): при
    #     скоррелированных остатках независимых точек не n, а n/tau, и без
    #     этой поправки кривизна выходит «значимой» на 5-6 sigma просто
    #     потому, что соседние бины повторяют друг друга.
    tau, rho1 = autocorr_time(res)
    X2 = np.column_stack([np.ones_like(x), x, x ** 2])
    N2 = X2.T @ (w[:, None] * X2)
    p2 = np.linalg.solve(N2, X2.T @ (w * y))
    C2 = np.linalg.inv(N2)
    rss2 = float(np.sum(w * (y - X2 @ p2) ** 2))
    scale2 = max(rss2 / max(len(x) - 3, 1), 1.0)
    sig_q = float(np.sqrt(max(C2[2, 2], 0.0) * scale2 * tau))

    out.update({
        "ok": True, "z": p, "cov_wls": C, "cov_birge": C_birge,
        "cov_boot": C_boot, "cov_stat": C_stat,
        "eta0": float(np.exp(p[0])), "eta_exp": float(p[1]),
        "chi2_red": float(chi2_red), "block": block,
        "acf1": float(rho1), "tau_int": float(tau),
        "curvature": float(p2[2]),
        "curvature_sigma": float(sig_q),
        "curvature_nsig": float(abs(p2[2]) / sig_q) if sig_q > 0 else np.nan,
        "resid": res, "x": x, "y": y, "w": w,
    })
    return out


def resid_report(res, w):
    """Статистика остатков регрессии — гардрейл CLAUDE.md.

    Тренд по частоте здесь тождественно нулевой (он и подгоняется), поэтому
    структура читается по КРИВИЗНЕ и автокорреляции, а не по тренду.
    """
    z = res * np.sqrt(w)
    z = z[np.isfinite(z)]
    if len(z) < 4:
        return {"n": int(len(z))}
    tau, rho1 = autocorr_time(res)
    out = {"n": int(len(z)), "mean": float(np.mean(z)), "std": float(np.std(z)),
           "skew": float(sps.skew(z)), "kurtosis": float(sps.kurtosis(z)),
           "acf1": float(rho1), "tau_int": float(tau),
           "jarque_bera_p": float(sps.jarque_bera(z).pvalue)}
    try:
        out["shapiro_p"] = float(sps.shapiro(z).pvalue)
    except Exception:                                            # noqa: BLE001
        out["shapiro_p"] = float("nan")
    return out


# ================================================================== измерение
def _band_mask(freqs, band, band_frac=(0.0, 0.0)):
    lo, hi = band
    if band_frac != (0.0, 0.0):
        span = hi - lo
        lo, hi = lo + band_frac[0] * span, hi - band_frac[1] * span
    return (np.asarray(freqs) >= lo) & (np.asarray(freqs) <= hi), (lo, hi)


def _one(ds, band=COMMON_BAND, kappa=fl.KAPPA_ABS_TO_VAR, eps_mult=None,
         drift=False, policy=None, snr_min=9.0, band_frac=(0.0, 0.0),
         angles_limit="ds", n_boot=600, rng=None, ch=None):
    """Одно измерение (eta0, a) при заданном варианте анализа."""
    al = ds.angles_limit if angles_limit == "ds" else angles_limit
    if ch is None:
        ch = channels(ds.name.split("@")[0], kappa=kappa, eps_mult=eps_mult,
                      drift=drift, policy=policy or ds.policy, angles_limit=al)
    f, Cov2, Tl = ch["freqs"], ch["Cov2"], ch["T_par"]
    bmask, actual = _band_mask(f, band, band_frac)
    snr = np.abs(Tl) ** 2 / np.maximum(Cov2[:, 1, 1], 1e-300)
    valid = bmask & (snr > snr_min)
    pw = powerlaw(f, ch["ratio"], ch["var_lnr"], valid, n_boot=n_boot, rng=rng)
    pw["band"] = actual
    pw["n_band"] = int(bmask.sum())
    return pw, ch


def _eta_band(ch, band):
    """Широкополосная eta по ОБЩЕЙ полосе, со снятым шумовым смещением.

    E|T_hat|^2 = |T|^2 + Var, поэтому наивное отношение сумм квадратов завышает
    тёмный канал. Вычитаются диагонали Cov2 — это ровно снимаемое смещение.
    """
    m, _ = _band_mask(ch["freqs"], band)
    p2 = np.abs(ch["T_perp"][m]) ** 2 - ch["Cov2"][m, 0, 0]
    l2 = np.abs(ch["T_par"][m]) ** 2 - ch["Cov2"][m, 1, 1]
    naive = float(np.sqrt(np.sum(np.abs(ch["T_par"][m]) ** 2)
                          / np.sum(np.abs(ch["T_perp"][m]) ** 2)))
    deb = float(np.sqrt(max(np.sum(l2), 0.0) / max(np.sum(p2), 1e-300)))
    return deb, naive


def measure(ds, band=COMMON_BAND, n_boot=600, jackknife=True, systematics=True,
            rng=None, verbose=True):
    """Полное измерение одного прогона: z = (ln eta0, a) с бюджетом ошибки."""
    rng = rng or np.random.default_rng(abs(hash(ds.name)) % (2 ** 31))
    base, ch = _one(ds, band=band, n_boot=n_boot, rng=rng)
    if not base.get("ok"):
        return {"dataset": ds.name, "ok": False,
                "reason": "бинов в общей полосе %d — мало" % base.get("n_used", 0)}

    z0 = base["z"]
    C_stat = base["cov_stat"]

    # --- складной нож ПО УГЛАМ: качество углового плана и одиночная трасса
    C_jack = np.zeros((2, 2))
    jack_n = 0
    if jackknife:
        ang = ch["angles"]
        uniq = np.unique(ang)
        Z = []
        for a in uniq:
            keep = uniq[uniq != a]
            lim = (float(keep.min()) - 0.5, float(keep.max()) + 0.5)
            if len(keep) < 4:
                continue
            try:
                # угловой диапазон вырезать нельзя (он бы срезал край), поэтому
                # выкидываем угол вручную, пересобрав матрицу
                sub = _jack_channels(ch, ang != a)
                pw, _ = _one(ds, band=band, n_boot=0, ch=sub, rng=rng)
            except Exception:                                    # noqa: BLE001
                continue
            if pw.get("ok"):
                Z.append(pw["z"])
        if len(Z) >= 3:
            Z = np.array(Z)
            nj = len(Z)
            d = Z - Z.mean(axis=0)
            C_jack = (nj - 1.0) / nj * (d.T @ d)
            jack_n = nj

    # --- систематика: разброс по вариантам анализа
    C_syst = np.zeros((2, 2))
    var_rows = []
    if systematics:
        variants = list(VARIANTS)
        if ds.name in a22.DIR_DATASETS:
            variants += _POLICY_VARIANTS
        n_v = 0
        for label, kw, in_budget in variants:
            try:
                pw, _ = _one(ds, band=band, n_boot=0, rng=rng, **kw)
            except Exception as exc:                             # noqa: BLE001
                var_rows.append({"variant": label, "error": str(exc)})
                continue
            if not pw.get("ok"):
                continue
            d = pw["z"] - z0
            var_rows.append({"variant": label, "in_budget": bool(in_budget),
                             "d_ln_eta0": float(d[0]), "d_a": float(d[1]),
                             "eta0": pw["eta0"], "eta_exp": pw["eta_exp"]})
            if in_budget:
                C_syst = C_syst + np.outer(d, d)
                n_v += 1
        if n_v:
            C_syst = C_syst / n_v

    C_tot = C_stat + C_jack + C_syst
    C_tight = C_stat + C_jack               # без систематики — контроль теста
    eta_deb, eta_naive = _eta_band(ch, band)
    eps0 = next((v["d_ln_eta0"] for v in var_rows
                 if v.get("variant") == "eps_mult=0" and "d_ln_eta0" in v), None)

    res = {
        "dataset": ds.name, "element": ds.element, "tract": ds.tract,
        "rank1_ok": ds.rank1_ok, "frozen": ds.frozen, "note": ds.note,
        "ok": True,
        "n_angles": int(len(np.unique(ch["angles"]))),
        "n_traces": int(len(ch["angles"])),
        "cond_design": ch["cond_design"],
        "theta0_deg": ch["theta0_deg"],
        "sigma_theta0_deg": ch["sigma_theta0_deg"],
        "lambda2_over_lambda1": ch["lambda2_over_lambda1"],
        "eps_mult": ch["eps_mult"], "chi2_red_bound": ch["chi2_red_bound"],
        "band": base["band"], "n_bins_band": base["n_band"],
        "n_bins_used": base["n_used"],
        "eta0": base["eta0"], "eta_exp": base["eta_exp"],
        "z": z0,
        "cov_wls": base["cov_wls"], "cov_birge": base["cov_birge"],
        "cov_boot": base["cov_boot"], "cov_stat": C_stat,
        "cov_jack": C_jack, "cov_syst": C_syst, "cov_total": C_tot,
        "cov_tight": C_tight, "jack_n": jack_n,
        "sensitivity_additive_noise": (None if eps0 is None
                                       else float(np.exp(eps0))),
        "sigma_ln_eta0": {
            "wls": float(np.sqrt(max(base["cov_wls"][0, 0], 0))),
            "birge": float(np.sqrt(max(base["cov_birge"][0, 0], 0))),
            "boot": float(np.sqrt(max(base["cov_boot"][0, 0], 0))),
            "jack": float(np.sqrt(max(C_jack[0, 0], 0))),
            "syst": float(np.sqrt(max(C_syst[0, 0], 0))),
            "total": float(np.sqrt(max(C_tot[0, 0], 0)))},
        "sigma_a": {
            "wls": float(np.sqrt(max(base["cov_wls"][1, 1], 0))),
            "birge": float(np.sqrt(max(base["cov_birge"][1, 1], 0))),
            "boot": float(np.sqrt(max(base["cov_boot"][1, 1], 0))),
            "jack": float(np.sqrt(max(C_jack[1, 1], 0))),
            "syst": float(np.sqrt(max(C_syst[1, 1], 0))),
            "total": float(np.sqrt(max(C_tot[1, 1], 0)))},
        "chi2_red_powerlaw": base["chi2_red"],
        "block": base["block"], "acf1": base["acf1"], "tau_int": base["tau_int"],
        "curvature_nsig": base["curvature_nsig"],
        "eta_band_debiased": eta_deb, "eta_band_naive": eta_naive,
        "resid_stats": resid_report(base["resid"], base["w"]),
        "variants": var_rows,
        "_curve": {"freqs": ch["freqs"], "eta_nu": np.abs(ch["ratio"]),
                   "sigma_ln": np.sqrt(np.maximum(ch["var_lnr"], 0.0))},
    }
    if verbose:
        _print_row(res)
    return res


def _jack_channels(ch, keep):
    """Пересобирает каналы, выбросив часть трасс. Веса и tnoise — прежние."""
    ang, Y, s2 = ch["angles"][keep], ch["Y"][keep], ch["sigma2"][keep]
    prof = fl.theta0_profile(ang, Y, s2)
    t0 = float(prof["theta0_deg"])
    bound = fl.solve_bins(ang, Y, s2, theta0_deg=t0)
    Tp, Tl, Cov2 = fl.channels_from_z(bound)
    r, var = fl.ratio_and_var(Tp, Tl, Cov2)
    out = dict(ch)
    out.update({"angles": ang, "Y": Y, "sigma2": s2, "theta0_deg": t0,
                "T_perp": Tp, "T_par": Tl, "Cov2": Cov2, "ratio": r,
                "var_lnr": var, "bound": bound})
    return out


def _print_row(r):
    flag = "" if r["rank1_ok"] else "  [eta ЭФФЕКТИВНАЯ: ранг 1 провален]"
    print("  %-16s %-6s eta0 = %.5f x/: %.3f   a = %.3f +- %.3f   "
          "бинов %d  cond %7.1f%s"
          % (r["dataset"], r["element"], r["eta0"],
             float(np.exp(r["sigma_ln_eta0"]["total"])),
             r["eta_exp"], r["sigma_a"]["total"], r["n_bins_used"],
             r["cond_design"], flag))


# ====================================================================== тесты
def homogeneity(entries, key="cov_total"):
    """chi^2 равенства z = (ln eta0, a) внутри группы.

        chi^2 = sum_i (z_i - z_bar)^T C_i^-1 (z_i - z_bar),  dof = 2(N-1)

    z_bar — взвешенное среднее с C_bar = inv(sum inv(C_i)). Тест на ПОЛНОМ
    векторе, а не на одной eta0: ln eta0 и a скоррелированы (общий рычаг ln nu),
    и по одним диагоналям вывод был бы неверен.
    """
    zs = np.array([e["z"] for e in entries])
    Cs = [np.asarray(e[key], dtype=np.float64) for e in entries]
    N = len(entries)
    if N < 2:
        return {"n": N}
    Ws = [np.linalg.inv(C) for C in Cs]
    Wsum = sum(Ws)
    zbar = np.linalg.solve(Wsum, sum(W @ z for W, z in zip(Ws, zs)))
    chi2 = float(sum((z - zbar) @ W @ (z - zbar) for z, W in zip(zs, Ws)))
    dof = 2 * (N - 1)

    # скалярная версия — она читается прямо и идёт в текст
    v = np.array([C[0, 0] for C in Cs])
    x = zs[:, 0]
    mu = float(np.sum(x / v) / np.sum(1.0 / v))
    chi2_s = float(np.sum((x - mu) ** 2 / v))
    return {
        "n": N, "z_mean": zbar.tolist(),
        "chi2": chi2, "dof": dof, "p_value": float(sps.chi2.sf(chi2, dof)),
        "chi2_ln_eta0": chi2_s, "dof_ln_eta0": N - 1,
        "p_ln_eta0": float(sps.chi2.sf(chi2_s, N - 1)),
        "eta0_mean": float(np.exp(mu)),
        "sigma_eta0_mean_rel": float(np.sqrt(1.0 / np.sum(1.0 / v))),
        "max_spread_rel": float(np.exp(x.max() - x.min()) - 1.0),
        "members": [e.get("dataset", "?") for e in entries],
        "eta0": [float(np.exp(z[0])) for z in zs],
        # Широкополосная eta — НЕ инвариант: она взвешена спектром прибора
        # (сумма |T_perp|^2 по полосе), а он у разных трактов разный. Разброс
        # по ней кладётся рядом именно для контраста с eta0 = eta(1 ТГц).
        "eta_band": [e.get("eta_band_debiased") for e in entries],
        "eta_band_spread_rel": (
            float(max(b for b in bands) / min(b for b in bands) - 1.0)
            if (bands := [e.get("eta_band_debiased") for e in entries])
            and all(b for b in bands) else None),
    }


def variance_decomposition(results):
    """Разброс ln eta0 ВНУТРИ элемента против разброса МЕЖДУ элементами.

    Это и есть проверяемое утверждение в одном числе: если eta — свойство
    элемента, внутригрупповой разброс обязан быть заметно меньше межгруппового.
    Тест Фишера по ln eta0; `PW-d1` и `PW-d2` считаются РАЗНЫМИ элементами
    (разные участки нерегулярной решётки) — иначе положительный контроль
    попал бы во внутригрупповой разброс и размыл его.
    """
    by_el = {}
    for r in results:
        if not r.get("ok") or r["element"] not in ELEMENTS_FOR_SPREAD:
            continue
        by_el.setdefault(r["element"], []).append(r)

    within_ss, within_df, means = 0.0, 0, []
    per_el = []
    for el, rs in sorted(by_el.items()):
        x = np.array([r["z"][0] for r in rs])
        means.append(float(x.mean()))
        if len(x) > 1:
            within_ss += float(np.sum((x - x.mean()) ** 2))
            within_df += len(x) - 1
        per_el.append({"element": el, "n": len(x),
                       "eta0_mean": float(np.exp(x.mean())),
                       "spread_rel": float(np.exp(x.max() - x.min()) - 1.0)
                       if len(x) > 1 else None,
                       "members": [r["dataset"] for r in rs]})
    means = np.array(means)
    s_within = float(np.sqrt(within_ss / within_df)) if within_df else np.nan
    s_between = float(np.std(means, ddof=1)) if len(means) > 1 else np.nan
    out = {"per_element": per_el, "n_elements": int(len(means)),
           "s_within_ln": s_within, "s_between_ln": s_between,
           "within_rel_pct": 100.0 * (np.exp(s_within) - 1.0),
           "between_factor": float(np.exp(s_between)),
           "ratio": float(s_between / s_within) if s_within > 0 else np.nan,
           "within_df": within_df}
    if within_df and len(means) > 1:
        F = (s_between ** 2) / (s_within ** 2)
        out["F"] = float(F)
        out["F_dof"] = [len(means) - 1, within_df]
        out["F_p"] = float(sps.f.sf(F, len(means) - 1, within_df))
    return out


def conditioning_control(band=COMMON_BAND, n_boot=200, verbose=True):
    """Не фабрикует ли ПЛОХОЙ угловой план сдвиг `eta`?

    Это контроль под purewave: его второй день снят в окне 82..110 град, где
    cond(X^T X) = 2138 при гардрейле проекта 50. Если узкое окно на ЧИСТОМ
    образце воспроизводит полнодиапазонную `eta`, разницу между днями purewave
    нельзя списать на угловой план.
    """
    rows = []
    for name, lims in (("test_grid_33_11", [None, (80.0, 100.0), (-100.0, 80.0)]),
                       ("specac", [None, (60.0, 100.0)])):
        ref = None
        for lim in lims:
            ds = DS(name, "ctl", str(lim), angles_limit=lim)
            pw, ch = _one(ds, band=band, n_boot=n_boot)
            if not pw.get("ok"):
                continue
            if ref is None:
                ref = pw["eta0"]
            rows.append({"dataset": name, "angles_limit": lim,
                         "n_angles": int(len(np.unique(ch["angles"]))),
                         "cond": ch["cond_design"], "eta0": pw["eta0"],
                         "eta_exp": pw["eta_exp"],
                         "rel_to_full": float(pw["eta0"] / ref - 1.0)})
            if verbose:
                print("  %-16s %-16s углов %2d  cond %8.1f  eta0 %.5f  (%+.1f %%)"
                      % (name, str(lim), rows[-1]["n_angles"], rows[-1]["cond"],
                         rows[-1]["eta0"], 100 * rows[-1]["rel_to_full"]))
    return rows


# ======================================================================= прогон
def run(band=COMMON_BAND, n_boot=600, skip_frozen=False, fast=False,
        verbose=True):
    rng = np.random.default_rng(20260817)
    sel = [d for d in DATASETS if not (skip_frozen and d.frozen)]
    if verbose:
        print("=" * 86)
        print("A25 — eta как свойство элемента. Общая полоса %.2f-%.2f ТГц, "
              "опорная частота %.1f ТГц" % (band[0], band[1], NU_REF_THZ))
        print("=" * 86)
    results = []
    for ds in sel:
        r = measure(ds, band=band, n_boot=0 if fast else n_boot,
                    jackknife=not fast, systematics=not fast, rng=rng,
                    verbose=verbose)
        results.append(r)

    ok = [r for r in results if r.get("ok")]
    by_name = {r["dataset"]: r for r in ok}

    groups = {}
    for el, names in TRACT_GROUPS.items():
        ent = [by_name[n] for n in names if n in by_name]
        if len(ent) >= 2:
            g = homogeneity(ent)
            # Контроль: тот же тест на ошибке БЕЗ систематики. Согласие, которое
            # держится только на широком интервале, — не согласие.
            g["tight"] = {k: v for k, v in homogeneity(ent, key="cov_tight").items()
                          if k in ("chi2", "dof", "p_value", "chi2_ln_eta0",
                                   "dof_ln_eta0", "p_ln_eta0")}
            groups[el] = g

    vdec = variance_decomposition(ok)

    # положительный контроль: тот же тракт, ДРУГОЙ участок решётки
    pc = None
    if "purewave@d1" in by_name and "purewave@d2" in by_name:
        a, b = by_name["purewave@d1"], by_name["purewave@d2"]
        d = b["z"][0] - a["z"][0]
        sd = float(np.sqrt(a["cov_total"][0, 0] + b["cov_total"][0, 0]))
        pc = {"eta0_d1": a["eta0"], "eta0_d2": b["eta0"],
              "ratio": float(np.exp(d)), "d_ln": float(d),
              "sigma_ln": sd, "n_sigma": float(abs(d) / sd) if sd > 0 else np.nan,
              "cond_d1": a["cond_design"], "cond_d2": b["cond_design"]}

    ctl = conditioning_control(band=band, n_boot=0 if fast else 200,
                               verbose=verbose)

    out = {"band_thz": list(band), "nu_ref_thz": NU_REF_THZ,
           "n_datasets": len(ok), "datasets": results,
           "tract_groups": groups, "variance_decomposition": vdec,
           "patch_control_purewave": pc, "conditioning_control": ctl,
           "skip_frozen": bool(skip_frozen)}
    if verbose:
        _report(out)
    return out


def _report(out):
    print()
    print("-" * 86)
    print("ТЕСТ ОДНОРОДНОСТИ: один элемент, РАЗНЫЙ тракт")
    print("-" * 86)
    for el, g in out["tract_groups"].items():
        print("  %-7s n = %d (%s)" % (el, g["n"], ", ".join(g["members"])))
        print("        eta0 = %s" % ", ".join("%.5f" % v for v in g["eta0"]))
        print("        размах %.1f %%; среднее %.5f; chi^2 = %.2f при dof %d, "
              "p = %.3g" % (100 * g["max_spread_rel"], g["eta0_mean"],
                            g["chi2"], g["dof"], g["p_value"]))
        print("        только ln eta0: chi^2 = %.2f при dof %d, p = %.3g"
              % (g["chi2_ln_eta0"], g["dof_ln_eta0"], g["p_ln_eta0"]))
        t = g.get("tight", {})
        if t:
            print("        БЕЗ систематики в ошибке: chi^2 = %.2f при dof %d, "
                  "p = %.3g (только ln eta0: p = %.3g)"
                  % (t["chi2"], t["dof"], t["p_value"], t["p_ln_eta0"]))
        if g.get("eta_band_spread_rel") is not None:
            print("        для сравнения — ШИРОКОПОЛОСНАЯ eta: %s, размах %.1f %% "
                  "(она взвешена спектром прибора и инвариантом НЕ является)"
                  % (", ".join("%.5f" % v for v in g["eta_band"]),
                     100 * g["eta_band_spread_rel"]))

    v = out["variance_decomposition"]
    print("-" * 86)
    print("РАЗЛОЖЕНИЕ ДИСПЕРСИИ ln eta0 (элементов %d, внутригрупп. dof %d)"
          % (v["n_elements"], v["within_df"]))
    for pe in v["per_element"]:
        sp = ("размах %.1f %%" % (100 * pe["spread_rel"])
              if pe["spread_rel"] is not None else "один прогон")
        print("  %-7s n=%d  eta0 = %.5f  (%s)" % (pe["element"], pe["n"],
                                                  pe["eta0_mean"], sp))
    print("  внутри элемента: %.1f %%   между элементами: x%.2f   отношение %.1f"
          % (v["within_rel_pct"], v["between_factor"], v["ratio"]))
    if "F" in v:
        print("  F = %.1f при dof (%d, %d), p = %.3g"
              % (v["F"], v["F_dof"][0], v["F_dof"][1], v["F_p"]))

    pc = out.get("patch_control_purewave")
    if pc:
        print("-" * 86)
        print("ПОЛОЖИТЕЛЬНЫЙ КОНТРОЛЬ purewave — тот же тракт, ДРУГОЙ участок:")
        print("  eta0 день1 %.5f (cond %.1f) против день2 %.5f (cond %.1f): "
              "x%.2f, %.1f sigma" % (pc["eta0_d1"], pc["cond_d1"], pc["eta0_d2"],
                                     pc["cond_d2"], pc["ratio"], pc["n_sigma"]))


# ====================================================================== графики
def make_plots(out, results):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    OUT.mkdir(parents=True, exist_ok=True)
    ok = [r for r in results if r.get("ok")]
    els = sorted({r["element"] for r in ok})
    cmap = plt.get_cmap("tab10")
    col = {e: cmap(i % 10) for i, e in enumerate(els)}
    band = out["band_thz"]

    fig, ax = plt.subplots(2, 2, figsize=(13, 9))

    a = ax[0, 0]
    for r in ok:
        c = r["_curve"]
        m = (c["freqs"] >= band[0]) & (c["freqs"] <= band[1])
        ls = "-" if r["rank1_ok"] else "--"
        a.plot(c["freqs"][m], c["eta_nu"][m], ls, lw=1.1, color=col[r["element"]],
               alpha=.85, label="%s (%s)" % (r["dataset"], r["element"]))
    a.set_xlabel("частота, ТГц"); a.set_ylabel(r"$\eta=|T_\parallel/T_\perp|$")
    a.set_yscale("log"); a.grid(alpha=.3)
    a.set_title("утечка без модели: цвет = ЭЛЕМЕНТ, штрих = ранг 1 провален")
    a.legend(fontsize=6.5, ncol=2)

    a = ax[0, 1]
    ys = np.arange(len(ok))[::-1]
    for y, r in zip(ys, ok):
        s = r["sigma_ln_eta0"]["total"]
        a.errorbar(r["eta0"], y, xerr=[[r["eta0"] * (1 - np.exp(-s))],
                                       [r["eta0"] * (np.exp(s) - 1)]],
                   fmt="o", ms=5, color=col[r["element"]], capsize=3)
        a.text(r["eta0"], y + 0.28, r["dataset"], fontsize=6.5, ha="center")
    for el, g in out["tract_groups"].items():
        a.axvline(g["eta0_mean"], color=col.get(el, "k"), ls=":", lw=.9)
    a.set_yticks([]); a.set_xscale("log")
    a.set_xlabel(r"$\eta_0=\eta(1\,$ТГц$)$ с ПОЛНОЙ ошибкой")
    a.set_title("одинаковый цвет = один элемент")
    a.grid(alpha=.3, axis="x")

    a = ax[1, 0]
    for r in ok:
        rs = r["resid_stats"]
        a.scatter(r["acf1"], r["chi2_red_powerlaw"], s=40,
                  color=col[r["element"]])
        a.annotate(r["dataset"], (r["acf1"], r["chi2_red_powerlaw"]),
                   fontsize=6, xytext=(3, 3), textcoords="offset points")
    a.axhline(1.0, ls="--", c="k", lw=.8)
    a.set_xlabel(r"автокорреляция остатков $\rho_1$")
    a.set_ylabel(r"$\chi^2$/dof степенного закона")
    a.set_title("почему наивная ошибка не годится")
    a.grid(alpha=.3)

    a = ax[1, 1]
    labels = ["ВМНК", "Бирдж", "бутстрап", "нож по углам", "варианты", "ИТОГ"]
    keys = ["wls", "birge", "boot", "jack", "syst", "total"]
    w = 0.8 / max(len(ok), 1)
    for i, r in enumerate(ok):
        vals = [100 * (np.exp(r["sigma_ln_eta0"][k]) - 1) for k in keys]
        a.bar(np.arange(len(keys)) + i * w - 0.4, vals, width=w,
              color=col[r["element"]], alpha=.85)
    a.set_xticks(np.arange(len(keys))); a.set_xticklabels(labels, fontsize=8)
    a.set_ylabel(r"вклад в $\sigma(\eta_0)$, %")
    a.set_yscale("log"); a.grid(alpha=.3, axis="y")
    a.set_title("бюджет ошибки")

    fig.suptitle("A25 — утечка как свойство элемента: %.2f-%.2f ТГц"
                 % (band[0], band[1]))
    fig.tight_layout()
    p = OUT / "a25_eta_invariance.png"
    fig.savefig(p, dpi=115); plt.close(fig)
    return p


# =================================================================== selftest
def selftest(verbose=True):
    checks = []

    def chk(name, ok, detail):
        checks.append({"check": name, "ok": bool(ok), "detail": detail})
        if verbose:
            print(f"  [{'OK ' if ok else 'FAIL'}] {name}: {detail}")

    rng = np.random.default_rng(20260817)
    nu = np.linspace(0.4, 1.36, 90)
    eta0_t, a_t = 0.042, 1.05

    # --- 1. регрессия восстанавливает закон на бесшумном случае
    r = eta0_t * nu ** a_t * np.exp(1j * 0.3)
    pw = powerlaw(nu, r, np.full(len(nu), 1e-10), np.ones(len(nu), bool), n_boot=0)
    chk("степенной закон восстановлен точно",
        abs(pw["eta0"] - eta0_t) < 1e-10 and abs(pw["eta_exp"] - a_t) < 1e-10,
        "eta0 = %.12f, a = %.12f" % (pw["eta0"], pw["eta_exp"]))

    # --- 2. `_wls` совпадает с fl.logratio_regression (без дублирования логики)
    lr = fl.logratio_regression(nu, r, np.full(len(nu), 1e-10))
    chk("_wls совпадает с fl.logratio_regression",
        abs(lr["eta0"] - pw["eta0"]) < 1e-12 and abs(lr["eta_exp"] - pw["eta_exp"]) < 1e-12,
        "d(eta0) = %.1e, d(a) = %.1e" % (abs(lr["eta0"] - pw["eta0"]),
                                         abs(lr["eta_exp"] - pw["eta_exp"])))

    # --- 3. ln|r| несмещена при SNR >= 3 (нетривиально: смещения
    #        |T|^2 и логарифма гасят друг друга в первом порядке)
    th = np.arange(-90.0, 90.0, 10.0)
    Tp_t, Tl_t = 0.93 + 0.11j, 0.031 * np.exp(0.7j)
    E = fl._synth_field(th, Tp_t, Tl_t, 0.0)
    for snr in (3.0, 5.0):
        s = abs(Tl_t) / snr * np.sqrt(len(th) / 2.0)
        acc = []
        for _ in range(1500):
            n = (rng.normal(0, s, len(th)) + 1j * rng.normal(0, s, len(th))) / np.sqrt(2)
            s2 = np.full((len(th), 1), s ** 2)
            fb = fl.solve_bins(th, (E + n)[:, None], s2, theta0_deg=0.0)
            tp, tl, _ = fl.channels_from_z(fb)
            acc.append(np.log(abs(tl[0] / tp[0])))
        bias = float(np.mean(acc) - np.log(abs(Tl_t / Tp_t)))
        se = float(np.std(acc) / np.sqrt(len(acc)))
        chk("ln eta несмещена при SNR тени %.0f" % snr, abs(bias) < 4 * se + 0.01,
            "смещение %+.4f +- %.4f (то есть %+.1f %%)" % (bias, se,
                                                           100 * (np.exp(bias) - 1)))

    # --- 4. var_lnr откалибрована (Монте-Карло)
    s = abs(Tl_t) / 5.0 * np.sqrt(len(th) / 2.0)
    s2 = np.full((len(th), 1), s ** 2)
    acc = []
    for _ in range(2000):
        n = (rng.normal(0, s, len(th)) + 1j * rng.normal(0, s, len(th))) / np.sqrt(2)
        fb = fl.solve_bins(th, (E + n)[:, None], s2, theta0_deg=0.0)
        tp, tl, _ = fl.channels_from_z(fb)
        acc.append(np.log(abs(tl[0] / tp[0])))
    fb0 = fl.solve_bins(th, E[:, None], s2, theta0_deg=0.0)
    tp0, tl0, C0 = fl.channels_from_z(fb0)
    _, v0 = fl.ratio_and_var(tp0, tl0, C0)
    ratio = float(np.var(acc) / v0[0])
    chk("var(ln eta) откалибрована (МК)", abs(ratio - 1.0) < 0.15,
        "эмпирическая/модельная = %.3f" % ratio)

    # --- 5. ГЛАВНОЕ методическое: при скоррелированных остатках наивный ВМНК
    #        занижает ошибку, блочный бутстрап — нет
    n_rep, cov_naive, cov_boot = 300, 0, 0
    rho, sig = 0.6, 0.05
    x = np.log(nu)
    w = np.full(len(nu), 1.0 / sig ** 2)
    truth = np.array([np.log(eta0_t), a_t])
    for _ in range(n_rep):
        e = np.empty(len(nu))
        e[0] = rng.normal(0, sig)
        for i in range(1, len(nu)):
            e[i] = rho * e[i - 1] + rng.normal(0, sig * np.sqrt(1 - rho ** 2))
        y = truth[0] + truth[1] * x + e
        p, C, rss, res = _wls(x, y, w)
        Cb, _ = block_bootstrap_cov(x, y, w, res, n_boot=150, rng=rng)
        d = abs(p[1] - truth[1])
        cov_naive += d < 1.96 * np.sqrt(C[1, 1])
        cov_boot += d < 1.96 * np.sqrt(max(Cb[1, 1], 1e-300))
    pn, pb = cov_naive / n_rep, cov_boot / n_rep
    chk("блочный бутстрап чинит покрытие при rho=0.6 (наивный ВМНК — нет)",
        pn < 0.80 and pb > pn + 0.10,
        "покрытие 95 %%: наивно %.2f, блочно %.2f" % (pn, pb))

    # --- 6. длина блока растёт с корреляцией
    e1 = rng.normal(0, 1, 400)
    e2 = np.empty(400); e2[0] = 0
    for i in range(1, 400):
        e2[i] = 0.85 * e2[i - 1] + rng.normal(0, 1)
    t1, _ = autocorr_time(e1)
    t2, _ = autocorr_time(e2)
    chk("время автокорреляции различает белый шум и AR(1)",
        t1 < 2.0 and t2 > 5.0, "белый %.2f, AR(0.85) %.2f" % (t1, t2))

    # --- 7. chi^2 однородности имеет верное распределение под H0
    n_mc, C = 400, np.array([[4e-4, 1e-4], [1e-4, 9e-4]])
    L = np.linalg.cholesky(C)
    stats_ = []
    for _ in range(n_mc):
        ent = [{"z": truth + L @ rng.normal(size=2), "cov_total": C}
               for _ in range(3)]
        stats_.append(homogeneity(ent)["chi2"])
    m = float(np.mean(stats_))
    chk("chi^2 однородности: среднее = dof под H0", abs(m - 4.0) < 0.6,
        "среднее %.2f при dof 4" % m)

    # --- 8. тест ловит РАЗНЫЕ элементы
    ent = [{"z": truth + L @ rng.normal(size=2), "cov_total": C},
           {"z": truth + np.array([0.5, 0.0]) + L @ rng.normal(size=2),
            "cov_total": C}]
    h = homogeneity(ent)
    chk("тест отвергает равенство при разнице ln eta0 = 0.5",
        h["p_ln_eta0"] < 1e-6, "p = %.2e" % h["p_ln_eta0"])

    # --- 9. ОТРИЦАТЕЛЬНЫЙ КОНТРОЛЬ и одновременно основание гипотезы:
    #        общий комплексный множитель тракта не двигает ни eta0, ни a
    g = 0.61 * np.exp(-1.3j)
    Y0 = np.stack([fl._synth_field(th, Tp_t * (1 + 0.2 * f), Tl_t * (1 + 0.2 * f) * f,
                                   1.35) for f in nu], axis=1)
    def _meas(Y):
        s2_ = np.full(Y.shape, 1e-12)
        t0 = fl.theta0_profile(th, Y, s2_)["theta0_deg"]
        fb = fl.solve_bins(th, Y, s2_, theta0_deg=t0)
        tp, tl, C2 = fl.channels_from_z(fb)
        rr, vv = fl.ratio_and_var(tp, tl, C2)
        return powerlaw(nu, rr, np.maximum(vv, 1e-14), np.ones(len(nu), bool),
                        n_boot=0)
    m0, m1 = _meas(Y0), _meas(Y0 * g)
    chk("общий множитель тракта не двигает eta0 и a (основание гипотезы)",
        abs(m1["eta0"] / m0["eta0"] - 1) < 1e-10
        and abs(m1["eta_exp"] - m0["eta_exp"]) < 1e-10,
        "d(eta0)/eta0 = %.1e, d(a) = %.1e"
        % (abs(m1["eta0"] / m0["eta0"] - 1), abs(m1["eta_exp"] - m0["eta_exp"])))

    # --- 10. а поляризационно-ЗАВИСИМЫЙ тракт двигает — иначе тест пустой
    Y2 = np.stack([fl._synth_field(th, Tp_t * (1 + 0.2 * f),
                                   Tl_t * (1 + 0.2 * f) * f + 0.004 * Tp_t, 1.35)
                   for f in nu], axis=1)
    m2 = _meas(Y2)
    chk("поляризационно-зависимая добавка eta0 ДВИГАЕТ (тест не пустой)",
        abs(m2["eta0"] / m0["eta0"] - 1) > 0.05,
        "сдвиг eta0 на %+.1f %%" % (100 * (m2["eta0"] / m0["eta0"] - 1)))

    # --- 11. `channels` на умолчаниях совпадает с a23.stage1
    try:
        import a23_stage2_blanco as a23
        ch = channels("test_grid_33_11")
        s1 = a23.stage1("test_grid_33_11")
        d_t0 = abs(ch["theta0_deg"] - s1["theta0_deg"])
        d_eta = float(np.max(np.abs(np.abs(ch["ratio"])
                                    - np.abs(s1["T_par"] / s1["T_perp"]))))
        chk("channels() совпадает с a23.stage1 на умолчаниях",
            d_t0 < 1e-12 and d_eta < 1e-12,
            "d(theta0) = %.1e, max d(eta) = %.1e" % (d_t0, d_eta))
    except Exception as exc:                                     # noqa: BLE001
        chk("channels() совпадает с a23.stage1", False,
            "%s: %s" % (type(exc).__name__, exc))

    n_ok = sum(c["ok"] for c in checks)
    if verbose:
        print("\n  Итог: %d/%d проверок пройдено" % (n_ok, len(checks)))
    return {"n_ok": n_ok, "n_total": len(checks), "checks": checks}


# ======================================================================== main
def _jsonable(o):
    if isinstance(o, dict):
        return {k: _jsonable(v) for k, v in o.items() if not k.startswith("_")}
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
    ap = argparse.ArgumentParser(description="A25: eta как свойство элемента")
    ap.add_argument("--selftest", action="store_true")
    ap.add_argument("--run", action="store_true")
    ap.add_argument("--fast", action="store_true",
                    help="без ножа, бутстрапа и систематики — только точки")
    ap.add_argument("--skip-frozen", action="store_true",
                    help="исключить test_grid_40_20 и purewave (широкое чтение В-21)")
    ap.add_argument("--n-boot", type=int, default=600)
    ap.add_argument("--no-plots", action="store_true")
    args = ap.parse_args()

    if not (args.selftest or args.run):
        ap.print_help()
        return 0

    OUT.mkdir(parents=True, exist_ok=True)
    rc = 0
    if args.selftest:
        print("A25 eta-инвариантность — самопроверка")
        st = selftest(verbose=True)
        (OUT / "a25_eta_selftest.json").write_text(
            json.dumps(_jsonable(st), ensure_ascii=False, indent=1),
            encoding="utf-8")
        rc = 0 if st["n_ok"] == st["n_total"] else 1

    if args.run:
        out = run(n_boot=args.n_boot, skip_frozen=args.skip_frozen,
                  fast=args.fast, verbose=True)
        p = OUT / "a25_eta_invariance.json"
        p.write_text(json.dumps(_jsonable(out), ensure_ascii=False, indent=1),
                     encoding="utf-8")
        print("\nартефакт: %s" % p)
        if not args.no_plots:
            print("график:   %s" % make_plots(out, out["datasets"]))
    return rc


if __name__ == "__main__":
    raise SystemExit(main())
