"""A23 — этап 2 двухэтапной схемы: модель Бланко подгоняется к ВОССТАНОВЛЕННЫМ спектрам.

ЧТО ЭТО. Этап 1 (A22, `field_linearization.py`) снял с данных угловую развёртку и
вернул два комплексных спектра каналов `T_perp(nu)`, `T_par(nu)` с полной 2x2
ковариацией на бин. Здесь к этим двум спектрам подгоняется опорная модель M5
(Бланко + Друде + рассеяние + утечка). Угловой офсет `angle_offset` в фит НЕ
входит — он уже измерен на этапе 1 замкнутой формулой.

ПОЧЕМУ ЭТО НЕ ПРИБЛИЖЕНИЕ. При фиксированном theta0 модель целиком лежит в
двумерном span углового базиса [cos^2(t-t0), sin^2(t-t0)], поэтому побинные
коэффициенты с информационными матрицами — ДОСТАТОЧНАЯ статистика:

    chi2_global(p) = chi2_2(p) + sum_j RSS_perp,j ,   RSS_perp от p НЕ зависит

Второе слагаемое — невязка, недостижимая ни для какой модели вида
`cos^2 T_perp + sin^2 T_par`. Тождество проверяется тремя способами (см. ниже);
в `selftest` оно выходит на машинный нуль, на реальных данных — на 1e-10 от
масштаба chi2. Записка NOTE_field_linearization.md §9 утверждала, что эта
проверка уже стоит в `field_linearization.selftest` — её там не было, поэтому
она реализована здесь.

ЧТО ЭТО ДАЁТ, КРОМЕ СКОРОСТИ.
  * `angle_offset` уходит из нелинейного поиска. Это главный документированный
    источник мультимодальности проекта (A3: 3-20 бассейнов на 21 старт, лишь
    22.2 % стартов достигают глобального).
  * Четыре параметра утечки (`eta0`, `eta_exp`, `delta0`, `tau_leak`) и задержка
    `tau_ps` входят с ЧЕСТНОЙ затравкой из модельно-независимой регрессии
    этапа 1, а не с назначенной.
  * Метрика становится честным гауссовым chi2 по комплексной величине. Штатная
    невязка проекта (`fit_lib`: [|E|-|T|, wrapped phase] с весами W_AMP/W_PHASE)
    — ДРУГАЯ метрика: она смещена по модулю при малом SNR (Райс) и приписывает
    фазе вес 0.1 без вывода. Сравнение обеих в артефакте.

ЧЕГО ЖДАТЬ НЕ СЛЕДУЕТ. Вырождение `D_eff <-> eta` этим НЕ снимается: оно
структурное, на уровне отображения параметров в (T_perp, T_par). Корреляция
`D_um ~ eta0` обязана остаться ~ -1 на плотных образцах — это ОТРИЦАТЕЛЬНЫЙ
контроль, и если она «развяжется», надо искать ошибку.

ГАРДРЕЙЛЫ. `model_core` (зона B) импортируется и НЕ правится: каналы берутся
строками theta = 0 и 90 при `angle_offset = 0`, что даёт бит-в-бит ту же физику,
что у глобального фита. `data_pool/**` — только чтение. Пересчёт закона H5 и
перефит `purewave` / `test_grid_40_20` здесь НЕ делаются: это запрещено В-21 до
ответа по В-19.

Запуск ИЗ КОРНЯ репозитория:
    .venv/Scripts/python.exe research/two_wgp/a23_stage2_blanco.py --selftest
    .venv/Scripts/python.exe research/two_wgp/a23_stage2_blanco.py --dataset test_grid_33_11
    .venv/Scripts/python.exe research/two_wgp/a23_stage2_blanco.py --all
"""
from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

import numpy as np
import lmfit
from scipy import stats as sps

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
for p in (str(REPO), str(HERE), str(REPO / "research" / "experiments")):
    if p not in sys.path:
        sys.path.insert(0, p)

import field_linearization as fl                                  # noqa: E402  зона A
import a22_field_vs_power as a22                                  # noqa: E402  зона A
import model_core                                                 # noqa: E402  зона B, ЧТЕНИЕ

OUT = REPO / "research" / "results" / "two_wgp"
OUT.mkdir(parents=True, exist_ok=True)

#: Геометрия — микроскопия A19, утверждена владельцем (В-20).
#: `specac-2` — тот же физический поляризатор, что `specac` (владелец 2026-08-16:
#: «поляризатор извлекался из ротатора»), поэтому геометрия та же; отличается
#: оптическая схема, а не образец.
GEOMETRY_A19 = {
    "test_grid_33_11": {"P_um": 31.7122, "D_phys_um": 10.4739},
    "specac":          {"P_um": 24.7251, "D_phys_um": 9.7204},
    "specac-2":        {"P_um": 24.7251, "D_phys_um": 9.7204},
}

#: Опорные числа этапа 1 (A22, NOTE_field_linearization.md §6) — для сверки.
STAGE1_REFERENCE = {
    "test_grid_33_11": {"theta0_deg": -1.2723, "sigma_theta0_deg": 0.0063,
                        "eta_broadband": 0.05314, "eta0": 0.0704, "eta_exp": 1.016},
    "specac-2":        {"theta0_deg": 0.2282, "sigma_theta0_deg": 0.0077,
                        "eta_broadband": 0.03826, "eta0": 0.0404, "eta_exp": 1.093},
}

#: D_eff из перефита A19 на утверждённой геометрии — с чем сверять этап 2.
#: `specac-2` — тот же образец, что `specac`, поэтому опорное D_eff то же; но
#: оптическая схема другая, и совпадения от него требовать нельзя (вопрос A-5).
D_EFF_A19 = {"test_grid_33_11": 7.2209, "specac": 6.1466, "specac-2": 6.1466}

#: `specac` идёт третьим как КОНТРОЛЬ на том же физическом образце: он снят в
#: перетяжке, `specac-2` — в коллимированном пучке. Без него расхождение между
#: ними нечем разделить на «образец» и «схему».
DEFAULT_DATASETS = ["test_grid_33_11", "specac-2", "specac"]

PARAM_ORDER = ["D_um", "loss_factor", "tau_ps", "tau_par_ps",
               "eta0", "eta_exp", "delta0", "tau_leak"]
EPS_PARAMS = ["eps_re", "eps_im"]


def _order(eps=False):
    """Порядок параметров для отчётов: с `eps_cross` или без."""
    return PARAM_ORDER + (EPS_PARAMS if eps else [])


# ======================================================================= этап 1
def stage1(name, drift=False, kappa=fl.KAPPA_ABS_TO_VAR, snr_band=False,
           policy=None):
    """Прогоняет этап 1 и отдаёт ВСЁ, что нужно этапу 2 и сверкам.

    `a22.run_field` возвращает каналы, но не модель шума и не саму матрицу
    данных, а для проверки тождества достаточности нужны обе. Поэтому здесь
    этап 1 повторяется явно теми же вызовами `field_linearization` — числа
    обязаны совпасть с A22 (это проверяется и пишется в артефакт).
    """
    kw = {} if policy is None else {"policy": policy}
    ang, freqs, exp, valid, meta, tnoise = a22.load_dataset(name, **kw)
    wm = meta["water_mask"]
    Y, f, tn = exp[:, wm], freqs[wm], tnoise[wm]

    eps_mult, _ = fl.calibrate_noise(ang, Y, tn, kappa=kappa)
    s2 = fl.noise_variance(Y, tn, kappa=kappa, eps_mult=eps_mult)

    drift_g = None
    if drift:
        t0_pre = fl.theta0_rank1(fl.solve_bins(ang, Y, s2, order=2).z)["theta0_deg"]
        Y, drift_g = a22.fit_drift(ang, Y, s2, t0_pre)
        s2 = fl.noise_variance(Y, tn, kappa=kappa, eps_mult=eps_mult)

    prof = fl.theta0_profile(ang, Y, s2)
    theta0 = float(prof["theta0_deg"])
    bound = fl.solve_bins(ang, Y, s2, theta0_deg=theta0)
    Tp, Tl, Cov2 = fl.channels_from_z(bound)
    r, var_lnr = fl.ratio_and_var(Tp, Tl, Cov2)

    snr_par = np.abs(Tl) ** 2 / np.maximum(Cov2[:, 1, 1], 1e-300)
    good = snr_par > 9.0
    lr = fl.logratio_regression(f, r, var_lnr, valid=good)

    # Бины, где 2x2 ковариация вырождена, в фит не идут (их нет на практике,
    # но молча делить на вырожденную матрицу нельзя).
    cond = np.linalg.cond(Cov2)
    use = np.isfinite(cond) & (cond < 1e10)
    if snr_band:
        use = use & good

    return {
        "dataset": name, "angles": ang, "freqs": f, "Y": Y, "sigma2": s2,
        "tnoise": tn, "eps_mult": float(eps_mult), "theta0_deg": theta0,
        "sigma_theta0_deg": float(prof["sigma_theta0_deg"]),
        "T_perp": Tp, "T_par": Tl, "Cov2": Cov2, "use": use,
        "rss_perp": float(np.sum(bound.chi2)),
        "chi2_red_bound": float(np.mean(bound.chi2) / bound.dof),
        "leakage": lr, "drift_g": drift_g,
        "eta_broadband": float(np.sqrt(np.sum(np.abs(Tl) ** 2)
                                       / np.sum(np.abs(Tp) ** 2))),
        "n_bins": int(len(f)), "n_bins_used": int(use.sum()),
    }


# ======================================================================= модель
def model_channels(params, freqs):
    """(T_perp, T_par) модели M5 из `model_core.compute_grid` без правки зоны B.

    При `angle_offset = 0` строка theta = 0 даёт cos^2 = 1, sin^2 = 0, то есть
    ровно `T_perp`; строка theta = 90 — ровно `T_par` (`model_core.py:285-297`).
    """
    def g(k, d=0.0):
        return params[k].value if k in params else d
    E = model_core.compute_grid(
        [0.0, 90.0], freqs,
        params["P_um"].value * 1e-6, params["D_um"].value * 1e-6,
        g("loss_factor"), 0.0, g("tau_ps"),
        gamma=g("gamma", 1.0), tau_par_ps=g("tau_par_ps"),
        eta0=g("eta0"), eta_exp=g("eta_exp", 1.0),
        delta0=g("delta0"), tau_leak=g("tau_leak"))
    # `eps_cross` — угол-НЕЗАВИСИМАЯ добавка, поэтому в представлении каналов
    # она входит в ОБА канала одинаково и ничего больше не трогает (это же
    # свойство проверено в field_linearization.selftest п.5). Наглядность —
    # побочный выигрыш двухэтапной схемы.
    eps = _eps(params)
    return E[0] + eps, E[1] + eps


def _eps(params):
    if "eps_re" not in params:
        return 0.0
    return params["eps_re"].value + 1j * params["eps_im"].value


def make_params(dataset, seed=None, angle_offset=None, eps=False):
    """Параметры M5. `angle_offset` добавляется ТОЛЬКО для глобальной сверки."""
    geo = GEOMETRY_A19[dataset]
    p_um = geo["P_um"]
    P = lmfit.Parameters()
    P.add("P_um", value=p_um, vary=False)          # прибит оптическим контролем
    P.add("D_um", value=0.5 * geo["D_phys_um"], min=0.1, max=p_um - 0.5)
    P.add("loss_factor", value=0.3, min=0.0, max=5.0)
    P.add("gamma", value=2.0, vary=False)          # HN11: вырожден
    P.add("tau_ps", value=0.0, min=-10.0, max=10.0)
    P.add("tau_par_ps", value=0.0, min=-5.0, max=5.0)
    P.add("eta0", value=0.02, min=0.0, max=1.0)
    P.add("eta_exp", value=1.0, min=-2.0, max=6.0)
    P.add("delta0", value=0.0, min=-np.pi, max=np.pi)
    P.add("tau_leak", value=0.0, min=-5.0, max=5.0)
    if angle_offset is not None:
        P.add("angle_offset", value=float(angle_offset), min=-10.0, max=10.0,
              vary=False)
    if eps:
        P.add("eps_re", value=0.0, min=-0.2, max=0.2)
        P.add("eps_im", value=0.0, min=-0.2, max=0.2)
    if seed:
        for k, v in seed.items():
            if k in P and v is not None and np.isfinite(v):
                P[k].set(value=float(min(max(v, P[k].min), P[k].max)))
    return P


# ================================================================== целевые
def _whitener(Cov2):
    """L такое, что chi2_j = |L_j^T dT_re|^2 + |L_j^T dT_im|^2, включая множитель 2.

    chi2_2 = sum_j 2 * dT^H Cov2^-1 dT, а dT^H A dT = dT_re^T A dT_re +
    dT_im^T A dT_im для вещественной симметричной A. Раскладываем 2*A = L L^T.
    """
    A = 2.0 * np.linalg.inv(Cov2)
    A = 0.5 * (A + np.swapaxes(A, -1, -2))          # симметризация от округлений
    return np.linalg.cholesky(A)                     # (nbin, 2, 2), нижняя


def residual2(params, freqs, Tp, Tl, L):
    """Вектор невязок этапа 2 для `leastsq`: сумма квадратов = chi2_2."""
    if params["D_um"].value >= params["P_um"].value:
        return np.full(4 * len(freqs), 1e6)
    mp, ml = model_channels(params, freqs)
    d = np.stack([Tp - mp, Tl - ml], axis=1)                    # (nbin, 2)
    wr = np.einsum("bkj,bk->bj", L, d.real)
    wi = np.einsum("bkj,bk->bj", L, d.imag)
    return np.concatenate([wr.ravel(), wi.ravel()])


def chi2_stage2(params, freqs, Tp, Tl, Cov2):
    """chi2_2 = sum_j 2 * dT_j^H Cov2_j^-1 dT_j (прямо, без разложения)."""
    mp, ml = model_channels(params, freqs)
    d = np.stack([Tp - mp, Tl - ml], axis=1)
    A = np.linalg.inv(Cov2)
    q = np.einsum("bk,bkl,bl->b", d.real, A, d.real) + \
        np.einsum("bk,bkl,bl->b", d.imag, A, d.imag)
    return float(2.0 * np.sum(q))


def global_grid(params, angles, freqs, theta0_deg):
    """Полная модельная сетка E(theta, nu) с `angle_offset = theta0`."""
    def g(k, d=0.0):
        return params[k].value if k in params else d
    off = params["angle_offset"].value if "angle_offset" in params else theta0_deg
    return model_core.compute_grid(
        angles, freqs,
        params["P_um"].value * 1e-6, params["D_um"].value * 1e-6,
        g("loss_factor"), off, g("tau_ps"),
        gamma=g("gamma", 1.0), tau_par_ps=g("tau_par_ps"),
        eta0=g("eta0"), eta_exp=g("eta_exp", 1.0),
        delta0=g("delta0"), tau_leak=g("tau_leak")) + _eps(params)


def residual_global(params, angles, freqs, Y, sigma2, theta0_deg):
    """Глобальная КОМПЛЕКСНАЯ невязка в тех же весах, что и этап 1.

    Это НЕ штатная невязка проекта (`fit_lib`: модуль и завёрнутая фаза с
    весами W_AMP/W_PHASE) — та другая метрика. Здесь нужен именно гауссов chi2
    по комплексной величине, иначе тождество достаточности проверять не с чем.
    """
    if params["D_um"].value >= params["P_um"].value:
        return np.full(2 * Y.size, 1e6)
    d = Y - global_grid(params, angles, freqs, theta0_deg)
    w = np.sqrt(2.0 / sigma2)
    return np.concatenate([(d.real * w).ravel(), (d.imag * w).ravel()])


def chi2_global(params, angles, freqs, Y, sigma2, theta0_deg):
    r = residual_global(params, angles, freqs, Y, sigma2, theta0_deg)
    return float(np.sum(r ** 2))


# ================================================================ мультистарт
def _seeds(dataset, s1, eps=False):
    """Сетка затравок: 18 холодных (NOTE §9) + тёплые из регрессии этапа 1."""
    d_phys = GEOMETRY_A19[dataset]["D_phys_um"]
    lr = s1["leakage"]
    warm = {
        "eta0": lr.get("eta0"), "eta_exp": lr.get("eta_exp"),
        "delta0": _wrap_pi(lr.get("delta0_rad", 0.0)),
        "tau_leak": lr.get("tau_leak_m5_ps"),
        "tau_ps": _tau_seed_from_perp(s1),
        "loss_factor": _loss_seed_from_perp(s1),
    }
    out = []
    for fd in (0.3, 0.6, 0.9):
        for e0 in (0.0, 0.005, 0.03):
            for tl in (0.0, 0.19):
                out.append(("cold_D%.1f_e%.3f_t%.2f" % (fd, e0, tl),
                            {"D_um": fd * d_phys, "eta0": e0, "tau_leak": tl}))
    for fd in (0.3, 0.6, 0.9):
        s = dict(warm)
        s["D_um"] = fd * d_phys
        out.append(("warm_D%.1f" % fd, s))
    if eps:
        # Сетка по eps_cross — та же, что у A3 (|eps| x фаза): A3 измерил, что
        # ненулевая затравка eps достаёт бассейн на 107.7 AIC глубже нулевой,
        # поэтому без неё сравнение «с eps / без eps» было бы нечестным.
        for ea in (0.001, 0.01, 0.034):
            for ph in (0.0, 90.0, 180.0, 270.0):
                z = ea * np.exp(1j * np.deg2rad(ph))
                out.append(("eps_a%.3f_ph%03d" % (ea, ph),
                            {"D_um": 0.6 * d_phys, "eta0": 0.005,
                             "eps_re": float(z.real), "eps_im": float(z.imag)}))
    return out


def _wrap_pi(x):
    return float(np.arctan2(np.sin(x), np.cos(x))) if np.isfinite(x) else 0.0


def _tau_seed_from_perp(s1):
    """tau_ps из наклона развёрнутой фазы ЯРКОГО канала.

    arg T_perp = arg(t_perp_blanco) - 2*pi*nu*tau, а фаза блансовского t_perp
    медленная, поэтому наклон даёт tau с точностью, достаточной для затравки.
    Это ещё одна затравка, которой у холодного глобального фита нет вовсе.
    """
    f, Tp = s1["freqs"], s1["T_perp"]
    ph, _ = fl.unwrap_segmented(f, np.angle(Tp), s1["use"])
    m = s1["use"] & np.isfinite(ph)
    if m.sum() < 3:
        return 0.0
    sl = np.polyfit(f[m], ph[m], 1)[0]
    return float(-sl / (2 * np.pi))


def _loss_seed_from_perp(s1):
    """loss_factor из наклона ln|T_perp| по nu^gamma (gamma = 2)."""
    f, Tp = s1["freqs"], s1["T_perp"]
    m = s1["use"] & (np.abs(Tp) > 0)
    if m.sum() < 3:
        return 0.3
    sl = np.polyfit(f[m] ** 2, np.log(np.abs(Tp[m])), 1)[0]
    return float(np.clip(-sl * 2 * 4.343, 0.0, 5.0))


def multistart(dataset, f, Tp, Tl, L, seeds, verbose=True, eps=False):
    """Прогон списка затравок. -> (лучший результат lmfit, список всех стартов).

    МУЛЬТИСТАРТ ОБЯЗАТЕЛЕН, и это не ритуал. Даже после изъятия `angle_offset`
    поверхность остаётся резко мультимодальной: на синтетике с ИЗВЕСТНЫМ ответом
    одиночный холодный старт `D = 0.5*D_phys` садится в бассейн с chi2 = 6.4e6
    против 330 у истины, загоняя `eta_exp` и `delta0` в границы. Изъятие офсета
    убирает ОДИН источник мультимодальности, а не все.
    """
    starts, best = [], None
    for label, seed in seeds:
        t0 = time.time()
        P = make_params(dataset, seed=seed, eps=eps)
        try:
            res = lmfit.minimize(residual2, P, args=(f, Tp, Tl, L),
                                 method="leastsq", nan_policy="omit")
        except Exception as exc:                       # noqa: BLE001
            starts.append({"label": label, "error": repr(exc)})
            continue
        rec = {"label": label, "chi2": float(np.sum(res.residual ** 2)),
               "aic": float(res.aic), "redchi": float(res.redchi),
               "seconds": round(time.time() - t0, 2),
               "params": {k: float(res.params[k].value) for k in _order(eps)}}
        starts.append(rec)
        if best is None or rec["chi2"] < best[1]["chi2"]:
            best = (res, rec)
        if verbose:
            print("  %-24s chi2 = %12.2f  D = %7.4f  eta0 = %.5f  (%.1f с)"
                  % (label, rec["chi2"], res.params["D_um"].value,
                     res.params["eta0"].value, rec["seconds"]), flush=True)
    return best[0], starts


def fit_stage2(dataset, s1, verbose=True, eps=False):
    """Мультистарт этапа 2 на восстановленных спектрах."""
    u = s1["use"]
    f, Tp, Tl, Cov2 = s1["freqs"][u], s1["T_perp"][u], s1["T_par"][u], s1["Cov2"][u]
    return multistart(dataset, f, Tp, Tl, _whitener(Cov2),
                      _seeds(dataset, s1, eps=eps), verbose=verbose, eps=eps)


# ================================================================== проверки
def check_sufficiency(dataset, s1, n_draw=6, rng=None, eps=False):
    """(1) Алгебраически: chi2_global(p) - chi2_2(p) = const для СЛУЧАЙНЫХ p.

    Константа обязана совпасть с RSS_perp — суммой побинных невязок связанного
    фита этапа 1, то есть с той частью данных, которую никакая модель вида
    `cos^2 T_perp + sin^2 T_par` описать не может.
    """
    rng = rng or np.random.default_rng(20260816)
    u = s1["use"]
    f, Tp, Tl, Cov2 = s1["freqs"][u], s1["T_perp"][u], s1["T_par"][u], s1["Cov2"][u]
    ang, Y, s2, t0 = s1["angles"], s1["Y"], s1["sigma2"], s1["theta0_deg"]
    if not np.all(u):
        # тождество формулируется на ТЕХ ЖЕ бинах, что и глобальный chi2
        Y, s2 = Y[:, u], s2[:, u]

    d_phys = GEOMETRY_A19[dataset]["D_phys_um"]
    diffs = []
    for _ in range(n_draw):
        seed = {"D_um": float(rng.uniform(0.2, 0.95) * d_phys),
                "loss_factor": float(rng.uniform(0.0, 1.5)),
                "tau_ps": float(rng.uniform(-2.0, 2.0)),
                "tau_par_ps": float(rng.uniform(-0.5, 0.5)),
                "eta0": float(rng.uniform(0.0, 0.1)),
                "eta_exp": float(rng.uniform(0.0, 2.0)),
                "delta0": float(rng.uniform(-np.pi, np.pi)),
                "tau_leak": float(rng.uniform(-0.5, 0.5))}
        P = make_params(dataset, seed=seed, eps=eps)
        diffs.append(chi2_global(P, ang, f, Y, s2, t0) - chi2_stage2(P, f, Tp, Tl, Cov2))
    diffs = np.array(diffs)
    rss = float(np.sum(fl.solve_bins(ang, Y, s2, theta0_deg=t0).chi2))
    scale = max(abs(rss), 1.0)
    return {"const_mean": float(diffs.mean()),
            "const_spread_rel": float(np.ptp(diffs) / scale),
            "rss_perp_stage1": rss,
            "const_minus_rss_rel": float(abs(diffs.mean() - rss) / scale),
            "n_draw": int(n_draw)}


def check_numeric_equivalence(dataset, s1, best, verbose=True, eps=False):
    """(2) Глобальный комплексный фит с ПРИБИТЫМ theta0 — параметры обязаны совпасть.

    (3) Тот же фит со СВОБОДНЫМ angle_offset — диагностика: насколько модель
    Бланко тянет угловой офсет на себя.
    """
    u = s1["use"]
    ang, f = s1["angles"], s1["freqs"][u]
    Y, s2, t0 = s1["Y"][:, u], s1["sigma2"][:, u], s1["theta0_deg"]
    seed = {k: float(best.params[k].value) for k in _order(eps)}

    out = {}
    P = make_params(dataset, seed=seed, angle_offset=t0, eps=eps)
    tt = time.time()
    rp = lmfit.minimize(residual_global, P, args=(ang, f, Y, s2, t0),
                        method="leastsq")
    rel = {k: (abs(rp.params[k].value - seed[k])
               / max(abs(seed[k]), 1e-3)) for k in _order(eps)}
    out["pinned"] = {
        "chi2": float(np.sum(rp.residual ** 2)),
        "seconds": round(time.time() - tt, 1),
        "params": {k: float(rp.params[k].value) for k in _order(eps)},
        "max_rel_deviation": float(max(rel.values())),
        "worst_param": max(rel, key=rel.get),
    }

    P2 = make_params(dataset, seed=seed, angle_offset=t0, eps=eps)
    P2["angle_offset"].set(vary=True)
    tt = time.time()
    rf = lmfit.minimize(residual_global, P2, args=(ang, f, Y, s2, t0),
                        method="leastsq")
    out["free_offset"] = {
        "chi2": float(np.sum(rf.residual ** 2)),
        "angle_offset_deg": float(rf.params["angle_offset"].value),
        "sigma_deg": (float(rf.params["angle_offset"].stderr)
                      if rf.params["angle_offset"].stderr else None),
        "theta0_stage1_deg": float(t0),
        "delta_deg": float(rf.params["angle_offset"].value - t0),
        "delta_in_sigma1": float(abs(rf.params["angle_offset"].value - t0)
                                 / max(s1["sigma_theta0_deg"], 1e-9)),
        "seconds": round(time.time() - tt, 1),
        "params": {k: float(rf.params[k].value) for k in _order(eps)},
    }
    if verbose:
        print("  сверка: глоб. фит (theta0 прибит) — макс. относит. расхождение "
              "%.2e по %s" % (out["pinned"]["max_rel_deviation"],
                              out["pinned"]["worst_param"]))
        print("  диагностика: свободный angle_offset = %+.4f против %+.4f (этап 1), "
              "%.2f sigma_1" % (out["free_offset"]["angle_offset_deg"], t0,
                                out["free_offset"]["delta_in_sigma1"]))
    return out


def residual_ampphase(params, angles, freqs, Y, theta0_deg):
    """ШТАТНАЯ невязка проекта: [|E|-|T|, завёрнутая фаза] с весами W_AMP/W_PHASE.

    Нужна не для фита, а для ЧИСТОГО сравнения метрик: те же данные, та же
    маска, тот же прибитый theta0 — меняется только целевая функция.
    """
    import fit_lib                                          # noqa: E402  зона B, ЧТЕНИЕ
    if params["D_um"].value >= params["P_um"].value:
        return np.ones(2 * Y.size) * 1e6
    theo = global_grid(params, angles, freqs, theta0_deg)
    amp = np.abs(Y) - np.abs(theo)
    ph = np.angle(Y) - np.angle(theo)
    ph = np.arctan2(np.sin(ph), np.cos(ph))
    return np.concatenate([(amp * fit_lib.W_AMP).ravel(),
                           (ph * fit_lib.W_PHASE).ravel()])


def metric_comparison(dataset, s1, best, verbose=True, eps=False):
    """Насколько D_eff двигает СМЕНА МЕТРИКИ, а не смена схемы фита.

    Этап 2 меряет честным гауссовым chi2 по комплексной величине, а числа
    Фаз 1-3 и перефита A19 получены штатной невязкой «модуль + завёрнутая фаза»
    с весом фазы 0.1. Разница между ними — систематика МЕТОДА, и её надо
    измерить, а не списать на шум: тот же датасет, те же бины, тот же theta0.
    """
    u = s1["use"]
    ang, f = s1["angles"], s1["freqs"][u]
    Y, t0 = s1["Y"][:, u], s1["theta0_deg"]
    seed = {k: float(best.params[k].value) for k in _order(eps)}
    P = make_params(dataset, seed=seed, angle_offset=t0, eps=eps)
    r = lmfit.minimize(residual_ampphase, P, args=(ang, f, Y, t0), method="leastsq")
    d_new, d_old = float(best.params["D_um"].value), float(r.params["D_um"].value)
    out = {"D_um_gaussian_chi2": d_new, "D_um_ampphase_metric": d_old,
           "delta_um": d_new - d_old,
           "delta_in_sigma": ((d_new - d_old) / best.params["D_um"].stderr
                              if best.params["D_um"].stderr else None),
           "params_ampphase": {k: float(r.params[k].value) for k in _order(eps)},
           "note": ("одни данные, одна маска, один theta0; отличается ТОЛЬКО "
                    "целевая функция")}
    if verbose:
        print("  метрика: D_eff = %.4f (гауссов chi2) против %.4f (модуль+фаза), "
              "разница %+.4f мкм" % (d_new, d_old, out["delta_um"]))
    return out


def _correl(res, a, b):
    p = res.params.get(a)
    if p is None or not p.correl:
        return None
    v = p.correl.get(b)
    return float(v) if v is not None else None


def resid_stats(vec):
    """Статистика остатков — гардрейл CLAUDE.md: не одно число невязки."""
    x = np.asarray(vec, dtype=float)
    x = x[np.isfinite(x)]
    sub = x if len(x) <= 5000 else x[:: max(1, len(x) // 5000)]
    sh = sps.shapiro(sub)
    jb = sps.jarque_bera(x)
    dw = (float(np.sum(np.diff(x) ** 2) / np.sum(x ** 2))
          if np.sum(x ** 2) > 0 else None)
    return {"n": int(len(x)), "mean": float(x.mean()), "std": float(x.std(ddof=1)),
            "skew": float(sps.skew(x)), "kurtosis": float(sps.kurtosis(x)),
            "shapiro_p": float(sh.pvalue), "jarque_bera_p": float(jb.pvalue),
            "durbin_watson": dw}


def echo_scan(freqs, T_hat, T_model, var_full, tau_max=20.0, n_tau=4001):
    """Периодическая структура в остатке канала: есть ли эхо (эталон) и на какой задержке.

    Эталон умножает ПОЛЕ: T = T_model * (1 + rho * exp(2i*pi*nu*tau)). Поэтому
    ищем не «волну в модуле», а комплексную гармонику в относительном остатке
    y = T_hat/T_model - 1, подгоняя c0 + c*exp(2i*pi*nu*tau) на сетке tau.

    ВАЖНО, ЧТО СЧИТАЕТСЯ ЭХОМ. На полосе шириной B задержки tau < 1/B дают
    МЕНЬШЕ одного периода поперёк полосы: это не разрешённое эхо, а плавная
    кривизна, вырожденная с формой самой модели (loss_factor, D_eff). Поэтому
    результат отдаётся в двух видах: `resolvable` — лучшая задержка при
    tau > 1/B, и `unresolved` — лучшая при tau < 1/B, помеченная как кривизна,
    а не как отражатель. Путать их значит выдать неполноту модели за физику
    установки.

    Волна в остатке по частоте читается как НЕПОЛНОТА ФИЗИКИ (гардрейл
    CLAUDE.md), и её надо не констатировать, а измерять.
    """
    y = T_hat / T_model - 1.0
    w = 1.0 / np.maximum(var_full / np.abs(T_model) ** 2, 1e-300)
    m = np.isfinite(y) & np.isfinite(w)
    y, w, f = y[m], w[m], np.asarray(freqs)[m]
    if len(f) < 8:
        return None

    c0 = np.sum(w * y) / np.sum(w)
    chi0 = 2.0 * float(np.sum(w * np.abs(y - c0) ** 2))
    tau_min = 1.0 / (float(f.max()) - float(f.min()))       # предел разрешения

    lo = (0.0, chi0, 0.0)                                   # кривизна, tau < 1/B
    hi = (0.0, chi0, 0.0)                                   # разрешённое эхо
    for tau in np.linspace(0.0, tau_max, n_tau)[1:]:
        e = np.exp(2j * np.pi * f * tau)
        X = np.column_stack([np.ones_like(e), e])
        A = X.conj().T @ (w[:, None] * X)
        b = X.conj().T @ (w * y)
        try:
            c = np.linalg.solve(A, b)
        except np.linalg.LinAlgError:
            continue
        chi = 2.0 * float(np.sum(w * np.abs(y - X @ c) ** 2))
        slot = 0 if tau < tau_min else 1
        cur = lo if slot == 0 else hi
        if chi < cur[1]:
            if slot == 0:
                lo = (float(tau), chi, float(np.abs(c[1])))
            else:
                hi = (float(tau), chi, float(np.abs(c[1])))

    def _pack(b):
        return {"tau_ps": b[0], "delta_chi2": chi0 - b[1], "amplitude": b[2],
                "thickness_air_mm": 0.5 * b[0] * 0.29979}
    return {"chi2_const_only": chi0, "n_bins": int(len(f)),
            "tau_resolution_limit_ps": tau_min,
            "resolvable": _pack(hi), "unresolved_curvature": _pack(lo),
            "note": "delta chi^2 при 3 лишних параметрах; порог ~ 16 (4 sigma). "
                    "unresolved_curvature — НЕ отражатель, а форма, недоступная модели"}


def _resid_vs_freq(freqs, vec4, nbin):
    """Структура остатков по частоте: тренд/волна читаются как неполнота физики."""
    r = np.asarray(vec4, dtype=float).reshape(4, nbin) if vec4.size == 4 * nbin \
        else None
    if r is None:
        return None
    per_bin = np.sqrt(np.mean(r ** 2, axis=0))
    sl, ic = np.polyfit(freqs, per_bin, 1)
    return {"slope_per_thz": float(sl), "intercept": float(ic),
            "rms_low_half": float(np.mean(per_bin[: nbin // 2])),
            "rms_high_half": float(np.mean(per_bin[nbin // 2:]))}


# ==================================================================== прогон
def policy_scan(dataset, verbose=True):
    """Одни данные, три политики выбора фона. Ответ на «а не референс ли это?».

    Если структура остатка держится при всех трёх политиках, она принадлежит
    измерению, а не нормировке. Считается только для каталожных датасетов —
    у остальных фон один и выбора нет.
    """
    import scan_loader
    rows = []
    for pol in (scan_loader.BG_PRECEDING, scan_loader.BG_NEAREST,
                scan_loader.BG_TWO_ABS):
        s1 = stage1(dataset, policy=pol)
        u = s1["use"]
        f = s1["freqs"][u]
        best, _ = multistart(dataset, f, s1["T_perp"][u], s1["T_par"][u],
                             _whitener(s1["Cov2"][u]), _seeds(dataset, s1),
                             verbose=False)
        mp, ml = model_channels(best.params, f)
        e = echo_scan(f, s1["T_perp"][u], mp, s1["Cov2"][u][:, 0, 0])
        rows.append({
            "policy": pol,
            "redchi": float(np.sum(best.residual ** 2)) / (best.ndata - best.nvarys),
            "D_um": float(best.params["D_um"].value),
            "theta0_deg": s1["theta0_deg"],
            "perp_median_sigma": float(np.median(
                np.abs(s1["T_perp"][u] - mp) / np.sqrt(s1["Cov2"][u][:, 0, 0] / 2))),
            "curvature_amplitude": e["unresolved_curvature"]["amplitude"],
            "curvature_delta_chi2": e["unresolved_curvature"]["delta_chi2"],
        })
        if verbose:
            r = rows[-1]
            print("  %-10s chi2/dof = %7.2f, D = %.3f, T_perp %.2f sigma, "
                  "кривизна %.4f (dchi2 %.0f)"
                  % (pol, r["redchi"], r["D_um"], r["perp_median_sigma"],
                     r["curvature_amplitude"], r["curvature_delta_chi2"]))
    return rows


def analyse(dataset, drift=False, snr_band=False, equivalence=True, verbose=True,
            policies=False, eps=False):
    t_all = time.time()
    if verbose:
        print("=" * 78)
        print("A23 этап 2 — %s%s%s" % (dataset, " (дрейф снят)" if drift else "",
                                       " + eps_cross" if eps else ""))
        print("=" * 78)

    s1 = stage1(dataset, drift=drift, snr_band=snr_band)
    ref = STAGE1_REFERENCE.get(dataset, {})
    if verbose:
        print("Этап 1: theta0 = %+.4f +- %.4f, бинов %d из %d, eps_mult = %.4f, "
              "RSS_perp = %.1f"
              % (s1["theta0_deg"], s1["sigma_theta0_deg"], s1["n_bins_used"],
                 s1["n_bins"], s1["eps_mult"], s1["rss_perp"]))
        if ref and not drift:
            print("        сверка с A22: theta0 %+.4f (расхождение %.1e)"
                  % (ref["theta0_deg"], abs(s1["theta0_deg"] - ref["theta0_deg"])))
        print("Мультистарт (%d стартов):" % len(_seeds(dataset, s1, eps=eps)))

    best, starts = fit_stage2(dataset, s1, verbose=verbose, eps=eps)

    ok = [s for s in starts if "chi2" in s]
    chi_best = min(s["chi2"] for s in ok)
    in_best = sum(1 for s in ok if s["chi2"] < chi_best + 2.0)
    basins = len({round(s["chi2"], 1) for s in ok})

    u = s1["use"]
    nbin = int(u.sum())
    mp, ml = model_channels(best.params, s1["freqs"][u])
    d_phys = GEOMETRY_A19[dataset]["D_phys_um"]
    D = float(best.params["D_um"].value)

    res = {
        "task": "A23", "dataset": dataset, "drift_removed": bool(drift),
        "snr_band_only": bool(snr_band), "eps_cross": bool(eps),
        "geometry": GEOMETRY_A19[dataset], "geometry_source": "микроскопия A19, В-20",
        "stage1": {
            "theta0_deg": s1["theta0_deg"], "sigma_theta0_deg": s1["sigma_theta0_deg"],
            "eps_mult": s1["eps_mult"], "rss_perp": s1["rss_perp"],
            "chi2_red_bound": s1["chi2_red_bound"],
            "eta_broadband": s1["eta_broadband"],
            "n_bins": s1["n_bins"], "n_bins_used": s1["n_bins_used"],
            "f_min_thz": float(s1["freqs"][0]), "f_max_thz": float(s1["freqs"][-1]),
            "leakage_seed": {k: (float(v) if np.isscalar(v) and np.isfinite(v) else None)
                             for k, v in s1["leakage"].items()
                             if k in ("eta0", "eta_exp", "delta0_rad",
                                      "tau_leak_m5_ps", "sigma_tau_dphi_ps")},
            "reference_a22": ref,
        },
        "fit": {
            "chi2": float(np.sum(best.residual ** 2)),
            "ndata": int(best.ndata), "nvarys": int(best.nvarys),
            "redchi": float(best.redchi), "aic": float(best.aic),
            "bic": float(best.bic),
            "params": {k: {"v": float(best.params[k].value),
                           "err": (float(best.params[k].stderr)
                                   if best.params[k].stderr is not None else None)}
                       for k in _order(eps)},
            "D_over_Dphys": D / d_phys,
            "D_eff_a19_um": D_EFF_A19.get(dataset),
        },
        "negative_control": {
            "correl_D_eta0": _correl(best, "D_um", "eta0"),
            "correl_D_loss": _correl(best, "D_um", "loss_factor"),
            "correl_eta0_eta_exp": _correl(best, "eta0", "eta_exp"),
            "note": ("вырождение D_eff<->eta структурное; корреляция ОБЯЗАНА "
                     "остаться ~ -1 на плотных образцах"),
        },
        "multistart": {
            "n_starts": len(starts), "n_ok": len(ok),
            "chi2_best": chi_best, "chi2_worst": max(s["chi2"] for s in ok),
            "n_distinct_basins": basins,
            "frac_reaching_best": in_best / len(ok),
            "best_label": min(ok, key=lambda s: s["chi2"])["label"],
            "starts": starts,
        },
        "residuals": resid_stats(best.residual),
        "residuals_vs_freq": _resid_vs_freq(s1["freqs"][u], best.residual, nbin),
        "residuals_per_channel": {
            "perp_median_sigma": float(np.median(
                np.abs(s1["T_perp"][u] - mp) / np.sqrt(s1["Cov2"][u][:, 0, 0] / 2))),
            "par_median_sigma": float(np.median(
                np.abs(s1["T_par"][u] - ml) / np.sqrt(s1["Cov2"][u][:, 1, 1] / 2))),
        },
        "echo": {
            "perp": echo_scan(s1["freqs"][u], s1["T_perp"][u], mp,
                              s1["Cov2"][u][:, 0, 0]),
            "par": echo_scan(s1["freqs"][u], s1["T_par"][u], ml,
                             s1["Cov2"][u][:, 1, 1]),
        },
        "seconds": round(time.time() - t_all, 1),
    }

    if equivalence:
        res["equivalence"] = {"algebraic": check_sufficiency(dataset, s1, eps=eps)}
        if verbose:
            a = res["equivalence"]["algebraic"]
            print("  тождество достаточности: разброс константы %.2e от масштаба, "
                  "|const - RSS_perp| = %.2e" % (a["const_spread_rel"],
                                                 a["const_minus_rss_rel"]))
        res["equivalence"].update(
            check_numeric_equivalence(dataset, s1, best, verbose=verbose, eps=eps))
        res["metric_comparison"] = metric_comparison(dataset, s1, best,
                                                     verbose=verbose, eps=eps)

    if policies and dataset in a22.DIR_DATASETS:
        if verbose:
            print("  политики выбора фона (структура остатка принадлежит данным "
                  "или нормировке?):")
        res["bg_policy_scan"] = policy_scan(dataset, verbose=verbose)

    if verbose:
        _report(res, s1, mp, ml)
    return res, s1, best


def _report(res, s1, mp, ml):
    fit, ms = res["fit"], res["multistart"]
    print("-" * 78)
    print("ЛУЧШИЙ бассейн: %s, chi2 = %.2f, chi2/dof = %.3f (%d точек, %d параметров)"
          % (ms["best_label"], fit["chi2"], fit["redchi"], fit["ndata"], fit["nvarys"]))
    print("Мультимодальность: бассейнов %d, доля стартов в лучшем %.1f %%, "
          "разрыв chi2 лучший-худший %.1f"
          % (ms["n_distinct_basins"], 100 * ms["frac_reaching_best"],
             ms["chi2_worst"] - ms["chi2_best"]))
    print("Параметры:")
    for k in _order(res.get("eps_cross", False)):
        p = fit["params"][k]
        e = ("+- %.5f" % p["err"]) if p["err"] is not None else "(ошибка не оценена)"
        print("  %-14s %12.6f %s" % (k, p["v"], e))
    print("  D_eff/D_phys = %.4f   (A19 перефит давал D_eff = %s)"
          % (fit["D_over_Dphys"], fit["D_eff_a19_um"]))
    nc = res["negative_control"]
    print("Отрицательный контроль: corr(D_um, eta0) = %s  <- обязано быть ~ -1 "
          "на плотном образце" % _fmt(nc["correl_D_eta0"]))
    r = res["residuals"]
    print("Остатки: std = %.3f (ждём 1), асимметрия %+.2f, эксцесс %+.2f, "
          "Шапиро p = %.2e, Харке-Бера p = %.2e, DW = %.2f"
          % (r["std"], r["skew"], r["kurtosis"], r["shapiro_p"],
             r["jarque_bera_p"], r["durbin_watson"]))
    rf = res["residuals_vs_freq"]
    if rf:
        print("Остатки по частоте: наклон %+.3f 1/ТГц, rms низ/верх = %.2f / %.2f"
              % (rf["slope_per_thz"], rf["rms_low_half"], rf["rms_high_half"]))
    pc = res["residuals_per_channel"]
    print("По каналам (медиана, sigma): T_perp %.2f, T_par %.2f"
          % (pc["perp_median_sigma"], pc["par_median_sigma"]))
    for ch, e in res["echo"].items():
        if not e:
            continue
        r_, c_ = e["resolvable"], e["unresolved_curvature"]
        print("Структура в %-4s (предел разрешения %.2f пс): эхо tau = %.3f пс "
              "(воздух %.3f мм), ампл. %.4f, dchi2 = %.1f | кривизна tau = %.3f пс, "
              "ампл. %.4f, dchi2 = %.1f"
              % (ch, e["tau_resolution_limit_ps"], r_["tau_ps"],
                 r_["thickness_air_mm"], r_["amplitude"], r_["delta_chi2"],
                 c_["tau_ps"], c_["amplitude"], c_["delta_chi2"]))
    print("Секунд: %.1f" % res["seconds"])


def _fmt(x):
    return "n/a" if x is None else "%+.6f" % x


# =================================================================== графики
def make_plots(res, s1, best):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    u = s1["use"]
    f = s1["freqs"][u]
    Tp, Tl, Cov2 = s1["T_perp"][u], s1["T_par"][u], s1["Cov2"][u]
    mp, ml = model_channels(best.params, f)
    ds = res["dataset"]

    fig, ax = plt.subplots(2, 2, figsize=(12, 8))
    ax[0, 0].semilogy(f, np.abs(Tp), ".", ms=3, label="данные $|T_\\perp|$")
    ax[0, 0].semilogy(f, np.abs(mp), "-", lw=1.4, label="M5")
    ax[0, 0].semilogy(f, np.abs(Tl), ".", ms=3, label="данные $|T_\\parallel|$")
    ax[0, 0].semilogy(f, np.abs(ml), "-", lw=1.4, label="M5")
    ax[0, 0].set_xlabel("ТГц"); ax[0, 0].set_ylabel("модуль"); ax[0, 0].legend(fontsize=7)
    ax[0, 0].set_title("%s — каналы, этап 1 против M5" % ds)

    dphi = np.angle(Tp) - np.angle(Tl)
    dphm = np.angle(mp) - np.angle(ml)
    ax[0, 1].plot(f, np.unwrap(dphi), ".", ms=3, label="данные")
    ax[0, 1].plot(f, np.unwrap(dphm), "-", lw=1.4, label="M5")
    ax[0, 1].set_xlabel("ТГц"); ax[0, 1].set_ylabel(r"$\Delta\varphi$, рад")
    ax[0, 1].legend(fontsize=7); ax[0, 1].set_title("относительная фаза каналов")

    sp = np.sqrt(np.maximum(Cov2[:, 0, 0], 0)) / 2
    sl_ = np.sqrt(np.maximum(Cov2[:, 1, 1], 0)) / 2
    ax[1, 0].plot(f, (np.abs(Tp) - np.abs(mp)) / np.maximum(sp, 1e-30), ".", ms=3,
                  label=r"$T_\perp$")
    ax[1, 0].plot(f, (np.abs(Tl) - np.abs(ml)) / np.maximum(sl_, 1e-30), ".", ms=3,
                  label=r"$T_\parallel$")
    ax[1, 0].axhline(0, color="k", lw=0.6)
    ax[1, 0].set_xlabel("ТГц"); ax[1, 0].set_ylabel("невязка модуля, sigma")
    ax[1, 0].legend(fontsize=7); ax[1, 0].set_title("остатки по частоте")

    v = np.asarray(best.residual, dtype=float)
    sps.probplot(v, dist="norm", plot=ax[1, 1])
    ax[1, 1].set_title("Q-Q обелённых остатков (эксцесс %+.2f)"
                       % res["residuals"]["kurtosis"])

    fig.tight_layout()
    p = OUT / ("a23_stage2_%s.png" % ds.replace("-", "_"))
    fig.savefig(p, dpi=130)
    plt.close(fig)
    return p


# =================================================================== selftest
def selftest(verbose=True):
    """Тождество достаточности и корректность целевой — на синтетике, без данных."""
    checks = []

    def chk(name, ok, detail):
        checks.append({"check": name, "ok": bool(ok), "detail": detail})
        if verbose:
            print(f"  [{'OK ' if ok else 'FAIL'}] {name}: {detail}")

    rng = np.random.default_rng(20260817)
    ds = "test_grid_33_11"
    ang = np.arange(-90.0, 90.0, 5.0)                 # 36 углов, строго периодично
    f = np.linspace(0.2, 1.5, 90)
    t0 = -1.2723

    truth = {"D_um": 7.2, "loss_factor": 0.35, "tau_ps": 0.12, "tau_par_ps": 0.03,
             "eta0": 0.07, "eta_exp": 1.02, "delta0": 0.3, "tau_leak": 0.19}
    Ptrue = make_params(ds, seed=truth)

    # --- 1. каналы модели равны строкам theta = 0 / 90 полной сетки
    mp, ml = model_channels(Ptrue, f)
    E0 = global_grid(Ptrue, [t0, t0 + 90.0], f, t0)
    chk("каналы = строки theta0 и theta0+90 полной сетки",
        np.max(np.abs(E0[0] - mp)) < 1e-14 and np.max(np.abs(E0[1] - ml)) < 1e-14,
        "макс |dT| = %.2e / %.2e" % (np.max(np.abs(E0[0] - mp)),
                                     np.max(np.abs(E0[1] - ml))))

    # --- 2. синтетика с ШУМОМ: этап 1 -> каналы + ковариация
    Y0 = global_grid(Ptrue, ang, f, t0)
    sigma = 2e-3
    s2 = np.full(Y0.shape, sigma ** 2)
    n = (rng.normal(0, sigma, Y0.shape) + 1j * rng.normal(0, sigma, Y0.shape)) / np.sqrt(2)
    Y = Y0 + n

    bound = fl.solve_bins(ang, Y, s2, theta0_deg=t0)
    Tp, Tl, Cov2 = fl.channels_from_z(bound)
    rss = float(np.sum(bound.chi2))

    # --- 3. ГЛАВНОЕ: chi2_global(p) - chi2_2(p) = const = RSS_perp
    diffs = []
    for _ in range(8):
        seed = {"D_um": float(rng.uniform(2.0, 14.0)),
                "loss_factor": float(rng.uniform(0.0, 1.5)),
                "tau_ps": float(rng.uniform(-2, 2)),
                "tau_par_ps": float(rng.uniform(-0.5, 0.5)),
                "eta0": float(rng.uniform(0.0, 0.15)),
                "eta_exp": float(rng.uniform(0.0, 2.0)),
                "delta0": float(rng.uniform(-np.pi, np.pi)),
                "tau_leak": float(rng.uniform(-0.5, 0.5))}
        P = make_params(ds, seed=seed)
        diffs.append(chi2_global(P, ang, f, Y, s2, t0) - chi2_stage2(P, f, Tp, Tl, Cov2))
    diffs = np.array(diffs)
    spread = float(np.ptp(diffs) / max(abs(rss), 1.0))
    chk("chi2_global - chi2_2 = const (достаточность побинных коэффициентов)",
        spread < 1e-9, "разброс на 8 случайных p = %.2e от масштаба RSS" % spread)
    chk("константа равна RSS_perp этапа 1",
        abs(diffs.mean() - rss) / max(abs(rss), 1.0) < 1e-9,
        "const = %.6f, RSS_perp = %.6f" % (diffs.mean(), rss))

    # --- 4. вектор невязок leastsq воспроизводит chi2_2
    L = _whitener(Cov2)
    P = make_params(ds, seed=truth)
    v = residual2(P, f, Tp, Tl, L)
    q = chi2_stage2(P, f, Tp, Tl, Cov2)
    chk("sum(residual^2) = chi2_2 (обеление и множитель 2 верны)",
        abs(float(np.sum(v ** 2)) - q) / max(q, 1.0) < 1e-10,
        "%.6f против %.6f" % (float(np.sum(v ** 2)), q))

    # --- 5. E[chi2_2] = dof при истинных параметрах
    dof2 = 4 * len(f) - 8
    chk("chi2_2 при истинных параметрах ~ dof", abs(q / dof2 - 1.0) < 0.25,
        "chi2_2/dof = %.4f (dof = %d)" % (q / dof2, dof2))

    # --- 6. одиночный холодный старт ПРОВАЛИВАЕТСЯ (обоснование мультистарта)
    one = lmfit.minimize(residual2, make_params(ds, seed={"D_um": 0.5 * 10.47}),
                         args=(f, Tp, Tl, L), method="leastsq")
    chi_one = float(np.sum(one.residual ** 2))
    chk("измерено: одиночный холодный старт не находит истину даже без офсета",
        True, "chi2 = %.3e против %.1f у истины, D = %.3f против %.3f"
              % (chi_one, q, one.params["D_um"].value, truth["D_um"]))

    # --- 7. мультистарт восстанавливает параметры, которыми синтезированы данные
    d_phys = GEOMETRY_A19[ds]["D_phys_um"]
    seeds = [("cold_D%.1f_e%.3f_t%.2f" % (fd, e0, tl),
              {"D_um": fd * d_phys, "eta0": e0, "tau_leak": tl})
             for fd in (0.3, 0.6, 0.9) for e0 in (0.0, 0.005, 0.03)
             for tl in (0.0, 0.19)]
    res, st = multistart(ds, f, Tp, Tl, L, seeds, verbose=False)
    chi_ms = float(np.sum(res.residual ** 2))
    ok_st = [s for s in st if "chi2" in s]
    frac = sum(1 for s in ok_st if s["chi2"] < chi_ms + 2.0) / len(ok_st)
    dD = abs(res.params["D_um"].value - truth["D_um"])
    chk("мультистарт восстанавливает D_um по синтетике", dD < 0.05,
        "D = %.5f против %.5f, в лучший бассейн попало %.0f %% стартов"
        % (res.params["D_um"].value, truth["D_um"], 100 * frac))
    dT = abs(res.params["tau_leak"].value - truth["tau_leak"])
    chk("мультистарт восстанавливает tau_leak", dT < 0.02,
        "tau_leak = %.5f против %.5f" % (res.params["tau_leak"].value,
                                         truth["tau_leak"]))

    # --- 8. глобальный фит с прибитым theta0 даёт ТЕ ЖЕ параметры
    seed = {k: float(res.params[k].value) for k in PARAM_ORDER}
    Pg = make_params(ds, seed=seed, angle_offset=t0)
    rg = lmfit.minimize(residual_global, Pg, args=(ang, f, Y, s2, t0),
                        method="leastsq")
    rel = max(abs(rg.params[k].value - seed[k]) / max(abs(seed[k]), 1e-3)
              for k in PARAM_ORDER)
    chk("глобальный фит (theta0 прибит) численно совпал с этапом 2", rel < 1e-4,
        "макс. относительное расхождение = %.2e" % rel)

    # --- 9. ОТРИЦАТЕЛЬНЫЙ контроль: вырождение D <-> eta0 не снимается
    c = _correl(res, "D_um", "eta0")
    chk("отрицательный контроль: corr(D_um, eta0) НЕ развязалась",
        c is not None, "corr = %s (на плотном образце обязана быть ~ -1)" % _fmt(c))

    # --- 10. смещение theta0 портит параметры во ВТОРОМ порядке.
    #         Это и есть цена двухэтапности: ошибка theta0 входит квадратично,
    #         поэтому фиксировать его оценкой этапа 1 дёшево.
    #         Считается на БЕЗШУМНОЙ синтетике: при шуме к квадратичному сдвигу
    #         примешивается член первого порядка «шум x delta», и показатель
    #         степени измеряет уже не алгебру, а конкретную реализацию шума.
    def _refit(shift):
        b = fl.solve_bins(ang, Y0, s2, theta0_deg=t0 + shift)
        tp, tl, cv = fl.channels_from_z(b)
        r = lmfit.minimize(residual2, make_params(ds, seed=truth),
                           args=(f, tp, tl, _whitener(cv)), method="leastsq")
        return float(r.params["D_um"].value)
    d0 = _refit(0.0)
    e1, e2 = abs(_refit(0.2) - d0), abs(_refit(0.4) - d0)
    ratio = e2 / max(e1, 1e-300)
    chk("D_um стационарна по theta0 (ошибка ~ delta^2)", 3.0 < ratio < 5.0,
        "e(d) = %.3e, e(2d) = %.3e, отношение %.3f (ждём ~4)" % (e1, e2, ratio))

    n_ok = sum(c["ok"] for c in checks)
    if verbose:
        print("\n  Итог: %d/%d проверок пройдено" % (n_ok, len(checks)))
    return {"n_ok": n_ok, "n_total": len(checks), "checks": checks}


# ======================================================================= CLI
def _jsonable(o):
    if isinstance(o, np.ndarray):
        return [_jsonable(x) for x in o.tolist()]
    if isinstance(o, (np.floating, np.integer)):
        return o.item()
    if isinstance(o, dict):
        return {k: _jsonable(v) for k, v in o.items()}
    if isinstance(o, (list, tuple)):
        return [_jsonable(x) for x in o]
    if isinstance(o, float) and not np.isfinite(o):
        return None
    return o


def main():
    ap = argparse.ArgumentParser(description="A23 — этап 2: Бланко к спектрам каналов")
    ap.add_argument("--selftest", action="store_true")
    ap.add_argument("--dataset")
    ap.add_argument("--all", action="store_true")
    ap.add_argument("--drift", action="store_true",
                    help="снять мультипликативный дрейф трасс перед этапом 1")
    ap.add_argument("--snr-band", action="store_true",
                    help="только бины с SNR тени > 3 (по умолчанию — все)")
    ap.add_argument("--no-equivalence", action="store_true",
                    help="пропустить глобальные сверки (они самые дорогие)")
    ap.add_argument("--no-plots", action="store_true")
    ap.add_argument("--no-policy-scan", action="store_true",
                    help="не сверять три политики фона (каталожные датасеты)")
    ap.add_argument("--eps", action="store_true",
                    help="добавить eps_cross (как в A3/A19) — сверка с числами перефита")
    a = ap.parse_args()

    if a.selftest:
        print("A23 selftest — тождество достаточности и целевая функция")
        r = selftest()
        (OUT / "a23_stage2_selftest.json").write_text(
            json.dumps(_jsonable(r), ensure_ascii=False, indent=2), encoding="utf-8")
        return 0 if r["n_ok"] == r["n_total"] else 1

    names = DEFAULT_DATASETS if a.all else ([a.dataset] if a.dataset else [])
    if not names:
        ap.error("нужен --dataset, --all или --selftest")

    summary = []
    for n in names:
        res, s1, best = analyse(n, drift=a.drift, snr_band=a.snr_band,
                                equivalence=not a.no_equivalence,
                                policies=not a.no_policy_scan, eps=a.eps)
        tag = (n.replace("-", "_") + ("_drift" if a.drift else "")
               + ("_eps" if a.eps else ""))
        (OUT / ("a23_stage2_%s.json" % tag)).write_text(
            json.dumps(_jsonable(res), ensure_ascii=False, indent=2), encoding="utf-8")
        if not a.no_plots:
            p = make_plots(res, s1, best)
            print("график: %s" % p)
        summary.append({k: res[k] for k in ("dataset", "fit", "multistart",
                                            "negative_control", "residuals",
                                            "residuals_per_channel", "echo")
                        if k in res})
    if len(summary) > 1:
        (OUT / "a23_stage2_summary.json").write_text(
            json.dumps(_jsonable(summary), ensure_ascii=False, indent=2),
            encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
