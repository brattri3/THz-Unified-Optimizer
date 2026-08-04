"""A3 — является ли `eps_cross` КОНСТАНТОЙ УСТАНОВКИ или per-dataset поглотителем невязки?

ВОПРОС (переформулирован после A2). В опорной модели M5 к комплексной амплитуде добавляется
угол-независимый член `eps_cross = eps_re + i*eps_im` — феноменологическая «кросс-поляризация
детектора». A2 показал, что физический однопараметрический кандидат на его роль (рассогласование
оси второго WGP, `phi2`) НЕ выдерживает мультистарта: бассейн `eps_cross` глубже на 110.5 AIC.
Остался вопрос об ИНТЕРПРЕТАЦИИ уцелевшего члена:

  * `eps_cross` СТАБИЛЕН по датасетам  -> это реальная константа тракта (детектор/оптика), её
    следует ПРИКОЛОТИТЬ или разделить между датасетами. Это разорвало бы вырождение `D_eff <-> eta`
    (корреляция +-1.0), которое блокирует цель №5 проекта.
  * `eps_cross` ГУЛЯЕТ  -> это per-dataset поглотитель остаточной структуры (паттерн HN10,
    переподгонка), и его нельзя трактовать физически.

ДИЗАЙН. Две группы датасетов с РАЗНОЙ оптической схемой (метаданные владельца 2026-07-27,
`data_pool/two_wgp_attenuator/README.md`):
  * АТТЕНЮАТОР (два реальных WGP): `att-11-16-356`, `att-11-16-s1`, `att-11-16-s2`,
    `att-11-16-s3` — это ОДИН прибор ATT-11-16-CA85 в четырёх условиях съёмки
    (переименование 2026-08-04; прежние имена `356att`, `series1..3`);
  * ОДИНОЧНЫЙ WGP между почти идеальными плёнками TYDEX: `specac`, `test_grid_40_20`.
Разделение существенно: если `eps_cross` — свойство ДЕТЕКТОРА, он обязан быть одинаков в обеих
группах; если свойство конкретной сборки — группы разойдутся. Это делает тест диагностическим,
а не просто описательным.

МУЛЬТИСТАРТ ОБЯЗАТЕЛЕН, И ПО СЕТКЕ, А НЕ В ДВЕ ТОЧКИ. Урок A2 обычно цитируют как «тёплый старт
лучше холодного», но это не так: A2 показал лишь, что ОДИН старт любого рода ненадёжен. Здесь
это воспроизведено количественно на 356att: старт `cold_eps0` (ровно затравка A2) даёт
AIC = -22684.4 и |eps| = 0.0345 — БИТ-В-БИТ опубликованное число A2; тёплый старт из `M5_1wgp`
даёт -22665.6; а старт с ненулевой затравкой `eps` достаёт -22792.1, т.е. на 107.7 AIC ГЛУБЖЕ
обоих. Поэтому стартов много и они покрывают сетку:
  * холодные: `eps` = 0; затем |eps| в {0.001, 0.01, 0.034} x фаза в {0, 90, 180, 270} град;
  * холодные по D: D_seed в {0.3, 0.6, 0.9} * D_phys (вырождение D <-> eta разводит бассейны
    именно по D, поэтому одной затравки D мало);
  * тёплые: из оптимума `M5_1wgp`, `eps` = 0 и |eps| = 0.034 x фаза в {0, 90, 180, 270} град.
Берётся ГЛУБОЧАЙШИЙ бассейн по AIC; в артефакт пишутся ВСЕ старты, разрыв между лучшим и худшим
и то, какое семейство победило — сами по себе меры мультимодальности задачи.

ТЕСТ ОДНОРОДНОСТИ. `eps_cross` — комплексное число, поэтому проверяется гипотеза равенства ПАРЫ
(Re, Im), а не модуля: chi2 = sum_i (z_i - z_bar)^T C_i^{-1} (z_i - z_bar), z_bar — взвешенное
среднее по обратным ковариациям, dof = 2*(N-1). Тест по модулю |eps| был бы смещён (модуль
неотрицателен, его распределение при малом SNR не гауссово).

ДИАПАЗОН УГЛОВ. Базовый прогон — на продакшн-диапазоне `config.ANGLES_LIMIT_2D=(0,90)`, чтобы числа
были прямо сравнимы с A2 и Фазами 1-3. Дополнительно, для датасетов с более широким покрытием, —
прогон на ПОЛНОМ диапазоне: это заодно закрывает пункт, явно оставленный незакрытым в A6
(«full-range M5 refit to see whether D_eff and eta move»).

СВЕРКА. Рядом печатается модельно-НЕЗАВИСИМОЕ `eta` из гармонического разбора A0
(`malus_linearization.channel_params`) — оценка утечки, не зависящая ни от Бланко, ни от Друде.
Расхождение между ней и `eta0` из M5 само по себе информативно.

Гардрейлы: оригиналы `data_pool/` и `unified_optimizer/` — только чтение; ядро (`fit_lib`,
`model_core`, `model_m5`, `t6_hn*`) импортируется, не правится; артефакты — в
`research/results/two_wgp/`.

Запуск ИЗ КОРНЯ репозитория:
    .venv/Scripts/python.exe research/two_wgp/a3_eps_stability.py
    ... --datasets att-11-16-356 att-11-16-s1 att-11-16-s2 att-11-16-s3 specac test_grid_40_20
    ... --no-fullrange        # только продакшн-диапазон (быстрее)
"""
from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

import numpy as np
import lmfit

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(REPO / "research" / "experiments"))
sys.path.insert(0, str(REPO))

import exp_builder                                          # noqa: E402  зона A
import malus_linearization as ml                            # noqa: E402  зона A
import fit_lib                                              # noqa: E402  зона B, ЧТЕНИЕ
import model_core                                           # noqa: E402  зона B, ЧТЕНИЕ
import t6_hn6                                               # noqa: E402  зона B, ЧТЕНИЕ
from t6_hn8 import compute_grid_hn8                         # noqa: E402  зона B, ЧТЕНИЕ
from unified_optimizer import config                        # noqa: E402  продакшн, ЧТЕНИЕ
from unified_optimizer.data_manager import DataManager       # noqa: E402
from unified_optimizer.optimizer_2d import get_transmission_spectra  # noqa: E402

OUT = REPO / "research" / "results" / "two_wgp"
OUT.mkdir(parents=True, exist_ok=True)
W_AMP, W_PHASE = fit_lib.W_AMP, fit_lib.W_PHASE

# Группировка по оптической схеме (метаданные владельца 2026-07-27).
GROUP = {
    "att-11-16-356": "attenuator", "att-11-16-s1": "attenuator",
    "att-11-16-s2": "attenuator", "att-11-16-s3": "attenuator",
    "specac": "single_wgp", "test_grid_40_20": "single_wgp",
    # purewave добавлен 2026-08-04 (A7), когда владелец домерил датасет; в прогоне A3
    # он ещё числился BLOCKED (файлы-заглушки нулевого размера).
    "purewave": "single_wgp",
}
DEFAULT_DATASETS = ["att-11-16-356", "att-11-16-s1", "att-11-16-s2", "att-11-16-s3",
                    "specac", "test_grid_40_20"]
# Затравка D_eff: известные значения M3 там, где есть; иначе половина D_phys
# (эмпирика «D_eff/D_phys = 0.43..0.64» из Фаз 1-3 — стартовая точка, НЕ пин).
SEED_D = {"att-11-16-356": 4.683, "test_grid_40_20": 11.447, "specac": 6.240}
EPS_SEED_ABS = 0.034            # |eps_cross| у 356att (A2) — масштаб альтернативных бассейнов


# ----------------------------------------------------------------- данные
def noise_floor(dataset, angles_limit="config"):
    """Порог шума по частотам — копия логики `t6_hn6.noise_floor` с ПАРАМЕТРОМ углов.

    Оригинал жёстко применяет `config.ANGLES_LIMIT_2D` к набору углов, по которым усредняется
    фон, поэтому на полном диапазоне он вернул бы порог, посчитанный по другому подмножеству
    трасс, чем сам фит. Здесь набор углов задаётся явно — ровно тем же, что идёт в фит.
    При `angles_limit="config"` результат совпадает с оригиналом (проверяется `--verify`).
    """
    data = DataManager(config.DATA_DIR).get_data_for_dataset(dataset)
    a_sorted = exp_builder.select_angles(data, angles_limit)
    bg, fc = [], None
    for a in a_sorted:
        t_s, E_s, t_b, E_b = data[a]
        f, s_s, s_b, tr = get_transmission_spectra(t_s, E_s, t_b, E_b)
        bg.append(s_b)
        fc = f if fc is None else fc
    bg_avg = np.mean(bg, axis=0)
    tail = np.where(fc >= 3.0)[0]
    if len(tail) == 0:
        tail = np.where(fc >= 2.5)[0]
    nfa = np.mean(bg_avg[tail])
    tnoise = 1.0 / ((bg_avg / nfa) ** 2)
    fmax = getattr(config, "F_MAX_2D", config.F_MAX)
    idx = np.where((fc >= config.F_MIN) & (fc <= fmax))[0]
    return tnoise[idx]


def build_masked(dataset, angles_limit="config"):
    """Сетка эксперимента + водяная маска + шумовая маска динамического диапазона.

    Эквивалент `model_m5.build_experiment_masked`, но через приватный билдер зоны A, чтобы
    угловой диапазон был параметром. `mask_angles=(0,90)` фиксирует водяную маску ровно как
    в Фазах 1-3 (находка A6: маска считается по среднему ПО УГЛАМ, т.е. молча зависит от их
    набора; фиксация отделяет эффект углов от эффекта маски).
    """
    ang, freqs, exp, valid, meta = exp_builder.load(
        dataset, angles_limit=angles_limit, mask_angles=(0.0, 90.0))
    tnoise = noise_floor(dataset, angles_limit)
    dr = np.abs(exp) > 1.5 * np.sqrt(np.broadcast_to(tnoise[None, :], exp.shape))
    return ang, freqs, exp, (valid & dr), meta


# ----------------------------------------------------------------- модель
def make_params(dataset, *, eps, seed=None):
    """Параметры M5 (+ опционально eps_cross). `gamma` приколочен (HN11: вырожден)."""
    geo = fit_lib.GEOMETRY[dataset]
    p_um = geo["P_um"]
    d_seed = SEED_D.get(dataset, 0.5 * geo["D_phys_um"])
    P = lmfit.Parameters()
    P.add("P_um", value=p_um, vary=False)
    P.add("D_um", value=d_seed, min=0.1, max=p_um - 0.5)
    P.add("loss_factor", value=0.3, min=0.0, max=5.0)
    P.add("gamma", value=2.0, vary=False)
    P.add("angle_offset", value=0.0, min=-10, max=10)
    P.add("tau_ps", value=0.0, min=-10, max=10)
    P.add("tau_par_ps", value=0.0, min=-5, max=5)
    P.add("eta0", value=0.02, min=0.0, max=1.0)
    P.add("eta_exp", value=1.0, min=-2, max=6)
    P.add("delta0", value=0.0, min=-np.pi, max=np.pi)
    P.add("tau_leak", value=0.0, min=-5, max=5)
    P.add("eps_re", value=0.0, min=-0.2, max=0.2, vary=eps)
    P.add("eps_im", value=0.0, min=-0.2, max=0.2, vary=eps)
    if seed:
        for k, v in seed.items():
            if k in P and v is not None:
                P[k].set(value=float(min(max(v, P[k].min), P[k].max)))
    return P


def model_grid(params, angles, freqs):
    """M5 (одиночный WGP, ядро `compute_grid_hn8`) + угол-независимый eps_cross."""
    E = compute_grid_hn8(angles, freqs,
                         params["P_um"].value * 1e-6, params["D_um"].value * 1e-6,
                         params["loss_factor"].value, params["angle_offset"].value,
                         params["tau_ps"].value, params["gamma"].value,
                         params["tau_par_ps"].value, params["eta0"].value,
                         params["eta_exp"].value, params["delta0"].value,
                         params["tau_leak"].value)
    if params["eps_re"].vary or params["eps_im"].vary:
        E = E + (params["eps_re"].value + 1j * params["eps_im"].value)
    return E


def residual(params, angles, freqs, exp, valid):
    if params["D_um"].value >= params["P_um"].value:
        return np.ones(int(np.sum(valid)) * 2) * 1e6
    theo = model_grid(params, angles, freqs)
    em, tm = exp[valid], theo[valid]
    amp = np.abs(em) - np.abs(tm)
    ph = np.angle(em) - np.angle(tm)
    ph = np.arctan2(np.sin(ph), np.cos(ph))
    return np.concatenate([amp * W_AMP, ph * W_PHASE])


def _metrics(res, dataset, angles, freqs, exp, valid, label, t0):
    theo = model_grid(res.params, angles, freqs)
    ampdb = np.where(valid, 20 * np.log10(np.maximum(np.abs(exp), 1e-12))
                     - 20 * np.log10(np.maximum(np.abs(theo), 1e-12)), np.nan)
    hi = angles >= np.percentile(angles, 75)          # «глубокая тень»
    amp = (np.abs(exp) - np.abs(theo))[valid]
    x = amp[np.isfinite(amp)]
    dw = float(np.sum(np.diff(x) ** 2) / np.sum(x ** 2)) if np.sum(x ** 2) > 0 else None
    D = float(res.params["D_um"].value)
    d_phys = fit_lib.GEOMETRY[dataset]["D_phys_um"]
    return {
        "label": label, "dataset": dataset,
        "nvarys": int(res.nvarys), "npts": int(np.sum(valid)),
        "aic": float(res.aic), "bic": float(res.bic), "redchi": float(res.redchi),
        "rmse_amp_db": float(np.sqrt(np.nanmean(ampdb[np.isfinite(ampdb)] ** 2))),
        "deepshadow_bias_db": float(np.nanmean(ampdb[hi, :])),
        "deepshadow_rmse_db": float(np.sqrt(np.nanmean(ampdb[hi, :] ** 2))),
        "durbin_watson": dw,
        "D_um": D, "D_phys_um": float(d_phys), "D_over_Dphys": D / float(d_phys),
        "seconds": round(time.time() - t0, 1),
        "params": {k: {"v": float(v.value),
                       "err": (float(v.stderr) if v.stderr is not None else None),
                       "vary": bool(v.vary)} for k, v in res.params.items()},
        "correl": {
            "D~eta0": _correl(res, "D_um", "eta0"),
            "eps_re~D": _correl(res, "eps_re", "D_um"),
            "eps_re~eta0": _correl(res, "eps_re", "eta0"),
            "eps_re~eps_im": _correl(res, "eps_re", "eps_im"),
        },
    }


def _correl(res, a, b):
    p = res.params.get(a)
    if p is None or not p.correl or b not in p.correl:
        return None
    return float(p.correl[b])


def run_fit(dataset, angles, freqs, exp, valid, *, eps, seed=None, label=""):
    t0 = time.time()
    P = make_params(dataset, eps=eps, seed=seed)
    res = lmfit.Minimizer(residual, P,
                          fcn_args=(angles, freqs, exp, valid)).minimize(method="leastsq")
    out = _metrics(res, dataset, angles, freqs, exp, valid, label, t0)
    out["_lmfit"] = res
    return out


def eps_with_error(row):
    """|eps|, arg(eps) и их sigma по дельта-методу из ковариации (Re, Im).

    Ковариация берётся из lmfit (масштабированная на redchi, как принято в пакете); если
    оценить её не удалось (вырождение), возвращаются None — такой датасет исключается из
    теста однородности, а не молча вносит нулевой вес.
    """
    pr, pi = row["params"]["eps_re"], row["params"]["eps_im"]
    re, im = pr["v"], pi["v"]
    mod = float(np.hypot(re, im))
    arg = float(np.degrees(np.arctan2(im, re)))
    cov = row.get("cov_eps")
    if cov is None or mod <= 0:
        return {"abs": mod, "abs_err": None, "arg_deg": arg, "arg_err_deg": None,
                "re": re, "im": im, "cov": None}
    C = np.asarray(cov, float)
    g_mod = np.array([re, im]) / mod
    g_arg = np.array([-im, re]) / mod ** 2
    var_mod = float(g_mod @ C @ g_mod)
    var_arg = float(g_arg @ C @ g_arg)
    return {"abs": mod, "abs_err": float(np.sqrt(max(var_mod, 0.0))),
            "arg_deg": arg, "arg_err_deg": float(np.degrees(np.sqrt(max(var_arg, 0.0)))),
            "re": re, "im": im, "cov": C.tolist()}


def homogeneity_test(entries):
    """Гипотеза «eps_cross одинаков во всех датасетах» для комплексной пары (Re, Im).

    chi2 = sum_i (z_i - z_bar)^T C_i^{-1} (z_i - z_bar), z_bar — взвешенное обратными
    ковариациями среднее, dof = 2*(N-1). Возвращает также «разброс/ошибка» — во сколько раз
    наблюдаемый разброс превышает заявленную точность.
    """
    from scipy import stats
    use = [e for e in entries if e["cov"] is not None]
    if len(use) < 2:
        return {"n": len(use), "verdict": "недостаточно датасетов с оценённой ковариацией"}
    Ws, Wz = np.zeros((2, 2)), np.zeros(2)
    for e in use:
        Ci = np.linalg.inv(np.asarray(e["cov"], float))
        z = np.array([e["re"], e["im"]])
        Ws += Ci
        Wz += Ci @ z
    zbar = np.linalg.solve(Ws, Wz)
    chi2 = 0.0
    per = []
    for e in use:
        Ci = np.linalg.inv(np.asarray(e["cov"], float))
        d = np.array([e["re"], e["im"]]) - zbar
        c = float(d @ Ci @ d)
        chi2 += c
        per.append({"dataset": e["dataset"], "chi2_contrib": c,
                    "mahalanobis_sigma": float(np.sqrt(c))})
    dof = 2 * (len(use) - 1)
    p = float(stats.chi2.sf(chi2, dof))
    mods = np.array([np.hypot(e["re"], e["im"]) for e in use])
    errs = np.array([e["abs_err"] if e["abs_err"] else np.nan for e in use])
    return {
        "n": len(use), "datasets": [e["dataset"] for e in use],
        "eps_mean_re": float(zbar[0]), "eps_mean_im": float(zbar[1]),
        "eps_mean_abs": float(np.hypot(*zbar)),
        "chi2": float(chi2), "dof": int(dof), "chi2_red": float(chi2 / dof), "p_value": p,
        "abs_spread_over_median_err": float(np.nanstd(mods, ddof=1) / np.nanmedian(errs))
        if np.isfinite(np.nanmedian(errs)) and np.nanmedian(errs) > 0 else None,
        "abs_min": float(mods.min()), "abs_max": float(mods.max()),
        "per_dataset": per,
        "verdict": ("СТАБИЛЕН (совместим с одной константой)" if p > 0.05 else
                    "НЕ стабилен (гипотеза общей константы отвергнута)"),
    }


# ------------------------------------------------------------ сетка стартов
def start_grid(dataset, warm):
    """Список (имя, seed) для мультистарта. Имя начинается с `cold`/`warm` — по нему
    определяется семейство. `seed=None` = штатная затравка `make_params` (затравка A2)."""
    d_phys = fit_lib.GEOMETRY[dataset]["D_phys_um"]
    d_max = fit_lib.GEOMETRY[dataset]["P_um"] - 0.5
    starts = [("cold_eps0", None)]
    for mod in (0.001, 0.01, EPS_SEED_ABS):
        for ph in (0, 90, 180, 270):
            starts.append((f"cold_eps{mod:g}_ph{ph}",
                           {"eps_re": mod * np.cos(np.deg2rad(ph)),
                            "eps_im": mod * np.sin(np.deg2rad(ph))}))
    for frac in (0.3, 0.6, 0.9):
        starts.append((f"cold_D{frac:g}Dphys",
                       {"D_um": min(frac * d_phys, d_max),
                        "eps_re": 0.001, "eps_im": 0.0}))
    starts.append(("warm_eps0", dict(warm, eps_re=0.0, eps_im=0.0)))
    for ph in (0, 90, 180, 270):
        starts.append((f"warm_eps_ph{ph}",
                       dict(warm,
                            eps_re=EPS_SEED_ABS * np.cos(np.deg2rad(ph)),
                            eps_im=EPS_SEED_ABS * np.sin(np.deg2rad(ph)))))
    return starts


# ------------------------------------------------------------- один датасет
def analyse(dataset, angles_limit="config", verbose=True):
    ang, freqs, exp, valid, meta = build_masked(dataset, angles_limit)
    if verbose:
        print(f"\n=== {dataset} ({GROUP.get(dataset, '?')}), углы "
              f"{ang.min():.0f}..{ang.max():.0f} (n={len(ang)}), частот {len(freqs)}, "
              f"валидных точек {int(np.sum(valid))} ===")

    # 1. подмодель без eps — база AIC и источник тёплого старта
    base = run_fit(dataset, ang, freqs, exp, valid, eps=False, label="M5_1wgp")
    warm = {k: v["v"] for k, v in base["params"].items()}

    # 2. мультистарт M5_eps по СЕТКЕ, оба семейства (см. докстринг: тёплый ≠ лучший).
    #    Холодный = штатная затравка make_params, как в A2 (`cold_eps0` её воспроизводит).
    starts = start_grid(dataset, warm)
    tried = []
    for name, seed in starts:
        r = run_fit(dataset, ang, freqs, exp, valid, eps=True, seed=seed,
                    label=f"M5_eps[{name}]")
        tried.append(r)

    # различные бассейны: старты, сошедшиеся в одну точку, склеиваются (AIC с точностью 0.5)
    basins = []
    for r in sorted(tried, key=lambda x: x["aic"]):
        for b in basins:
            if abs(b["aic"] - r["aic"]) < 0.5:
                b["starts"].append(r["label"].split("[")[1].rstrip("]"))
                break
        else:
            basins.append({"aic": r["aic"], "D_um": r["D_um"],
                           "eta0": r["params"]["eta0"]["v"],
                           "eps_abs": float(np.hypot(r["params"]["eps_re"]["v"],
                                                     r["params"]["eps_im"]["v"])),
                           "starts": [r["label"].split("[")[1].rstrip("]")]})
    # Разброс параметров среди СТАТИСТИЧЕСКИ НЕРАЗЛИЧИМЫХ бассейнов. Прямая мера
    # идентифицируемости: если внутри dAIC<=2 (общепринятый порог «модели неразличимы»)
    # |eps| меняется в разы, то данные его не фиксируют — и «значение eps_cross у датасета»
    # не определено ДО всякого сравнения датасетов между собой.
    best_aic = min(b["aic"] for b in basins)
    span = {}
    for thr in (2.0, 10.0):
        near = [b for b in basins if b["aic"] - best_aic <= thr]
        vals = {k: [b[k] for b in near] for k in ("eps_abs", "D_um", "eta0")}
        span[str(int(thr))] = {
            "n_basins": len(near),
            "eps_abs_min": min(vals["eps_abs"]), "eps_abs_max": max(vals["eps_abs"]),
            "eps_abs_ratio": (max(vals["eps_abs"]) / min(vals["eps_abs"])
                              if min(vals["eps_abs"]) > 0 else None),
            "D_um_min": min(vals["D_um"]), "D_um_max": max(vals["D_um"]),
            "eta0_min": min(vals["eta0"]), "eta0_max": max(vals["eta0"]),
        }

    if verbose:
        print(f"    найдено различных бассейнов: {len(basins)} из {len(tried)} стартов")
        for b in basins:
            print(f"      AIC = {b['aic']:>11.1f}  D = {b['D_um']:>7.3f}  "
                  f"eta0 = {b['eta0']:.5f}  |eps| = {b['eps_abs']:.5f}  "
                  f"<- {len(b['starts'])} стартов ({', '.join(b['starts'][:3])}"
                  f"{', …' if len(b['starts']) > 3 else ''})")
    best = min(tried, key=lambda r: r["aic"])
    basin_gap = float(max(t["aic"] for t in tried) - best["aic"])
    # какое семейство стартов достало более глубокий бассейн — методологически важно
    fam = {}
    for f in ("cold", "warm"):
        sub = [t["aic"] for t in tried if t["label"].startswith(f"M5_eps[{f}")]
        fam[f] = min(sub) if sub else None
    win = ("cold" if fam["warm"] is None or (fam["cold"] is not None
                                             and fam["cold"] <= fam["warm"]) else "warm")
    start_family = {"best_cold_aic": fam["cold"], "best_warm_aic": fam["warm"],
                    "winner": win,
                    "cold_minus_warm_aic": (None if None in fam.values()
                                            else float(fam["cold"] - fam["warm"]))}

    # ковариация (Re, Im) из лучшего бассейна — нужна тесту однородности
    res = best.pop("_lmfit")
    cov_eps = None
    if res.covar is not None:
        names = [n for n in res.params if res.params[n].vary]
        if "eps_re" in names and "eps_im" in names:
            i, j = names.index("eps_re"), names.index("eps_im")
            cov_eps = [[float(res.covar[i, i]), float(res.covar[i, j])],
                       [float(res.covar[j, i]), float(res.covar[j, j])]]
    best["cov_eps"] = cov_eps
    for r in tried:
        r.pop("_lmfit", None)
    base.pop("_lmfit", None)

    best["dAIC_vs_M5_1wgp"] = best["aic"] - base["aic"]
    e = eps_with_error(best)
    e["dataset"] = dataset

    # 3. модельно-независимая сверка: eta из гармоник A0 (не зависит от Бланко/Друде)
    try:
        harm = ml.analyse_dataset(dataset, angles_limit=angles_limit, verbose=False)
        eta_free = harm["parseval_order4"]["eta_amplitude"]
        eta_free_med = harm["per_freq_summary"]["eta_median"]
        theta0_free = harm["parseval_order4"]["theta0_deg_from_h2"]
        theta0_free_err = harm["parseval_order4"]["theta0_err_deg"]
    except Exception as exc:                                   # noqa: BLE001
        eta_free = eta_free_med = theta0_free = theta0_free_err = None
        if verbose:
            print(f"    [гармоники недоступны] {type(exc).__name__}: {exc}")

    out = {
        "dataset": dataset, "group": GROUP.get(dataset), "angles_limit": meta["angles_limit"],
        "n_angles": meta["n_angles"], "n_freqs": meta["n_freqs"],
        "n_valid_points": int(np.sum(valid)),
        "base": base, "best": best, "starts": tried, "basins": basins,
        "n_basins": len(basins), "eps_span_within_dAIC": span,
        "basin_gap_aic": basin_gap, "start_family": start_family, "eps": e,
        "model_free": {"eta_parseval": eta_free, "eta_per_freq_median": eta_free_med,
                       "theta0_deg": theta0_free, "theta0_err_deg": theta0_free_err},
    }
    if verbose:
        b = best
        print(f"  ЛУЧШИЙ бассейн: {b['label']}, разрыв между стартами {basin_gap:.1f} AIC; "
              f"победило семейство '{start_family['winner']}'"
              + (f" (холодный минус тёплый = {start_family['cold_minus_warm_aic']:+.1f} AIC)"
                 if start_family["cold_minus_warm_aic"] is not None else ""))
        print(f"  eps_cross = {e['re']:+.5f}{e['im']:+.5f}i  |eps| = {e['abs']:.5f}"
              + (f" +- {e['abs_err']:.5f}" if e['abs_err'] else " (ошибка не оценена)")
              + f"  arg = {e['arg_deg']:+.1f}"
              + (f" +- {e['arg_err_deg']:.1f} град" if e['arg_err_deg'] else " град"))
        print(f"  dAIC(M5_eps - M5_1wgp) = {b['dAIC_vs_M5_1wgp']:+.1f}   "
              f"D_eff = {b['D_um']:.3f} мкм (D/D_phys = {b['D_over_Dphys']:.3f})   "
              f"eta0 = {b['params']['eta0']['v']:.5f}")
        print(f"  tau_leak = {b['params']['tau_leak']['v']:+.4f} пс   "
              f"angle_offset = {b['params']['angle_offset']['v']:+.3f} град   "
              f"DW = {b['durbin_watson']:.3f}   DS_rmse = {b['deepshadow_rmse_db']:.3f} дБ")
        print(f"  корреляции: D~eta0 = {_fmt(b['correl']['D~eta0'])}, "
              f"eps_re~D = {_fmt(b['correl']['eps_re~D'])}, "
              f"eps_re~eta0 = {_fmt(b['correl']['eps_re~eta0'])}")
        for thr in ("2", "10"):
            s = span[thr]
            print(f"  внутри dAIC<={thr:>2}: {s['n_basins']} бассейнов, "
                  f"|eps| {s['eps_abs_min']:.5f}..{s['eps_abs_max']:.5f}"
                  + (f" (в {s['eps_abs_ratio']:.1f} раза)" if s['eps_abs_ratio'] else "")
                  + f", D {s['D_um_min']:.3f}..{s['D_um_max']:.3f}, "
                  f"eta0 {s['eta0_min']:.5f}..{s['eta0_max']:.5f}")
        if eta_free is not None:
            print(f"  модельно-независимо (гармоники A0): eta = {eta_free:.5f} "
                  f"(медиана по частотам {eta_free_med:.5f}), "
                  f"theta0 = {theta0_free:+.3f} +- {theta0_free_err:.3f} град")
    return out


def _fmt(x, fmt="{:+.3f}"):
    return fmt.format(x) if x is not None else "n/a"


# ---------------------------------------------------------------- проверки
def verify(verbose=True):
    """Сверка приватных копий с оригиналами (гарантия «не форкнул логику молча»)."""
    ok = True
    dm = DataManager(config.DATA_DIR)
    for ds in ("att-11-16-356", "att-11-16-s1", "specac"):
        if not dm.get_data_for_dataset(ds):
            continue
        mine = noise_floor(ds, "config")
        ref = t6_hn6.noise_floor(ds)
        same = np.array_equal(mine, ref)
        ok &= same
        if verbose:
            print(f"  noise_floor({ds:<16}) vs t6_hn6: "
                  f"{'бит-в-бит' if same else 'РАСХОЖДЕНИЕ'}")
    import model_m5 as m5
    for ds in ("att-11-16-356", "att-11-16-s1"):
        if not dm.get_data_for_dataset(ds):
            continue
        a1, f1, e1, v1, _ = build_masked(ds, "config")
        a2, f2, e2, v2 = m5.build_experiment_masked(ds)
        same = all(np.array_equal(x, y) for x, y in
                   ((a1, a2), (f1, f2), (e1, e2), (v1, v2)))
        ok &= same
        if verbose:
            print(f"  build_masked({ds:<16}) vs model_m5: "
                  f"{'бит-в-бит' if same else 'РАСХОЖДЕНИЕ'}")
    return ok


# -------------------------------------------------------------------- main
def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0],
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--datasets", nargs="*", default=DEFAULT_DATASETS)
    ap.add_argument("--no-fullrange", action="store_true",
                    help="только продакшн-диапазон ANGLES_LIMIT_2D=(0,90)")
    ap.add_argument("--verify", action="store_true",
                    help="только сверка приватных копий с оригиналами")
    ap.add_argument("--out", default=str(OUT / "a3_eps_stability.json"))
    args = ap.parse_args()

    if args.verify:
        print("Сверка приватных копий с оригиналами (допуск 0):")
        good = verify()
        print("ИТОГ:", "совпадение бит-в-бит" if good else "РАСХОЖДЕНИЕ!")
        sys.exit(0 if good else 1)

    print("Сверка приватных копий с оригиналами (допуск 0):")
    if not verify():
        print("ОСТАНОВ: приватные копии разошлись с оригиналами.")
        sys.exit(1)

    dm = DataManager(config.DATA_DIR)
    payload = {"task": "A3_refit_series", "session": "A", "nick": "A3",
               "question": "eps_cross — константа тракта или per-dataset поглотитель невязки?",
               "runs": {}, "skipped": {}}

    for mode, limit in (("config_0_90", "config"),
                        ("full_range", None)):
        if mode == "full_range" and args.no_fullrange:
            continue
        print(f"\n{'=' * 78}\nДИАПАЗОН УГЛОВ: {mode}\n{'=' * 78}")
        rows, entries = [], []
        for ds in args.datasets:
            data = dm.get_data_for_dataset(ds)
            if not data:
                payload["skipped"][ds] = "датасет не найден или пуст (нет непустых трасс)"
                print(f"\n[пропуск] {ds}: нет данных")
                continue
            if mode == "full_range":
                n_full = len(exp_builder.select_angles(data, None))
                n_cfg = len(exp_builder.select_angles(data, "config"))
                if n_full <= n_cfg:
                    print(f"\n[пропуск] {ds}: полный диапазон не даёт новых углов "
                          f"({n_full} = {n_cfg})")
                    continue
            try:
                r = analyse(ds, angles_limit=limit)
            except Exception as exc:                            # noqa: BLE001
                payload["skipped"][f"{ds}@{mode}"] = f"{type(exc).__name__}: {exc}"
                print(f"\n[пропуск] {ds}: {type(exc).__name__}: {exc}")
                continue
            rows.append(r)
            e = dict(r["eps"]); e["dataset"] = ds
            entries.append(e)

        hom_all = homogeneity_test(entries)
        hom_grp = {}
        for g in ("attenuator", "single_wgp"):
            sub = [e for e in entries if GROUP.get(e["dataset"]) == g]
            hom_grp[g] = homogeneity_test(sub)
        payload["runs"][mode] = {"datasets": rows,
                                 "homogeneity_all": hom_all,
                                 "homogeneity_by_group": hom_grp}
        _report(mode, rows, hom_all, hom_grp)

    # сравнение диапазонов: сдвинулись ли D_eff / eta / eps при снятии обрезки углов
    if "full_range" in payload["runs"]:
        cmp = {}
        by_cfg = {r["dataset"]: r for r in payload["runs"]["config_0_90"]["datasets"]}
        for r in payload["runs"]["full_range"]["datasets"]:
            ds = r["dataset"]
            if ds not in by_cfg:
                continue
            a, b = by_cfg[ds]["best"], r["best"]
            cmp[ds] = {
                "D_um": [a["D_um"], b["D_um"]],
                "D_over_Dphys": [a["D_over_Dphys"], b["D_over_Dphys"]],
                "eta0": [a["params"]["eta0"]["v"], b["params"]["eta0"]["v"]],
                "eps_abs": [by_cfg[ds]["eps"]["abs"], r["eps"]["abs"]],
                "eps_arg_deg": [by_cfg[ds]["eps"]["arg_deg"], r["eps"]["arg_deg"]],
                "angle_offset": [a["params"]["angle_offset"]["v"],
                                 b["params"]["angle_offset"]["v"]],
                "tau_leak": [a["params"]["tau_leak"]["v"], b["params"]["tau_leak"]["v"]],
                "n_angles": [by_cfg[ds]["n_angles"], r["n_angles"]],
            }
        payload["range_comparison"] = cmp
        print(f"\n{'=' * 78}\nЭФФЕКТ СНЯТИЯ ОБРЕЗКИ УГЛОВ (закрывает открытый пункт A6)\n{'=' * 78}")
        print(f"{'датасет':<18}{'углов':>9}{'D_eff':>17}{'D/Dphys':>15}"
              f"{'eta0':>17}{'|eps|':>17}")
        for ds, c in cmp.items():
            print(f"{ds:<18}{c['n_angles'][0]:>4}->{c['n_angles'][1]:<4}"
                  f"{c['D_um'][0]:>8.3f}->{c['D_um'][1]:<8.3f}"
                  f"{c['D_over_Dphys'][0]:>7.3f}->{c['D_over_Dphys'][1]:<7.3f}"
                  f"{c['eta0'][0]:>8.5f}->{c['eta0'][1]:<8.5f}"
                  f"{c['eps_abs'][0]:>8.5f}->{c['eps_abs'][1]:<8.5f}")

    Path(args.out).write_text(json.dumps(payload, ensure_ascii=False, indent=2),
                              encoding="utf-8")
    print(f"\nартефакт -> {args.out}")
    if payload["skipped"]:
        print("пропущено:", json.dumps(payload["skipped"], ensure_ascii=False))


def _report(mode, rows, hom_all, hom_grp):
    print(f"\n--- СВОДКА ({mode}) ---")
    print(f"{'датасет':<18}{'группа':<12}{'|eps|':>9}{'+-':>9}{'arg,град':>10}"
          f"{'dAIC':>10}{'D_eff':>9}{'D/Dphys':>9}{'eta0':>9}{'eta(A0)':>9}")
    for r in rows:
        e, b = r["eps"], r["best"]
        mf = r["model_free"]["eta_parseval"]
        print(f"{r['dataset']:<18}{r['group'] or '?':<12}{e['abs']:>9.5f}"
              f"{(e['abs_err'] if e['abs_err'] else float('nan')):>9.5f}"
              f"{e['arg_deg']:>+10.1f}{b['dAIC_vs_M5_1wgp']:>+10.1f}"
              f"{b['D_um']:>9.3f}{b['D_over_Dphys']:>9.3f}"
              f"{b['params']['eta0']['v']:>9.5f}"
              f"{(mf if mf is not None else float('nan')):>9.5f}")
    for name, h in (("ВСЕ датасеты", hom_all),
                    ("группа attenuator", hom_grp.get("attenuator", {})),
                    ("группа single_wgp", hom_grp.get("single_wgp", {}))):
        if "chi2" not in h:
            print(f"  однородность eps_cross [{name}]: {h.get('verdict', 'н/д')}")
            continue
        print(f"  однородность eps_cross [{name}]: chi2 = {h['chi2']:.1f} / {h['dof']} dof "
              f"(chi2_red = {h['chi2_red']:.1f}), p = {h['p_value']:.2e}  -> {h['verdict']}")
        print(f"      |eps| в диапазоне {h['abs_min']:.5f}..{h['abs_max']:.5f}, "
              f"взвеш. среднее |eps_bar| = {h['eps_mean_abs']:.5f}; "
              f"разброс/типичная ошибка = {_fmt(h['abs_spread_over_median_err'], '{:.1f}')}")


if __name__ == "__main__":
    main()
