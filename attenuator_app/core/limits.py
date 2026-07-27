"""Границы применимости (F1..F8) и консервативная оценка неопределённости.

Политика неопределённости прототипа — СОЗНАТЕЛЬНОЕ ЗАНИЖЕНИЕ ТОЧНОСТИ ради
надёжного покрытия 95 %:
  * угол: сигма отсчёта = деление/2 (вместо деление/sqrt(12)), плюс столько же
    на выставление нуля;
  * параметры модели: сигма из разброса между независимыми сеансами калибровки,
    а не из ковариации одного фита (для D_eff это в ~25 раз больше);
  * интервал по углу считается ПОДСТАНОВКОЙ границ в модель, а не линейным
    распространением ошибки: вблизи скрещённого положения A(θ) резко нелинейна,
    и линеаризация занижает интервал.
"""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from .forward import Setup, attenuation_db, attenuation_integral_db, output_polarization
from .passport import Passport

K95 = 1.959963985


@dataclass
class Check:
    code: str
    ok: bool
    message: str
    remedy: str = ""

    def __str__(self):
        mark = "OK " if self.ok else "!! "
        s = f"{mark}{self.code}: {self.message}"
        return s + (f"\n      -> {self.remedy}" if (self.remedy and not self.ok) else "")


@dataclass
class Uncertainty:
    sigma_db: float = 0.0
    lo_db: float = 0.0
    hi_db: float = 0.0
    parts: dict = field(default_factory=dict)


# ---------------------------------------------------------------------------
def _eval(theta, passport, setup, freqs, weight, ref_theta, integral, **kw):
    """Скалярное затухание в точке — ТА ЖЕ величина, что показывает api.attenuation.

    Спектральный режим усредняет по полосе (а не берёт первую частоту), иначе
    интервал неопределённости считался бы для одного числа, а показывался
    рядом с другим — и мог не содержать его вовсе.
    """
    if integral:
        return float(attenuation_integral_db(theta, freqs, weight, passport, setup,
                                             ref_theta1_deg=ref_theta, **kw)[0])
    db, _ = attenuation_db(theta, freqs, passport, setup, ref_theta1_deg=ref_theta, **kw)
    return float(np.mean(db))


def uncertainty(theta_deg, passport: Passport, setup: Setup, freqs, *,
                weight=None, ref_theta_deg=None, integral=False) -> Uncertainty:
    """Консервативный 95 % интервал затухания в точке theta_deg."""
    base = _eval(theta_deg, passport, setup, freqs, weight, ref_theta_deg, integral)
    parts, lo, hi = {}, base, base

    # (1) угловая шкала: подстановка границ, интервал асимметричен
    s_th = passport.scale.sigma_total_deg()
    if s_th > 0:
        vals = [_eval(theta_deg + d * K95 * s_th, passport, setup, freqs,
                      weight, ref_theta_deg, integral) for d in (-1.0, +1.0)]
        parts["angle_scale"] = (min(vals), max(vals))
        lo, hi = min(lo, *vals), max(hi, *vals)

    # (2) параметры модели: огибающая по D_eff, loss_factor, gamma
    for name, p, key, rng in (
        ("D_eff", passport.D_eff_um, "d_eff", None),
        ("loss", passport.loss_factor, "loss", None),
        ("gamma", passport.gamma, "gamma", passport.gamma_range),
    ):
        if rng is not None and rng[0] != rng[1]:
            cand = [rng[0], rng[1]]
        elif p.sigma > 0:
            cand = [p.value - K95 * p.sigma, p.value + K95 * p.sigma]
        else:
            continue
        if key == "d_eff":
            cand = [max(c, 0.05) for c in cand]
        if key == "loss":
            cand = [max(c, 0.0) for c in cand]
        vals = [_eval(theta_deg, passport, setup, freqs, weight, ref_theta_deg,
                      integral, **{key: c}) for c in cand]
        parts[name] = (min(vals), max(vals))
        lo, hi = min(lo, *vals), max(hi, *vals)

    # (3) паспортная RMSE аппроксимации модели — независимый вклад, который
    #     ковариация параметров не описывает (переотражения между решётками,
    #     дифракционное рассеяние на больших углах)
    rmse = getattr(passport, "model_rmse_db", 0.0)
    if rmse > 0:
        parts["model_rmse"] = (base - K95 * rmse, base + K95 * rmse)
        lo, hi = lo - K95 * rmse, hi + K95 * rmse

    # (4) штраф за экстраполяцию: расширение вне полосы калибровки
    nu = np.atleast_1d(np.asarray(freqs, dtype=float))
    f_lo, f_hi = passport.band_thz
    out = np.concatenate([nu[nu < f_lo] / max(f_lo, 1e-9), nu[nu > f_hi] / f_hi])
    if out.size:
        # одна октава за полосой -> +1 дБ (эмпирический консервативный штраф)
        oct_dist = float(np.max(np.abs(np.log2(np.maximum(out, 1e-9)))))
        pen = 1.0 * oct_dist
        parts["extrapolation"] = (base - pen, base + pen)
        lo, hi = lo - pen, hi + pen

    return Uncertainty(sigma_db=(hi - lo) / (2 * K95), lo_db=lo, hi_db=hi, parts=parts)


# ---------------------------------------------------------------------------
def slope_db_per_deg(theta_deg, passport: Passport, setup: Setup, freqs, *,
                     weight=None, ref_theta_deg=None, integral=False, h=0.25):
    """Крутизна dA/dθ — прямо определяет достижимую точность установки."""
    a = _eval(theta_deg - h, passport, setup, freqs, weight, ref_theta_deg, integral)
    b = _eval(theta_deg + h, passport, setup, freqs, weight, ref_theta_deg, integral)
    return (b - a) / (2 * h)


def attenuation_floor(passport: Passport, setup: Setup, freqs, *, weight=None,
                      ref_theta_deg=0.0, integral=False, span=(0.0, 180.0), n=721):
    """Пол затухания прибора: максимум A(θ) и угол, где он достигается.

    Ищется ЧИСЛЕННО. Это принципиально: измеренная точка экстинкции лежит не
    на 90 град (у одного из образцов проекта — около 84 град).
    """
    th = np.linspace(span[0], span[1], n)
    if integral:
        a = attenuation_integral_db(th, freqs, weight, passport, setup,
                                    ref_theta1_deg=ref_theta_deg)
    else:
        a, _ = attenuation_db(th, freqs, passport, setup, ref_theta1_deg=ref_theta_deg)
        a = a.mean(axis=1)
    i = int(np.argmax(a))
    return float(a[i]), float(th[i])


def dynamic_range_db(bg_path, freqs_thz, tail_from_thz=3.0):
    """DR(nu) по ВЧ-хвосту фонового измерения (там сигнала физически нет).

    Шумовой пол берётся по СОБСТВЕННОЙ сетке файла, а не по рабочей: рабочая
    сетка обычно 0.2-1.5 ТГц, и хвоста в ней просто нет.
    """
    from pathlib import Path

    raw = np.loadtxt(Path(bg_path))
    t, E = raw[:, 0], raw[:, 1]
    E = E - np.mean(E[:50])
    spec = np.abs(np.fft.rfft(E)) ** 2
    f_full = np.fft.rfftfreq(len(t), float(t[1] - t[0]))

    tail = spec[f_full >= tail_from_thz]
    if tail.size < 8:                       # запасной вариант для коротких сканов
        tail = spec[f_full >= 0.8 * f_full.max()]
    if tail.size < 4 or tail.mean() <= 0:
        return None
    dr_full = 10.0 * np.log10(np.maximum(spec, 1e-300) / tail.mean())
    return np.interp(np.asarray(freqs_thz, dtype=float), f_full, dr_full,
                     left=dr_full[0], right=0.0)


# ---------------------------------------------------------------------------
def run_checks(target_db, theta_deg, passport: Passport, setup: Setup, freqs, *,
               weight=None, ref_theta_deg=None, integral=False,
               tol_sigma_db=1.0, tol_spread_db=1.5, tol_azimuth_deg=2.0,
               dr_db=None, dr_margin_db=6.0,
               floor_db=None, floor_theta=None, achieved_db=None) -> list[Check]:
    """Восемь проверок применимости. Возвращает список Check (ok=False — отказ)."""
    out: list[Check] = []
    nu = np.atleast_1d(np.asarray(freqs, dtype=float))

    if theta_deg is None:
        out.append(Check("F0", False,
                         "обратная задача не имеет решения — угол не определён",
                         "см. F2: цель выше пола прибора"))

    # F1 — усиление невозможно
    out.append(Check("F1", target_db >= -1e-9,
                     f"цель {target_db:+.2f} дБ"
                     + ("" if target_db >= 0 else " — отрицательное затухание невозможно"),
                     "минимум = 0 дБ (опорные углы)"))

    # F2 — пол затухания прибора
    if floor_db is not None:
        ok = target_db <= floor_db + 1e-9
        out.append(Check("F2", ok,
                         f"пол прибора {floor_db:.1f} дБ при theta={floor_theta:.1f} град"
                         + ("" if ok else f"; цель {target_db:.1f} дБ недостижима"),
                         "решётка с меньшим D/P, либо каскад с пластиной HR-Si (−3.01 дБ)"))

    # F3 — динамический диапазон установки
    if dr_db is None:
        out.append(Check("F3", True, "динамический диапазон не задан — проверка пропущена",
                         "загрузить фоновое измерение (bg) для оценки DR(ν)"))
    else:
        dr = float(np.min(np.atleast_1d(dr_db)))
        ok = target_db <= dr - dr_margin_db
        out.append(Check("F3", ok, f"DR установки {dr:.1f} дБ (запас {dr_margin_db:.0f} дБ)"
                         + ("" if ok else "; результат будет НЕПРОВЕРЯЕМ — сигнал уйдёт под шум"),
                         "увеличить усреднение или сузить рабочую полосу"))

    # F4 — угловая разрешимость
    if theta_deg is not None:
        sl = abs(slope_db_per_deg(theta_deg, passport, setup, nu, weight=weight,
                                  ref_theta_deg=ref_theta_deg, integral=integral))
        s_a = sl * passport.scale.sigma_total_deg()
        ok = s_a <= tol_sigma_db
        need = tol_sigma_db / max(sl, 1e-9) * 2.0 / np.sqrt(2.0)
        out.append(Check("F4", ok,
                         f"крутизна {sl:.2f} дБ/град, шкала {passport.scale.division_deg:g} град "
                         f"-> sigma(угол) = {s_a:.2f} дБ (допуск {tol_sigma_db:.2f})",
                         f"нужна шкала с делением <= {need:.2f} град или моторизованный ротатор"))

    # F5 — спектральная неравномерность
    if theta_deg is not None and nu.size > 1:
        db, _ = attenuation_db(theta_deg, nu, passport, setup, ref_theta1_deg=ref_theta_deg)
        m = nu <= passport.zone_red_min_thz
        spread = float(db[0, m].max() - db[0, m].min()) if m.any() else 0.0
        ok = spread <= tol_spread_db
        out.append(Check("F5", ok,
                         f"неравномерность {spread:.2f} дБ p-p в {nu[m].min():.2f}–{nu[m].max():.2f} ТГц"
                         if m.any() else "полоса вне области применимости",
                         "сузить рабочую полосу или использовать пластины HR-Si"))

    # F6 — зона применимости решётки
    z = passport.zone(nu)
    n_red, n_amb = int((z == "red").sum()), int((z == "amber").sum())
    out.append(Check("F6", n_red == 0,
                     f"зоны: зелёная {(z=='green').sum()}, жёлтая {n_amb}, красная {n_red} "
                     f"(nu аномалии = {passport.nu_anomaly_thz:.1f} ТГц, "
                     f"зелёная до {passport.zone_green_max_thz:.2f} ТГц)",
                     "снизить верхнюю частоту: в красной зоне решётка работает как дифракционная"))

    # F7 — неопределённость экстраполяции
    if theta_deg is not None:
        u = uncertainty(theta_deg, passport, setup, nu, weight=weight,
                        ref_theta_deg=ref_theta_deg, integral=integral)
        ok = u.sigma_db <= tol_sigma_db
        f_lo, f_hi = passport.band_thz
        n_out = int(((nu < f_lo) | (nu > f_hi)).sum())
        out.append(Check("F7", ok,
                         f"полная σ = {u.sigma_db:.2f} дБ, 95 % [{u.lo_db:.2f}, {u.hi_db:.2f}] дБ"
                         + (f"; {n_out} частот вне полосы калибровки {f_lo}–{f_hi} ТГц" if n_out else ""),
                         "докалибровать прибор в нужной полосе"))

    # F8 — состояние поляризации на выходе
    if theta_deg is not None:
        az, el = output_polarization(theta_deg, nu, passport, setup)
        m = z != "red"
        if m.any():
            az_max = float(np.max(np.abs(az[0, m] - setup.det_deg)))
            el_max = float(np.max(np.abs(el[0, m])))
            ok = az_max <= tol_azimuth_deg
            i_bad = int(np.argmax(np.abs(az[0] - setup.det_deg) * m))
            out.append(Check("F8", ok,
                             f"азимут выхода отклонён на {az_max:.2f} град "
                             f"(на {nu[i_bad]:.2f} ТГц), эллиптичность {el_max:.2f} град",
                             "в глубокой тени основная компонента подавлена сильнее, чем "
                             "утечка выходного поляризатора, поэтому азимут уезжает; "
                             "не уходить так глубоко либо ставить выходной анализатор "
                             "с лучшей экстинкцией в этой части полосы"))
    return out
