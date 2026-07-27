"""Валидация модели на реальных измерениях: покрытие 95 % интервала.

Вопрос, на который отвечает скрипт: если оператор выставит угол theta, попадёт
ли ИЗМЕРЕННОЕ относительное затухание в интервал, который предсказала программа,
с заявленной доверительной вероятностью 95 %?

Метод — относительный режим (SET ZERO на theta = 0):
    A_изм(theta) = -10 log10 [ U_sig(theta)/U_bg(theta) : U_sig(0)/U_bg(0) ]
где U = интеграл |E(nu)|^2 по рабочей полосе (Парсеваль). Деление на фон снимает
дрейф лазера между точками — иначе он маскирует проверку.

Данные: один и тот же прибор, шесть сеансов.
    356att      — вращался ВТОРОЙ поляризатор (ближний к детектору, анализатор)
    series1..5  — вращался ПЕРВЫЙ поляризатор
По закону Малюса затухание зависит только от относительного угла, поэтому оба
набора должны ложиться на одну кривую. Расхождение между группами — прямая
проверка этого утверждения на железе.

Запуск:
    .venv\\Scripts\\python.exe -m attenuator_app.tools.validate
    ... --passport <файл> --band 0.2 1.5 --datasets 356att,series1
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
sys.path.insert(0, str(REPO))

from attenuator_app.api import Attenuator
from attenuator_app.core.forward import Setup
from attenuator_app.core.limits import uncertainty
from attenuator_app.core.passport import Passport

DATA = REPO / "data_pool"
ROTATES = {"356att": "второй (анализатор)", "series1": "первый", "series2": "первый",
           "series3": "первый", "series4": "первый", "series5": "первый"}


def energy(path, f1, f2):
    """Энергия импульса в полосе [f1, f2] по Парсевалю."""
    raw = np.loadtxt(path)
    t, E = raw[:, 0], raw[:, 1]
    E = E - np.mean(E[:50])
    spec = np.abs(np.fft.rfft(E)) ** 2
    f = np.fft.rfftfreq(len(t), float(t[1] - t[0]))
    m = (f >= f1) & (f <= f2)
    return float(np.trapezoid(spec[m], f[m]))


def measured(dataset, f1, f2):
    """{угол: нормированное на фон пропускание} для одного сеанса."""
    out = {}
    for sig in sorted(DATA.glob(f"{dataset}_*deg_rep*_sig.txt")):
        bg = Path(str(sig).replace("_sig.txt", "_bg.txt"))
        if not bg.exists():
            continue
        try:
            us, ub = energy(sig, f1, f2), energy(bg, f1, f2)
        except Exception:
            continue
        if ub <= 0:
            continue
        ang = float(sig.name.split("_")[1].replace("deg", ""))
        ang = ((ang + 180.0) % 360.0) - 180.0            # нормализация как в DataManager
        out.setdefault(ang, []).append(us / ub)
    return {a: float(np.mean(v)) for a, v in out.items()}


def run(passport_path, band, datasets, verbose=True):
    p = Passport.load(passport_path)
    f1, f2 = band
    nu = np.linspace(f1, f2, 96)
    att = Attenuator(p, Setup(scheme=p.fit_scheme, detector=p.fit_detector,
                              use_film=False), freqs_thz=nu)
    att.set_zero(0.0)
    att.set_weight("flat", f1=f1, f2=f2)

    print(f"=== валидация {p.serial} ===")
    print(f"полоса {f1}-{f2} ТГц, схема {att.setup.scheme}, детектор {att.setup.detector}, "
          f"шкала {p.scale.division_deg:g} град")
    print(f"относительный режим, SET ZERO на 0 град\n")

    rows, tot_in, tot_n = [], 0, 0
    for ds in datasets:
        meas = measured(ds, f1, f2)
        if 0.0 not in meas or len(meas) < 3:
            print(f"  {ds}: пропущен (нет опорной точки 0 град или мало углов)")
            continue
        ref = meas[0.0]
        angles = sorted(meas)
        n_in = 0
        if verbose:
            print(f"  --- {ds}  (вращался {ROTATES.get(ds, '?')}, {len(angles)} углов) ---")
            print(f"  {'theta':>7} {'изм,дБ':>9} {'модель,дБ':>10} {'95% интервал':>19} "
                  f"{'откл':>7}  ")
        for a in angles:
            a_meas = -10.0 * np.log10(max(meas[a] / ref, 1e-30))
            pred = att.attenuation(a, relative=True, integral=True)
            u = uncertainty(a, p, att.setup, nu, weight=att.weight,
                            ref_theta_deg=0.0, integral=True)
            inside = u.lo_db - 1e-9 <= a_meas <= u.hi_db + 1e-9
            n_in += inside
            tot_n += 1
            rows.append((ds, a, a_meas, pred["att_db"], u.lo_db, u.hi_db, inside))
            if verbose:
                print(f"  {a:7.1f} {a_meas:9.2f} {pred['att_db']:10.2f} "
                      f"[{u.lo_db:7.2f},{u.hi_db:7.2f}] {a_meas-pred['att_db']:+7.2f}"
                      f"  {'' if inside else '<-- ВНЕ ИНТЕРВАЛА'}")
        tot_in += n_in
        cov = 100.0 * n_in / len(angles)
        print(f"  {ds}: покрытие {n_in}/{len(angles)} = {cov:.0f} %\n")

    if not tot_n:
        print("нет данных для валидации")
        return 1

    cov = 100.0 * tot_in / tot_n
    dev = np.array([r[2] - r[3] for r in rows])
    print("=" * 62)
    print(f"ИТОГО покрытие 95 % интервала: {tot_in}/{tot_n} = {cov:.1f} %")
    print(f"смещение (изм - модель): среднее {dev.mean():+.2f} дБ, RMSE {np.sqrt((dev**2).mean()):.2f} дБ")

    # разбивка по тому, какой поляризатор вращался — проверка закона Малюса на железе
    groups = {}
    for grp, name in (("356att", "второй (анализатор)"), ("series", "первый")):
        sel = [r for r in rows if r[0].startswith(grp)]
        if not sel:
            continue
        d = np.array([r[2] - r[3] for r in sel])
        n_in_g = sum(r[6] for r in sel)
        groups[grp] = d
        print(f"  вращался {name:<20}: {n_in_g}/{len(sel)} = {100*n_in_g/len(sel):3.0f} %, "
              f"смещение {d.mean():+.2f} дБ, RMSE {np.sqrt((d**2).mean()):.2f} дБ")
    if len(groups) == 2:
        a, b = groups["356att"], groups["series"]
        print(f"  -> разница смещений между группами {a.mean()-b.mean():+.2f} дБ "
              f"(закон Малюса: затухание не должно зависеть от того, "
              f"какой поляризатор вращать)")

    # насколько интервал консервативен: во сколько раз он шире фактического разброса
    half = np.array([(r[5] - r[4]) / 2 for r in rows])
    need = 1.959963985 * float(np.sqrt((dev ** 2).mean()))
    print(f"  ширина интервала: полу-ширина в среднем {half.mean():.2f} дБ; "
          f"по фактическому разбросу хватило бы {need:.2f} дБ "
          f"-> запас x{half.mean()/max(need,1e-9):.2f}")

    # где именно промахи
    bad = [r for r in rows if not r[6]]
    if bad:
        print(f"\n  вне интервала {len(bad)} точек, углы: "
              + ", ".join(f"{r[1]:.0f}" for r in sorted(bad, key=lambda x: abs(x[1]))[:14]))
        big = max(bad, key=lambda r: abs(r[2] - r[3]))
        print(f"  худшая: {big[0]} {big[1]:.0f} град, изм {big[2]:.2f} vs модель {big[3]:.2f} дБ")
    print(f"\n  ВЕРДИКТ: {'заявленные 95 % подтверждаются' if cov >= 95 else 'покрытие НИЖЕ заявленных 95 % — интервал занижен'}")

    from attenuator_app.core import plots
    if plots.available():
        out = HERE.parent / "passports" / "validation.png"
        try:
            plots.validation(out, rows, footer=(
                f"{p.serial} · схема {att.setup.scheme} · детектор {att.setup.detector} · "
                f"полоса {f1}–{f2} ТГц · относительный режим, SET ZERO на 0° · "
                f"покрытие {tot_in}/{tot_n} = {cov:.1f} %"))
            print(f"  график: {out}")
        except Exception as e:
            print(f"  [график не построен: {type(e).__name__}: {e}]")
    return rows


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--passport",
                    default=str(HERE.parent / "passports" / "ATT-11-16-CA85_02721.json"))
    ap.add_argument("--band", nargs=2, type=float, default=[0.2, 1.5])
    ap.add_argument("--datasets", default="356att,series1,series2,series3,series4,series5")
    ap.add_argument("--quiet", action="store_true")
    a = ap.parse_args()
    run(a.passport, a.band, [d for d in a.datasets.split(",") if d], verbose=not a.quiet)


if __name__ == "__main__":
    main()
