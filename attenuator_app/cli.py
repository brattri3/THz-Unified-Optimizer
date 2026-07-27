"""Интерактивный CLI ТГц-аттенюатора (бесконечный цикл).

Запуск из корня THz-Unified-Optimizer:
    .venv\\Scripts\\python.exe -m attenuator_app.cli
    .venv\\Scripts\\python.exe -m attenuator_app.cli --passport attenuator_app/passports/<файл>.json

Все затухания — в дБ ПО МОЩНОСТИ (10*log10). Отношение по полю показывается
справочно рядом. Углы вводятся в градусах; прибор с ручной установкой угла,
поэтому решения привязываются к делениям шкалы из паспорта.
"""
from __future__ import annotations

import argparse
import shlex
import sys
import traceback
from pathlib import Path

import numpy as np

from .api import Attenuator
from .core.forward import DETECTORS, SCHEMES, SOURCES, ideal_cos4_db
from .core.passport import Passport
from .core import plots
from .core.limits import uncertainty
from .core.session import Session, ascii_plot
from .core.weights import PRESETS

HERE = Path(__file__).resolve().parent
DEFAULT_PASSPORT = HERE / "passports" / "ATT-11-16-CA85_02721.json"

HELP = """
Команды
  dev                          сводка прибора и текущей конфигурации
  scheme [S0|S1|S2|S3]         оптическая схема (без аргумента - список)
  det    [coherent|power]      тип детектора
  src    [linear|unpolarized]  тип источника
  psi <град> | dop <0..1>      азимут и степень поляризации источника
  band <f1> <f2> [N]           рабочая сетка частот, ТГц
  w <тип> [опции]              спектральный вес: preset <имя> | flat <f1> <f2>
                               | gauss <nu0> <dnu> | bg <файл>
  bg <файл>                    фоновое измерение -> динамический диапазон (проверка F3)

  zero [град]                  SET ZERO: объявить положение опорным (0 дБ)
  autozero                     AUTO-ZERO: найти максимум пропускания
  autocross                    AUTO-CROSS: найти точку экстинкции (она НЕ на 90 град)
  set <град>                   считать, что прибор сейчас на этом угле
  mode abs|rel                 режим отсчёта затухания

  a <град>                     прямая задача: затухание в точке
  ai <град>                    то же, интегрально по спектру источника
  spec <град>                  спектр затухания + зоны + поляризация выхода
  curve [шаг]                  таблица A(theta) рядом с идеальной cos^4
  solve <дБ> [-i]              обратная задача (-i = интегрально)
  check <дБ> [град]            восемь проверок применимости F1..F8
  sweep <дБ,дБ,...>            развёртка уровней, монотонный обход
  plan <дБ>                    решение + инструкция оператору
  pol <град>                   поляризация на выходе (азимут, эллиптичность)

  ascii on|off                 печатать ли ASCII-графики в терминал
  help | q
"""


def bar(v, lo, hi, width=34, ch="#"):
    if not np.isfinite(v):
        return ""
    x = 0.0 if hi <= lo else (v - lo) / (hi - lo)
    return ch * int(np.clip(x, 0, 1) * width)


class Shell:
    def __init__(self, att: Attenuator, sess=None):
        self.att = att
        self.relative = True
        self.sess = sess
        self.ascii_plots = True
        self.sync_choices()

    def sync_choices(self):
        """Записать текущее состояние выборов в метафайл сеанса."""
        if not self.sess:
            return
        a, s, p = self.att, self.att.setup, self.att.passport
        self.sess.choose(
            passport=f"{p.serial}", scheme=s.scheme, detector=s.detector,
            source=s.source, psi_deg=s.psi_deg, dop=s.dop, det_deg=s.det_deg,
            phi2_deg=s.phi2_deg, use_film=s.use_film,
            mode="relative" if self.relative else "absolute",
            zero_deg=a.zero_deg, theta_deg=a.theta_deg,
            band_thz=[float(a.freqs[0]), float(a.freqs[-1]), len(a.freqs)],
            weight=a.weight_desc, scale_division_deg=p.scale.division_deg,
            sigma_angle_deg=round(p.scale.sigma_total_deg(), 3),
            model_rmse_db=p.model_rmse_db)

    # -- вывод ---------------------------------------------------------
    def dev(self):
        s = self.att.device_summary()
        p = self.att.passport
        print(f"\nПрибор {s['serial']}  {s['model']}")
        print(f"  апертура     {s['aperture_mm']:g} мм -> шкала {s['scale_division_deg']:g} град, "
              f"sigma(угол) {s['sigma_angle_deg']:.2f} град")
        print(f"  P            {s['P_um']:.3f} мкм")
        print(f"  D_eff        {p.D_eff_um.value:.3f} +- {p.D_eff_um.sigma:.3f} мкм "
              f"(D_eff/D_phys = {s['D_eff_over_D_phys']:.3f}; это ЭФФЕКТИВНЫЙ, не геометрический)")
        print(f"  потери       {p.loss_factor.value:.3f} +- {p.loss_factor.sigma:.3f} дБ/ТГц^{p.gamma.value:g}"
              f"  (диапазон gamma {p.gamma_range[0]}..{p.gamma_range[1]} -> огибающая)")
        print(f"  offset       {p.angle_offset_deg.value:+.2f} +- {p.angle_offset_deg.sigma:.2f} град")
        print(f"  полоса калибровки {p.band_thz[0]}-{p.band_thz[1]} ТГц, схема {p.fit_scheme}, "
              f"детектор {p.fit_detector}")
        print(f"  nu аномалии  {s['nu_anomaly_thz']:.1f} ТГц; зелёная зона до {s['zone_green_max_thz']:.2f} ТГц")
        if p.film:
            print(f"  плёнки       {p.film.name}: K1={p.film.K1:.2f}, "
                  f"ER {p.film.er_db(0.5):.0f} дБ @0.5 ТГц / {p.film.er_db(2.0):.0f} дБ @2 ТГц")
        print(f"\nКонфигурация: схема {self.att.setup.scheme} ({SCHEMES[self.att.setup.scheme]})")
        print(f"  детектор {self.att.setup.detector}, источник {self.att.setup.source} "
              f"(psi={self.att.setup.psi_deg:g} град, DOP={self.att.setup.dop:g})")
        print(f"  сетка {self.att.freqs[0]:.2f}-{self.att.freqs[-1]:.2f} ТГц, {len(self.att.freqs)} точек")
        print(f"  вес источника: {self.att.weight_desc}")
        print(f"  режим: {'относительный' if self.relative else 'абсолютный'}; "
              f"ноль: {'не задан' if self.att.zero_deg is None else f'{self.att.zero_deg:+.2f} град'}; "
              f"текущий угол {self.att.theta_deg:+.2f} град")
        if p.calibration_note:
            print(f"\n  {p.calibration_note}")

    def show_att(self, theta, integral=False):
        r = self.att.attenuation(theta, relative=self.relative, integral=integral)
        kind = "интегральное" if integral else "среднее по полосе"
        print(f"\n  theta = {r['theta_deg']:+.2f} град   ({kind}, "
              f"{'относительно нуля' if r['relative'] else 'абсолютное'})")
        print(f"  ЗАТУХАНИЕ  {r['att_db']:.2f} дБ   "
              f"[95%: {r['ci95_db'][0]:.2f} ... {r['ci95_db'][1]:.2f}]  sigma={r['sigma_db']:.2f} дБ")
        print(f"  по мощности x{r['power_ratio']:.4g}   по полю x{r['field_ratio']:.4g}   "
              f"идеал cos^4: {r['ideal_cos4_db']:.2f} дБ")
        if r["parts"]:
            print("  вклад в неопределённость (95% границы, дБ):")
            for k, (lo, hi) in r["parts"].items():
                print(f"     {k:<14} {lo:7.2f} ... {hi:7.2f}")

    def spec(self, theta):
        s = self.att.spectrum(theta, relative=self.relative)
        db, nu, z = s["att_db"], s["freq_thz"], s["zone"]
        m = z != "red"
        lo, hi = (np.nanmin(db[m]), np.nanmax(db[m])) if m.any() else (0, 1)
        title = (f"спектр затухания при theta = {theta:+.2f} град "
                 f"({'относительный' if self.relative else 'абсолютный'} режим)")
        lines = [f"  {'ТГц':>6} {'дБ':>8} {'зона':>6} {'экстр':>6} {'азимут':>7} {'эллипт':>7}"]
        step = max(1, len(nu) // 24)
        for i in range(0, len(nu), step):
            lines.append(f"  {nu[i]:6.2f} {db[i]:8.2f} {str(z[i]):>6} "
                         f"{'да' if s['extrapolated'][i] else '  ':>6} "
                         f"{s['azimuth_deg'][i]:7.2f} {s['ellipticity_deg'][i]:7.2f}  "
                         f"{bar(db[i], lo, hi)}")
        if m.any():
            lines.append(f"  размах в применимой области: "
                         f"{db[m].max()-db[m].min():.2f} дБ p-p")
        n_ex = int(s["extrapolated"].sum())
        if n_ex:
            lines.append(f"  ВНИМАНИЕ: {n_ex} точек вне полосы калибровки "
                         f"{self.att.passport.band_thz[0]}-{self.att.passport.band_thz[1]} ТГц")
        body = "\n".join(lines)
        print(f"\n  {title}\n{body}")

        # коридор 95 % по частотам — для заливки на графике
        lo_db, hi_db = [], []
        for f in nu:
            u = uncertainty(theta, self.att.passport, self.att.setup, np.array([f]),
                            ref_theta_deg=self.att.zero_deg if self.relative else None)
            lo_db.append(u.lo_db)
            hi_db.append(u.hi_db)
        dr = None
        if self.att._dr is not None:
            dr = np.atleast_1d(self.att._dr)
            dr = dr if len(dr) == len(nu) else None

        self.emit(title, body, "spectrum",
                  ascii_plot(nu, db, xlabel="частота, ТГц", ylabel="затухание, дБ"),
                  {"freq_thz": nu, "att_db": db, "lo_db": lo_db, "hi_db": hi_db,
                   "zone": [str(v) for v in z],
                   "extrapolated": s["extrapolated"], "azimuth_deg": s["azimuth_deg"],
                   "ellipticity_deg": s["ellipticity_deg"]},
                  png=lambda path: plots.spectrum(
                      path, nu, db, lo_db=lo_db, hi_db=hi_db, zones=z,
                      band=self.att.passport.band_thz, extrapolated=s["extrapolated"],
                      dr_db=dr, theta_deg=theta,
                      mode="относительный" if self.relative else "абсолютный",
                      footer=self.footer()))

    def curve(self, step=10.0):
        from .core.limits import slope_db_per_deg
        ref = self.att.zero_deg or 0.0
        th = np.arange(ref, ref + 90.0 + 1e-9, step)
        cols = {"theta_deg": [], "att_db": [], "lo_db": [], "hi_db": [],
                "ideal_cos4_db": [], "slope_db_per_deg": []}
        lines = [f"  {'theta':>7} {'модель,дБ':>10} {'95% интервал':>18} "
                 f"{'идеал cos^4':>12} {'крутизна':>10}"]
        for t in th:
            r = self.att.attenuation(t, relative=True)
            ideal = float(ideal_cos4_db(t, ref))
            sl = slope_db_per_deg(t, self.att.passport, self.att.setup, self.att.freqs,
                                  ref_theta_deg=ref)
            lines.append(f"  {t:7.1f} {r['att_db']:10.2f} "
                         f"[{r['ci95_db'][0]:7.2f},{r['ci95_db'][1]:7.2f}] "
                         f"{ideal:12.2f} {sl:9.2f}/град")
            for k, v in zip(cols, (t, r["att_db"], r["ci95_db"][0], r["ci95_db"][1],
                                   ideal, sl)):
                cols[k].append(v)
        body = "\n".join(lines)
        print("\n" + body)

        # Таблица печатается с шагом пользователя, а кривая строится на плотной
        # сетке: при шаге 10 град график получается ломаным и врёт про форму.
        fine = np.arange(ref, ref + 90.0 + 1e-9, 2.0)
        f = {"theta_deg": fine, "att_db": [], "lo_db": [], "hi_db": [],
             "ideal_cos4_db": [], "slope_db_per_deg": []}
        for t in fine:
            r = self.att.attenuation(t, relative=True)
            f["att_db"].append(r["att_db"])
            f["lo_db"].append(r["ci95_db"][0])
            f["hi_db"].append(r["ci95_db"][1])
            f["ideal_cos4_db"].append(float(ideal_cos4_db(t, ref)))
            f["slope_db_per_deg"].append(
                slope_db_per_deg(t, self.att.passport, self.att.setup,
                                 self.att.freqs, ref_theta_deg=ref))

        self.emit("угловая кривая затухания", body, "curve",
                  ascii_plot(cols["theta_deg"], cols["att_db"],
                             xlabel="угол, град", ylabel="затухание, дБ",
                             ref=cols["ideal_cos4_db"], ref_label="идеал cos^4"),
                  cols,
                  png=lambda path: plots.angular_curve(
                      path, f["theta_deg"], f["att_db"], f["ideal_cos4_db"],
                      lo_db=f["lo_db"], hi_db=f["hi_db"],
                      slope=f["slope_db_per_deg"], footer=self.footer()))

    def emit(self, title, body, slug, plot=None, cols=None, png=None):
        """Показать результат и сохранить его в папку сеанса.

        ASCII-график печатается в терминал всегда (там его и читают), а
        matplotlib-версия — если библиотека доступна: она нужна для отчётов и
        приложений к паспорту, где псевдографика неуместна.
        """
        if plot and self.ascii_plots:
            print("\n" + plot)
        if not self.sess:
            return
        self.sync_choices()
        p = self.sess.save_plot(title, (body + "\n\n" + plot) if plot else body, slug)
        if cols:
            self.sess.save_csv(slug, cols)
        saved = [f"plots/{p.name}"]
        if png is not None and plots.available():
            try:
                fig_dir = self.sess.dir / "figures"
                fig_dir.mkdir(exist_ok=True)
                out = png(fig_dir / f"{p.stem}.png")
                saved.append(f"figures/{Path(out).name}")
                self.sess.meta["artifacts"].append(
                    {"file": f"figures/{Path(out).name}", "title": title, "kind": "png"})
            except Exception as e:
                print(f"  [график не построен: {type(e).__name__}: {e}]")
        print(f"\n  сохранено: {self.sess.dir.name}/{{{', '.join(saved)}}}")

    def footer(self):
        """Подпись под графиком: без неё число в дБ не интерпретируется."""
        p, s, a = self.att.passport, self.att.setup, self.att
        return (f"{p.serial} · схема {s.scheme} · детектор {s.detector} · "
                f"шкала {p.scale.division_deg:g}° (sigma {p.scale.sigma_total_deg():.2f}°) · "
                f"полоса калибровки {p.band_thz[0]}–{p.band_thz[1]} ТГц · "
                f"ноль {a.zero_deg if a.zero_deg is not None else '—'}° · "
                f"вес {a.weight_desc} · сеанс {self.sess.dir.name if self.sess else '—'}")

    def pol(self, theta):
        """Состояние поляризации на выходе: азимут и эллиптичность по спектру."""
        s = self.att.spectrum(theta, relative=self.relative)
        nu, az, el = s["freq_thz"], s["azimuth_deg"], s["ellipticity_deg"]
        lines = [f"  {'ТГц':>6} {'азимут':>9} {'эллипт':>9} {'зона':>7}"]
        step = max(1, len(nu) // 20)
        for i in range(0, len(nu), step):
            lines.append(f"  {nu[i]:6.2f} {az[i]:9.2f} {el[i]:9.2f} {str(s['zone'][i]):>7}")
        lines.append(f"  максимум |азимут| = {np.max(np.abs(az)):.2f} град, "
                     f"|эллиптичность| = {np.max(np.abs(el)):.2f} град")
        body = "\n".join(lines)
        print("\n  поляризация на выходе\n" + body)
        self.emit("поляризация на выходе", body, "polarization",
                  ascii_plot(nu, az, xlabel="частота, ТГц", ylabel="азимут, град"),
                  {"freq_thz": nu, "azimuth_deg": az, "ellipticity_deg": el,
                   "zone": [str(v) for v in s["zone"]]},
                  png=lambda path: plots.polarization(
                      path, nu, az, el, zones=s["zone"],
                      band=self.att.passport.band_thz, theta_deg=theta,
                      footer=self.footer()))

    def solve(self, target, integral=False):
        sols = self.att.solve(target, integral=integral)
        if not sols:
            print(f"\n  ЦЕЛЬ {target:.2f} дБ НЕДОСТИЖИМА")
            th_ext, floor = self.att.auto_cross()
            print(f"  максимум прибора: {floor:.2f} дБ при theta = {th_ext:.1f} град")
            print("  средства: пластина HR-Si (-3.01 дБ, широкополосно); "
                  "решётка с меньшим D/P; сузить полосу")
            return
        print(f"\n  решения для цели {target:.2f} дБ "
              f"({'интегрально' if integral else 'по полосе'}, "
              f"шкала {self.att.passport.scale.division_deg:g} град):")
        for s in sols:
            flag = "" if abs(s.achieved_db - target) < 0.5 else "  <- шкала не даёт точнее"
            print(f"   {s.rank}. {s}{flag}")
        best = sols[0]
        print(f"\n  РЕКОМЕНДАЦИЯ: выставить {best.theta_set_deg:+.2f} град "
              f"-> {best.achieved_db:.2f} дБ [95%: {best.lo_db:.2f} ... {best.hi_db:.2f}]")

    def check(self, target, theta=None):
        if theta is None:
            s = self.att.solve(target)
            theta = s[0].theta_set_deg if s else None
            if theta is not None:
                print(f"  (угол взят из решения обратной задачи: {theta:+.2f} град)")
        print()
        checks = self.att.checks(target, theta)
        for c in checks:
            print("  " + str(c).replace("\n", "\n  "))
        bad = [c.code for c in checks if not c.ok]
        print(f"\n  ИТОГ: {'ДОСТИЖИМО' if not bad else 'ОТКАЗ по ' + ', '.join(bad)}")

    def sweep(self, targets):
        rows = self.att.sweep(targets)
        print(f"\n  развёртка уровней (обход по возрастанию угла, без реверсов)")
        print(f"  {'цель,дБ':>8} {'theta':>8} {'факт,дБ':>9} {'95% интервал':>20} {'крутизна':>10}")
        for t, s in rows:
            if s is None:
                print(f"  {t:8.1f} {'--':>8} {'НЕДОСТИЖИМО':>9}")
            else:
                print(f"  {t:8.1f} {s.theta_set_deg:8.2f} {s.achieved_db:9.2f} "
                      f"  [{s.lo_db:7.2f} ...{s.hi_db:7.2f}] {s.slope_db_per_deg:9.2f}")

    def plan(self, target):
        sols = self.att.solve(target)
        if not sols:
            print("  недостижимо — см. solve")
            return
        s = sols[0]
        print(f"\n  ПЛАН для {target:.2f} дБ:")
        for c in self.att.motion_plan(s):
            print("   " + c.as_operator_instruction())
        print(f"   ожидаемое затухание {s.achieved_db:.2f} дБ "
              f"[95%: {s.lo_db:.2f} ... {s.hi_db:.2f}]")
        self.att.apply(s)
        print(f"   текущий угол теперь {self.att.theta_deg:+.2f} град")

    # -- цикл ----------------------------------------------------------
    def run(self):
        print(f"THz Attenuator CLI  (паспорт {self.att.passport.serial}).  "
              f"help - список команд, q - выход")
        while True:
            try:
                line = input("\natt> ").strip()
            except (EOFError, KeyboardInterrupt):
                print()
                return
            if not line:
                continue
            if self.sess:
                self.sess.log(line)
            try:
                res = self.dispatch(shlex.split(line))
                self.sync_choices()
                if res is False:
                    if self.sess:
                        self.sess.close()
                        print(f"  сеанс сохранён: {self.sess.dir}")
                    return
            except Exception as e:
                print(f"  ОШИБКА: {type(e).__name__}: {e}")
                if "--debug" in sys.argv:
                    traceback.print_exc()

    def dispatch(self, argv):
        c, a = argv[0].lower(), argv[1:]
        A, S = self.att, self.att.setup

        if c in ("q", "quit", "exit"):
            return False
        elif c in ("help", "?", "h"):
            print(HELP)
        elif c == "dev":
            self.dev()
        elif c == "scheme":
            if a:
                S.scheme = a[0].upper(); S.validate(); print(f"  схема {S.scheme}: {SCHEMES[S.scheme]}")
            else:
                for k, v in SCHEMES.items():
                    print(f"  {k}: {v}")
        elif c == "det":
            if a:
                S.detector = a[0].lower(); S.validate(); print(f"  {DETECTORS[S.detector]}")
            else:
                for k, v in DETECTORS.items():
                    print(f"  {k}: {v}")
        elif c == "src":
            if a:
                S.source = a[0].lower(); S.validate(); print(f"  {SOURCES[S.source]}")
            else:
                for k, v in SOURCES.items():
                    print(f"  {k}: {v}")
        elif c == "psi":
            S.psi_deg = float(a[0]); print(f"  psi = {S.psi_deg:g} град")
        elif c == "dop":
            S.dop = float(a[0]); S.validate(); print(f"  DOP = {S.dop:g}")
        elif c == "band":
            n = int(a[2]) if len(a) > 2 else 128
            A.freqs = np.linspace(float(a[0]), float(a[1]), n)
            A.weight = None; A.weight_desc = "не задан (сетка изменилась)"
            print(f"  сетка {A.freqs[0]:.3f}-{A.freqs[-1]:.3f} ТГц, {n} точек")
        elif c == "w":
            kind = a[0].lower()
            if kind == "preset":
                lo, hi = A.set_weight("preset", name=a[1] if len(a) > 1 else "pca")
            elif kind == "flat":
                lo, hi = A.set_weight("flat", f1=float(a[1]), f2=float(a[2]))
            elif kind == "gauss":
                lo, hi = A.set_weight("gauss", nu0=float(a[1]), dnu=float(a[2]))
            elif kind == "bg":
                lo, hi = A.set_weight("bg_file", path=a[1])
            else:
                print(f"  пресеты: {list(PRESETS)}"); return
            print(f"  вес: {A.weight_desc}; эффективная полоса {lo:.2f}-{hi:.2f} ТГц")
        elif c == "bg":
            dr = A.load_background(a[0])
            print("  DR не определён (нет ВЧ-хвоста)" if dr is None
                  else f"  динамический диапазон установки: min {dr:.1f} дБ по сетке")
        elif c == "zero":
            z = A.set_zero(float(a[0]) if a else None)
            print(f"  SET ZERO: опорный угол {z:+.2f} град, это 0 дБ по определению")
        elif c == "autozero":
            z = A.auto_zero(); print(f"  AUTO-ZERO: максимум пропускания на {z:+.2f} град")
        elif c == "autocross":
            th, fl = A.auto_cross()
            d = abs(th - (A.zero_deg or 0.0))
            print(f"  AUTO-CROSS: экстинкция на {th:.2f} град, пол {fl:.2f} дБ")
            print(f"  расстояние от нуля {d:.2f} град"
                  + ("" if abs(d - 90) <= 2 else
                     f"  <- отличается от 90 на {d-90:+.1f} град: перекос либо не-Малюсовость"))
        elif c == "set":
            A.theta_deg = float(a[0]); print(f"  текущий угол {A.theta_deg:+.2f} град")
        elif c == "mode":
            self.relative = a[0].lower().startswith("rel")
            print(f"  режим: {'относительный' if self.relative else 'абсолютный'}")
        elif c == "a":
            self.show_att(float(a[0]) if a else A.theta_deg)
        elif c == "ai":
            self.show_att(float(a[0]) if a else A.theta_deg, integral=True)
        elif c == "spec":
            self.spec(float(a[0]) if a else A.theta_deg)
        elif c == "ascii":
            self.ascii_plots = not a or a[0].lower() in ("on", "1", "да")
            print(f"  ASCII-графики в терминале: {'вкл' if self.ascii_plots else 'выкл'}"
                  + ("" if plots.available() else "  (matplotlib недоступен!)"))
        elif c == "pol":
            self.pol(float(a[0]) if a else A.theta_deg)
        elif c == "curve":
            self.curve(float(a[0]) if a else 10.0)
        elif c == "solve":
            self.solve(float(a[0]), integral=("-i" in a))
        elif c == "check":
            self.check(float(a[0]), float(a[1]) if len(a) > 1 else None)
        elif c == "sweep":
            self.sweep([float(x) for x in ",".join(a).replace(";", ",").split(",") if x])
        elif c == "plan":
            self.plan(float(a[0]))
        else:
            print(f"  неизвестная команда {c!r}; help - список")
        return True


def main():
    ap = argparse.ArgumentParser(description="ТГц аттенюатор — интерактивный расчёт")
    ap.add_argument("--passport", default=str(DEFAULT_PASSPORT))
    ap.add_argument("--scheme", default=None)
    ap.add_argument("--detector", default=None)
    ap.add_argument("--runs", default=None, help="корень папок сеансов (умолч. ./runs)")
    ap.add_argument("--tag", default="", help="суффикс имени папки сеанса")
    ap.add_argument("--no-session", action="store_true", help="не вести папку сеанса")
    ap.add_argument("--debug", action="store_true")
    args = ap.parse_args()

    p = Passport.load(args.passport)
    att = Attenuator(p)
    if args.scheme:
        att.setup.scheme = args.scheme.upper()
    if args.detector:
        att.setup.detector = args.detector
    att.setup.validate()
    att.set_zero(0.0)

    sess = None if args.no_session else Session(args.runs, args.tag)
    if sess:
        print(f"сеанс: {sess.dir}")
    Shell(att, sess).run()


if __name__ == "__main__":
    main()
