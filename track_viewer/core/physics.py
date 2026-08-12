# -*- coding: utf-8 -*-
"""Расчёты пропускания: интеграл во времени, спектр, окно, пик, выбор референса.

Главная величина — интегральное пропускание во временной области (ТЗ §3.1):

    T = ∫ (E_s(t) − DC_s)² dt  /  ∫ (E_r(t) − DC_r)² dt

Почему именно она, а не медиана по полосе: по равенству Парсеваля интеграл
квадрата поля по времени равен интегралу спектральной плотности по частоте, так
что выбирать границы полосы усреднения не требуется вовсе. У величины нет
свободного параметра, о котором можно было бы спорить постфактум.

Сверка с ядром
--------------
При `dc_mode = 'none'` формула совпадает с `T_time` из
`research/two_wgp/a21_transmission_summary.py`, и на `003_a0.txt`/`002_bg.txt`
даёт **0.98909** — число из ТЗ §7.2. Режим `'pre_pulse'` (по умолчанию, решение
владельца 2026-08-12) даёт **0.99422**: вычитание постоянного смещения детектора
поднимает `T` на +0.51 п.п. Это больше разброса между повторами, поэтому оба
числа инструмент показывает рядом, а не выбирает одно молча.

Отличие от `a21_transmission.json` (99.08 %) — не в этой формуле: там счёт идёт
через `DataManager`, который при двух повторах угла **когерентно усредняет сами
трассы** и лишь потом считает `T` (`unified_optimizer/data_manager.py:84`).
Порознь rep1 = 98.909 %, rep2 = 98.461 %.
"""
from __future__ import annotations

import math

import numpy as np

from .compat import rfft_freqs, trapezoid

DC_NONE = "none"
DC_PRE_PULSE = "pre_pulse"
DC_FULL_MEAN = "full_mean"

REF_EARLIER = "earlier"
REF_LATER = "later"
REF_AVERAGE = "average"

CENTER_OWN = "own_peak"
CENTER_REF = "ref_peak"

# FWHM, при которой окно шире любой реальной трассы ⇒ фактическое выключение (ТЗ §3.3).
WINDOW_OFF_PS = 1000.0


class Settings(object):
    """Настройки расчёта. Значения по умолчанию — решения владельца 2026-08-12."""

    __slots__ = ("dc_mode", "dc_fraction", "window_on", "window_fwhm_ps",
                 "window_center", "ref_mode", "band_lo", "band_hi", "dyn_range_db")

    def __init__(self, dc_mode=DC_PRE_PULSE, dc_fraction=0.15,
                 window_on=False, window_fwhm_ps=20.0, window_center=CENTER_OWN,
                 ref_mode=REF_EARLIER, band_lo=0.2, band_hi=3.0, dyn_range_db=40.0):
        self.dc_mode = dc_mode
        self.dc_fraction = dc_fraction
        self.window_on = window_on
        self.window_fwhm_ps = window_fwhm_ps
        self.window_center = window_center
        self.ref_mode = ref_mode
        self.band_lo = band_lo
        self.band_hi = band_hi
        self.dyn_range_db = dyn_range_db

    def copy(self):
        s = Settings()
        for name in self.__slots__:
            setattr(s, name, getattr(self, name))
        return s

    def describe(self):
        """Строка для панели данных: активные режимы обязаны быть видны (ТЗ §3.3)."""
        dc = {DC_NONE: u"нет (как в ядре)",
              DC_PRE_PULSE: u"предымпульс %.0f %%" % (100 * self.dc_fraction),
              DC_FULL_MEAN: u"полное среднее"}[self.dc_mode]
        if not self.window_on or self.window_fwhm_ps >= WINDOW_OFF_PS:
            win = u"окно выключено"
        else:
            centre = (u"свой пик" if self.window_center == CENTER_OWN
                      else u"общий центр по пику референса")
            win = u"окно гаусс FWHM %.3g пс, центр — %s" % (self.window_fwhm_ps, centre)
        ref = {REF_EARLIER: u"ближайший ранний",
               REF_LATER: u"ближайший поздний",
               REF_AVERAGE: u"среднее двух"}[self.ref_mode]
        return u"%s; DC — %s; референс — %s" % (win, dc, ref)


# --------------------------------------------------------------------- базовые
def dc_level(E, mode, fraction=0.15):
    """Постоянная составляющая трассы.

    По умолчанию — среднее по **предымпульсному** участку (первые `fraction`
    отсчётов), а не по всей трассе: хвост осцилляций смещает полное среднее и
    вносит систематику в интеграл (ТЗ §3.2).
    """
    if mode == DC_NONE:
        return 0.0
    if mode == DC_FULL_MEAN:
        return float(np.mean(E))
    n = max(4, int(round(fraction * len(E))))
    return float(np.mean(E[:n]))


def peak_position(t, E, dc_fraction=0.15):
    """Положение пика |E(t)| с параболической интерполяцией по трём точкам.

    Без интерполяции разрешение по задержке ограничено шагом сетки (0.104 пс на
    наших прогонах), и сдвиг образца в доли шага не разглядеть (ТЗ §3.4).

    **Постоянная составляющая вычитается здесь всегда**, независимо от режима DC
    для `T`. Причина измерена, а не предположена: у `011_a90.txt` (зона гашения)
    смещение детектора −0.095 БОЛЬШЕ амплитуды импульса 0.077, поэтому максимум
    сырого |E| садится на случайную точку базовой линии и даёт 14.062 пс вместо
    14.539 пс — промах почти на пять шагов сетки. В ярком положении разницы нет.
    """
    E = np.asarray(E, dtype=np.float64)
    y = np.abs(E - dc_level(E, DC_PRE_PULSE, dc_fraction))
    i = int(np.argmax(y))
    if i == 0 or i == len(y) - 1:
        return float(t[i]), float(y[i])

    y0, y1, y2 = float(y[i - 1]), float(y[i]), float(y[i + 1])
    denom = y0 - 2.0 * y1 + y2
    if denom == 0.0:
        return float(t[i]), y1
    delta = 0.5 * (y0 - y2) / denom
    # Вершина параболы обязана лежать между соседями; иначе тройка не описывает пик.
    if not (-1.0 < delta < 1.0):
        return float(t[i]), y1
    dt = float(t[i + 1] - t[i - 1]) * 0.5
    return float(t[i]) + delta * dt, y1 - 0.25 * (y0 - y2) * delta


def gaussian_window(t, center_ps, fwhm_ps):
    """Гауссово окно, ширина задана как **FWHM** (так подписано и в интерфейсе).

    σ = FWHM / (2√(2 ln 2)). Назначение двойное: отсечь эхо от переотражений и
    подавить шумовую дорожку от линий поглощения воды (ТЗ §3.3).
    """
    sigma = float(fwhm_ps) / (2.0 * math.sqrt(2.0 * math.log(2.0)))
    if sigma <= 0:
        raise ValueError(u"ширина окна должна быть положительной")
    z = (np.asarray(t, dtype=np.float64) - float(center_ps)) / sigma
    return np.exp(-0.5 * z * z)


def energy(t, E, settings):
    """∫ (E − DC)² dt — числитель и знаменатель пропускания во времени."""
    d = dc_level(E, settings.dc_mode, settings.dc_fraction)
    return trapezoid((np.asarray(E, dtype=np.float64) - d) ** 2, t)


def dc_share(t, E, fraction=0.15):
    """Доля сырого ∫E², которую убирает вычитание постоянной составляющей.

        share = 1 − ∫(E−DC)²dt / ∫E²dt

    Считается именно так, а не как `DC²·окно / ∫E²`: последнее не ограничено
    единицей, потому что перекрёстный член `2·DC·∫(E−DC)dt` отрицателен, и на
    тёмной трассе «доля» выходит 101 %.

    Диагностика для режима `dc_mode='none'`. Измеренные значения на
    `test_grid_33_11`: яркое положение `003_a0` — **18 %**, зона гашения
    `011_a90` — **99.3 %**. В ярком положении смещение одинаково давит на
    образец и на референс, и отношение почти не страдает (98.91 % против
    99.42 %). В тёмном образец состоит из смещения на 99 %, а референс всё те
    же 18 % — отношение теряет смысл и `T` завышается **в 121 раз**
    (18.26 % вместо 0.151 %).

    Отсюда вывод за пределами этого инструмента: `T_time` из
    `a21_transmission_summary.py` считается без вычитания DC и потому пригоден
    **только для яркого положения**, ради которого он там и написан. Переносить
    его на тёмные углы нельзя.
    """
    E = np.asarray(E, dtype=np.float64)
    total = trapezoid(E ** 2, t)
    if total <= 0:
        return 0.0
    ac = trapezoid((E - dc_level(E, DC_PRE_PULSE, fraction)) ** 2, t)
    return float(1.0 - ac / total)


def _prepared(t, E, settings, centre_ps):
    """Трасса с вычтенным DC и, если включено, наложенным окном — для БПФ."""
    E = np.asarray(E, dtype=np.float64) - dc_level(E, settings.dc_mode, settings.dc_fraction)
    if settings.window_on and settings.window_fwhm_ps < WINDOW_OFF_PS:
        E = E * gaussian_window(t, centre_ps, settings.window_fwhm_ps)
    return E


def power_spectrum(t, E, settings, centre_ps):
    """(частоты в ТГц, |FFT|²) для подготовленной трассы."""
    Ew = _prepared(t, E, settings, centre_ps)
    dt = float(t[1] - t[0])
    spec = np.fft.rfft(Ew)
    return rfft_freqs(len(Ew), dt), np.abs(spec) ** 2


def noise_power(t, E, settings, centre_ps, fraction=0.15):
    """Оценка шумового пола спектра трассы, одно число на всю полосу.

    Берётся **предымпульсный** участок — там сигнала заведомо нет, только шум
    приёмного тракта. Для белого шума с дисперсией σ² математическое ожидание
    `|FFT|²` от N отсчётов равно σ²·N, что и возвращается.

    Зачем это нужно: выше примерно 1.6 ТГц отношение `T(ν)` на наших прогонах
    состоит из шума, делённого на сигнал, и скачет от нуля до единиц. Без
    отрисованного пола оператор не отличает «образец так пропускает» от «здесь
    уже нечего измерять», а по одному виду кривой это неразличимо.

    Оценка сознательно грубая (белый шум, один уровень на полосу): её задача —
    показать порядок величины на графике, а не служить метрологией. Окно, если
    включено, домножает шум так же, как сигнал, поэтому берётся от той же
    подготовленной трассы.
    """
    Ew = _prepared(t, E, settings, centre_ps)
    n = max(4, int(round(fraction * len(Ew))))
    sigma2 = float(np.var(Ew[:n]))
    return sigma2 * len(Ew)


# ------------------------------------------------------------------- результат
class PointResult(object):
    """Всё, что панель данных показывает по одной сигнальной трассе (ТЗ §4.3)."""

    __slots__ = ("trace", "refs", "T", "T_db", "energy_s", "energy_r",
                 "peak_s", "peak_r", "freqs", "T_nu", "T_noise", "warnings",
                 "settings", "ref_divergence")

    def __init__(self):
        self.trace = None
        self.refs = []          # список Trace, один или два
        self.T = float("nan")
        self.T_db = float("nan")
        self.energy_s = float("nan")
        self.energy_r = []      # по одному на референс
        self.peak_s = float("nan")
        self.peak_r = []
        self.freqs = None
        self.T_nu = None
        self.T_noise = None     # шумовой пол в тех же единицах, что T_nu
        self.warnings = []
        self.settings = None
        self.ref_divergence = None   # dict для режима усреднения, иначе None

    @property
    def ref_names(self):
        return u" / ".join([r.name for r in self.refs]) if self.refs else u"—"

    def summary_lines(self):
        """Текст панели данных — тот же и в GUI, и в CLI, чтобы числа не разъехались."""
        tr = self.trace
        out = [
            u"файл %s   угол %+d°%s   № %d" % (
                tr.name, tr.angle, (u"  rep%d" % tr.rep) if tr.rep > 1 else u"", tr.seq),
            u"референс %s%s" % (
                self.ref_names,
                u"  (среднее двух)" if len(self.refs) > 1 else u""),
        ]
        if self.peak_r:
            peaks = u" / ".join([u"%.3f" % p for p in self.peak_r])
            dpk = u" / ".join([u"%+.3f" % (self.peak_s - p) for p in self.peak_r])
            out.append(u"пик образца %.3f пс · референс %s пс · разность %s пс"
                       % (self.peak_s, peaks, dpk))
        ens = u" / ".join([u"%.4g" % e for e in self.energy_r])
        out.append(u"∫E² образец %.4g · референс %s" % (self.energy_s, ens))
        out.append(u"T = %.2f %%  (%+.3f дБ)" % (100.0 * self.T, self.T_db))
        if self.ref_divergence is not None:
            d = self.ref_divergence
            out.append(u"расхождение референсов: энергии %.3f, пики %+.3f пс"
                       % (d["energy_ratio"], d["peak_diff_ps"]))
        out.append(self.settings.describe())
        for w in self.warnings:
            out.append(u"⚠ " + w)
        return out


def resolve_refs(scan, trace, settings):
    """Какие референсы обслуживают данную трассу. -> (список Trace, предупреждения).

    Краевой случай ТЗ §3.5: у первой трассы прогона нет раннего референса, у
    последней — позднего. Откат на единственный доступный делается, но **всегда
    с явной пометкой** — молча подменять референс нельзя, иначе точка на границе
    прогона выглядит как обычная.
    """
    earlier, later = scan.neighbours_bg(trace.seq)
    warns = []

    if settings.ref_mode == REF_EARLIER:
        if earlier is not None:
            return [earlier], warns
        if later is not None:
            warns.append(u"раннего референса нет (начало прогона) — взят поздний %s"
                         % later.name)
            return [later], warns
    elif settings.ref_mode == REF_LATER:
        if later is not None:
            return [later], warns
        if earlier is not None:
            warns.append(u"позднего референса нет (конец прогона) — взят ранний %s"
                         % earlier.name)
            return [earlier], warns
    else:  # REF_AVERAGE
        both = [r for r in (earlier, later) if r is not None]
        if len(both) == 2:
            return both, warns
        if len(both) == 1:
            warns.append(u"второго референса нет (край прогона) — усреднять не с чем, "
                         u"взят единственный %s" % both[0].name)
            return both, warns

    warns.append(u"референсов в прогоне нет — пропускание не определено")
    return [], warns


def compute_point(scan, trace, settings):
    """Полный расчёт по одной сигнальной трассе.

    Усреднение двух референсов (решение владельца O-3, 2026-08-12):
      * **во времени** — усредняются два ИНТЕГРАЛА, а не трассы: расхождение
        задержек между референсами дало бы интерференционное занижение;
      * **в спектре** — усредняются два `|FFT|²`, то есть по мощности.
        Комплексное усреднение отвергнуто по измерению: на паре `002_bg`/`006_bg`
        оно завышает `T` в 8.6 раза на 2.76 ТГц, где отношение сигнал/шум уже
        мертво (деструктивная интерференция двух шумовых спектров в знаменателе).
        В рабочей полосе разница между способами 0.15 %.
    """
    res = PointResult()
    res.trace = trace
    res.settings = settings

    refs, warns = resolve_refs(scan, trace, settings)
    res.refs = refs
    res.warnings.extend(warns)
    if not refs:
        return res

    t_s, E_s = trace.t, trace.E
    res.peak_s = peak_position(t_s, E_s, settings.dc_fraction)[0]
    res.peak_r = [peak_position(r.t, r.E, settings.dc_fraction)[0] for r in refs]

    # ---- временная область
    res.energy_s = energy(t_s, E_s, settings)
    res.energy_r = [energy(r.t, r.E, settings) for r in refs]
    denom = sum(res.energy_r) / float(len(res.energy_r))
    if denom > 0:
        res.T = res.energy_s / denom
        res.T_db = 10.0 * math.log10(res.T) if res.T > 0 else float("-inf")
    else:
        res.warnings.append(u"энергия референса равна нулю — пропускание не определено")
        return res

    if len(refs) == 2:
        e1, e2 = res.energy_r
        res.ref_divergence = {
            "energy_ratio": (e1 / e2) if e2 > 0 else float("inf"),
            "peak_diff_ps": res.peak_r[0] - res.peak_r[1],
        }
        # Порог 5 % выбран по наблюдаемому разбросу соседних референсов в
        # test_grid_33_11 (доли процента): 5 % — это уже не шум, а событие.
        if abs(res.ref_divergence["energy_ratio"] - 1.0) > 0.05:
            res.warnings.append(
                u"референсы расходятся по энергии на %.1f %% — усреднение сомнительно"
                % (100.0 * abs(res.ref_divergence["energy_ratio"] - 1.0)))

    if settings.dc_mode == DC_NONE:
        share = dc_share(t_s, E_s, settings.dc_fraction)
        # Порог 50 % лежит в широком зазоре между измеренными режимами: яркое
        # положение даёт 18 %, зона гашения 99 %. Предупреждение обязано молчать
        # там, где смещение безобидно, иначе его перестанут читать.
        if share > 0.50:
            res.warnings.append(
                u"вычитание DC выключено, а постоянное смещение даёт %.1f %% интеграла: "
                u"это измерение смещения детектора, а не пропускания. "
                u"Включите DC — предымпульс." % (100.0 * share))

    if res.T > 1.0:
        # Прецедент: test_grid_40_20 дал T = 105.7 % из-за перепутанных sig/bg,
        # и это едва не ушло в отчёт. Молчать здесь нельзя (ТЗ §4.1).
        res.warnings.append(
            u"T = %.1f %% > 100 %%: проверьте, не перепутаны ли образец и референс, "
            u"и не сместился ли образец между трассой и фоном" % (100.0 * res.T))

    # ---- спектральная область
    _spectral(res, t_s, E_s, refs, settings)
    return res


def _spectral(res, t_s, E_s, refs, settings):
    """T(ν) = |FFT(E_s)|² / ⟨|FFT(E_r)|²⟩, обрезанный по полосе и динамике."""
    lengths = [len(E_s)] + [len(r.E) for r in refs]
    n = min(lengths)
    if len(set(lengths)) > 1:
        res.warnings.append(
            u"длины трасс различаются (%s точек) — для спектра взяты первые %d"
            % (u"/".join([str(x) for x in lengths]), n))

    dts = [float(t_s[1] - t_s[0])] + [float(r.t[1] - r.t[0]) for r in refs]
    if max(dts) - min(dts) > 1e-9:
        res.warnings.append(u"шаг по времени у образца и референса различается "
                            u"(%s пс) — спектры несопоставимы" % u"/".join(
                                [u"%.6f" % x for x in dts]))

    # Центр окна: свой пик у каждой трассы либо общий по пику первого референса.
    if settings.window_center == CENTER_REF:
        centre_s = res.peak_r[0]
        centres_r = [res.peak_r[0]] * len(refs)
    else:
        centre_s = res.peak_s
        centres_r = list(res.peak_r)

    freqs, P_s = power_spectrum(t_s[:n], E_s[:n], settings, centre_s)
    P_r = None
    for r, c in zip(refs, centres_r):
        _, p = power_spectrum(r.t[:n], r.E[:n], settings, c)
        P_r = p if P_r is None else P_r + p
    P_r = P_r / float(len(refs))            # усреднение ПО МОЩНОСТИ (O-3)

    band = (freqs >= settings.band_lo) & (freqs <= settings.band_hi)
    # Порог динамического диапазона задан в дБ по АМПЛИТУДЕ относительно максимума
    # спектра референса; в мощности это вдвое больше, отсюда множитель 2.
    floor = np.max(P_r) * (10.0 ** (-settings.dyn_range_db * 2.0 / 10.0))
    good = band & (P_r > floor) & np.isfinite(P_r)

    res.freqs = freqs[good]
    res.T_nu = P_s[good] / P_r[good]
    # Шумовой пол образца, приведённый к тем же единицам, что и T(ν): выше него
    # кривая — измерение, на нём — шум, делённый на сигнал референса.
    res.T_noise = noise_power(t_s[:n], E_s[:n], settings, centre_s,
                              settings.dc_fraction) / P_r[good]
    if res.freqs.size == 0:
        res.warnings.append(u"в полосе %.2f…%.2f ТГц нет точек выше порога динамического "
                            u"диапазона %.0f дБ" % (settings.band_lo, settings.band_hi,
                                                    settings.dyn_range_db))


def compute_track(scan, settings, progress=None):
    """Расчёт по всем непустым сигнальным трассам. -> список PointResult.

    Порядок — по сквозному номеру: он же порядок слайдера и порядок съёмки,
    то есть ось «измерение №» на нижней шкале интерфейса.
    """
    out = []
    signals = scan.signals
    for i, tr in enumerate(signals):
        out.append(compute_point(scan, tr, settings))
        if progress is not None:
            progress(i + 1, len(signals))
    return out
