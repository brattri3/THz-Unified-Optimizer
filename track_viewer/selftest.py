# -*- coding: utf-8 -*-
"""Приёмка track_viewer — §7 ТЗ в виде автотестов.

Запуск из корня репозитория:
    .venv\\Scripts\\python.exe -m track_viewer.selftest

На целевой машине (Win7, комплект):
    python -m track_viewer.selftest

Своя мини-обвязка вместо pytest намеренно: на ПК спектрометра стоит голый
Python 3.8 с numpy и matplotlib, и приёмка обязана проходить там же, где
работает инструмент, а не только на машине разработчика.

Главный тест — `test_acceptance_number`: инструмент, расходящийся с ядром без
объяснения, вреднее отсутствия инструмента (ТЗ §7).
"""
from __future__ import annotations

import io
import math
import os
import shutil
import sys
import tempfile
import traceback

import numpy as np

from .core import fit_malus as fm
from .core import physics as ph
from .core.compat import Tee, setup_console, trapezoid
from .core.scan import Scan

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA = os.path.join(REPO, "data_pool")
GRID = os.path.join(DATA, "test_grid_33_11")
PUREWAVE = os.path.join(DATA, "purewave")

# Числа, посчитанные до написания кода и записанные в ACTIVITY.md 2026-08-12.
T_NO_DC = 0.989089        # 003_a0.txt / 002_bg.txt, без вычитания DC — число ТЗ §7.2
T_PRE_PULSE = 0.994231    # то же с DC по предымпульсу 15 %
T_REP2 = 0.984612         # 068_a0.txt / 066_bg.txt, без DC

_TESTS = []


def test(title):
    def deco(fn):
        _TESTS.append((title, fn))
        return fn
    return deco


def close(a, b, tol, what):
    if not abs(a - b) <= tol:
        raise AssertionError(u"%s: получено %.6f, ожидалось %.6f (допуск %g)"
                             % (what, a, b, tol))


# ---------------------------------------------------------------- §7.1 отбор
@test(u"§7.1 каталоги открываются, мета и нормализованные копии отсеяны")
def t_filtering(log):
    sc = Scan(GRID)
    c = sc.counts()
    log(u"test_grid_33_11: %d сигнальных, %d референсов, %d пустых, %d ошибок"
        % (c["sig_total"], c["bg_total"], c["empty"], c["errors"]))
    # 51 sig + 17 bg + 1 META = 69 = максимальный сквозной номер на диске
    # (069_a-40.txt). Это ПОЛНЫЙ план прогона, а не остановленный: те же
    # 51 + 17 стоят в MEASUREMENT_ORDER.md §3 и получаются из build_plan().
    if c["sig_total"] != 51 or c["bg_total"] != 17:
        raise AssertionError(u"ожидалось 51 сигнальная и 17 референсов")
    if c["empty"] or c["errors"]:
        raise AssertionError(u"в снятом прогоне не должно быть пустых и ошибочных")

    leaked = [t.name for t in sc.traces if "deg_rep" in t.name or "_sig" in t.name]
    if leaked:
        raise AssertionError(u"нормализованные копии просочились: %s" % leaked[:3])
    log(u"нормализованные копии (108 шт.), README.md, MEASUREMENT_ORDER.md, "
        u"000_META_TEMPLATE.txt, 001_META.txt, microscopy/ — отсеяны")

    sc2 = Scan(PUREWAVE)
    c2 = sc2.counts()
    log(u"purewave: %d сигнальных (снято %d), %d референсов (снято %d), %d пустых"
        % (c2["sig_total"], c2["sig_done"], c2["bg_total"], c2["bg_done"], c2["empty"]))
    if c2["sig_total"] != 51 or c2["bg_total"] != 18 or c2["empty"] != 29:
        raise AssertionError(u"purewave: ожидалось 51/18 и 29 пустых")
    if c2["errors"]:
        raise AssertionError(u"purewave: ошибок разбора быть не должно")
    log(u"недоснятая дорога обработана штатно: 29 пустышек пропущены молча")

    nxt = sc2.next_unmeasured()
    log(u"следующий незаполненный файл purewave: %s" % (nxt.name if nxt else u"—"))


@test(u"§7.1 старая схема NNN_META.txt и одиночный META.txt не считаются трассами")
def t_meta(log):
    d = tempfile.mkdtemp(prefix="tv_meta_")
    try:
        for name in ("001_META.txt", "META.txt", "000_META_TEMPLATE.txt",
                     "README.md", "MEASUREMENT_ORDER.md", "notes.txt"):
            with io.open(os.path.join(d, name), "w", encoding="utf-8") as fh:
                fh.write(u"seq: 001\nдата: 2026-08-12\n")
        shutil.copy(os.path.join(GRID, "003_a0.txt"), os.path.join(d, "003_a0.txt"))
        shutil.copy(os.path.join(GRID, "002_bg.txt"), os.path.join(d, "002_bg.txt"))
        sc = Scan(d)
        if len(sc.traces) != 2:
            raise AssertionError(u"ожидались ровно 2 трассы, найдено %d: %s"
                                 % (len(sc.traces), [t.name for t in sc.traces]))
        log(u"обе схемы META пропущены; архивный каталог открывается")
    finally:
        shutil.rmtree(d, ignore_errors=True)


# ------------------------------------------------------- §7.2 сверка с ядром
@test(u"§7.2 ГЛАВНЫЙ: T сходится с ядром и расхождение по DC объяснено")
def t_acceptance_number(log):
    sc = Scan(GRID)
    tr = [t for t in sc.signals if t.name == "003_a0.txt"][0]

    r_none = ph.compute_point(sc, tr, ph.Settings(dc_mode=ph.DC_NONE))
    r_pre = ph.compute_point(sc, tr, ph.Settings(dc_mode=ph.DC_PRE_PULSE))

    if r_none.ref_names != "002_bg.txt":
        raise AssertionError(u"референс должен быть 002_bg.txt, получен %s"
                             % r_none.ref_names)
    close(r_none.T, T_NO_DC, 1e-5, u"T без вычитания DC")
    close(r_pre.T, T_PRE_PULSE, 1e-5, u"T с DC по предымпульсу")

    log(u"без DC          T = %.4f %%  (%+.4f дБ)  ← число ТЗ §7.2"
        % (100 * r_none.T, r_none.T_db))
    log(u"предымпульс 15%% T = %.4f %%  (%+.4f дБ)  ← режим по умолчанию"
        % (100 * r_pre.T, r_pre.T_db))
    log(u"расхождение по трактовке DC: %+.3f п.п. (%+.4f дБ)"
        % (100 * (r_pre.T - r_none.T), r_pre.T_db - r_none.T_db))

    # Правило сопоставления: ядро использовало «ближайший ранний bg».
    # Проверено побайтно — test_grid_33_11_0deg_rep1_bg.txt == 002_bg.txt.
    tr2 = [t for t in sc.signals if t.name == "068_a0.txt"][0]
    r2 = ph.compute_point(sc, tr2, ph.Settings(dc_mode=ph.DC_NONE))
    close(r2.T, T_REP2, 1e-5, u"T для повтора 068_a0.txt")
    if r2.trace.rep != 2:
        raise AssertionError(u"068_a0.txt должен опознаваться как rep2, а не rep%d"
                             % r2.trace.rep)
    log(u"повтор того же угла: rep1 = %.4f %%, rep2 = %.4f %%, среднее по T = %.4f %%"
        % (100 * r_none.T, 100 * r2.T, 100 * (r_none.T + r2.T) / 2))
    log(u"расхождение с a21_transmission.json (99.08 %) объяснено: A21 идёт через")
    log(u"DataManager, а тот при двух повторах КОГЕРЕНТНО усредняет сами трассы")
    log(u"(unified_optimizer/data_manager.py:84) и лишь потом считает T. Здесь")
    log(u"усреднения нет вовсе — каждая трасса своя точка, как и нужно у прибора.")


@test(u"§7.2 своя трапеция совпадает с numpy на реальных данных")
def t_trapezoid(log):
    sc = Scan(GRID)
    tr = sc.signals[0]
    mine = trapezoid(tr.E ** 2, tr.t)
    ref = getattr(np, "trapezoid", None) or getattr(np, "trapz")
    theirs = float(ref(tr.E ** 2, tr.t))
    close(mine, theirs, abs(theirs) * 1e-12, u"трапеция")
    log(u"%.12e против %.12e (numpy %s) — совпало" % (mine, theirs, np.__version__))


@test(u"§3.2 без вычитания DC зона гашения завышается в 121 раз — инструмент предупреждает")
def t_dc_dominates_in_shadow(log):
    sc = Scan(GRID)
    tr = [t for t in sc.signals if t.name == "011_a90.txt"][0]
    r_none = ph.compute_point(sc, tr, ph.Settings(dc_mode=ph.DC_NONE))
    r_pre = ph.compute_point(sc, tr, ph.Settings(dc_mode=ph.DC_PRE_PULSE))
    share = ph.dc_share(tr.t, tr.E)

    log(u"011_a90.txt (гашение): постоянное смещение даёт %.1f %% сырого ∫E²"
        % (100 * share))
    log(u"T без DC = %.4f %%, T с предымпульсом = %.4f %% — завышение в %.0f раз"
        % (100 * r_none.T, 100 * r_pre.T, r_none.T / r_pre.T))
    if r_none.T / r_pre.T < 50:
        raise AssertionError(u"ожидалось завышение на два порядка")
    if not any(u"смещение" in w for w in r_none.warnings):
        raise AssertionError(u"инструмент обязан предупредить о доминировании DC")
    log(u"предупреждение выдано: %s"
        % [w for w in r_none.warnings if u"смещение" in w][0])

    # В ярком положении та же проверка молчит — предупреждение не должно быть шумом.
    bright = [t for t in sc.signals if t.name == "003_a0.txt"][0]
    rb = ph.compute_point(sc, bright, ph.Settings(dc_mode=ph.DC_NONE))
    if any(u"смещение" in w for w in rb.warnings):
        raise AssertionError(u"в ярком положении предупреждение излишне")
    log(u"в ярком положении доля смещения %.2f %% — предупреждения нет"
        % (100 * ph.dc_share(bright.t, bright.E)))
    log(u"следствие за пределами инструмента: T_time из a21_transmission_summary.py")
    log(u"считается без вычитания DC и пригоден ТОЛЬКО для яркого положения")


# ------------------------------------------------------------ §7.3 весь трек
@test(u"§7.3 слайдер проходит весь трек без исключений, спектры считаются")
def t_full_track(log):
    for path in (GRID, PUREWAVE):
        sc = Scan(path)
        pts = ph.compute_track(sc, ph.Settings())
        if not pts:
            raise AssertionError(u"%s: ни одной точки" % path)
        for p in pts:
            if not np.isfinite(p.T):
                raise AssertionError(u"%s: T не конечно" % p.trace.name)
            if p.freqs is None or p.freqs.size == 0:
                raise AssertionError(u"%s: пустой спектр" % p.trace.name)
            if not np.all(np.isfinite(p.T_nu)):
                raise AssertionError(u"%s: NaN в T(ν)" % p.trace.name)
        Ts = np.array([p.T for p in pts])
        log(u"%s: %d точек, T от %.4f %% до %.2f %%, спектр %d…%d отсчётов"
            % (os.path.basename(path), len(pts), 100 * Ts.min(), 100 * Ts.max(),
               min([p.freqs.size for p in pts]), max([p.freqs.size for p in pts])))


@test(u"§7.3 три режима референса и усреднение по двум")
def t_ref_modes(log):
    sc = Scan(GRID)
    tr = [t for t in sc.signals if t.name == "011_a90.txt"][0]
    out = {}
    for mode in (ph.REF_EARLIER, ph.REF_LATER, ph.REF_AVERAGE):
        r = ph.compute_point(sc, tr, ph.Settings(ref_mode=mode))
        out[mode] = r
        log(u"%-8s референс %-24s T = %.4f %%" % (mode, r.ref_names, 100 * r.T))
    avg = out[ph.REF_AVERAGE]
    if len(avg.refs) != 2:
        raise AssertionError(u"режим усреднения должен взять два референса")
    # Во времени усредняются ИНТЕГРАЛЫ, а не трассы (ТЗ §3.5).
    expect = avg.energy_s / ((avg.energy_r[0] + avg.energy_r[1]) / 2.0)
    close(avg.T, expect, 1e-12, u"усреднение двух интегралов")
    if avg.ref_divergence is None:
        raise AssertionError(u"в режиме усреднения нужна оценка расхождения референсов")
    log(u"расхождение референсов: энергии %.4f, пики %+.4f пс"
        % (avg.ref_divergence["energy_ratio"], avg.ref_divergence["peak_diff_ps"]))


@test(u"§3.5 край прогона: откат на единственный референс помечается явно")
def t_edge_reference(log):
    sc = Scan(GRID)
    last = sc.signals[-1]
    r = ph.compute_point(sc, last, ph.Settings(ref_mode=ph.REF_LATER))
    if not r.warnings:
        raise AssertionError(u"откат на ранний референс обязан быть помечен")
    log(u"%s, режим «поздний»: %s" % (last.name, r.warnings[0]))


# --------------------------------------------------------- §3.4 положение пика
@test(u"§3.4 параболическая интерполяция ловит сдвиг меньше шага сетки")
def t_parabola(log):
    dt = 0.104167
    t = np.arange(651) * dt
    shift = 0.03            # втрое меньше шага
    for true_peak in (14.0 + shift, 20.0 - shift):
        E = np.exp(-0.5 * ((t - true_peak) / 0.35) ** 2)
        got = ph.peak_position(t, E)[0]
        if abs(got - true_peak) > 0.005:
            raise AssertionError(u"пик %.4f вместо %.4f" % (got, true_peak))
    naive = t[int(np.argmax(np.exp(-0.5 * ((t - (14.0 + shift)) / 0.35) ** 2)))]
    log(u"истинный 14.030 пс: парабола %.4f, ближайший узел сетки %.4f"
        % (ph.peak_position(t, np.exp(-0.5 * ((t - 14.03) / 0.35) ** 2))[0], naive))


@test(u"§3.4 в зоне гашения пик ищется по трассе с вычтенным DC")
def t_peak_needs_dc(log):
    sc = Scan(GRID)
    tr = [t for t in sc.signals if t.name == "011_a90.txt"][0]
    dc = ph.dc_level(tr.E, ph.DC_PRE_PULSE, 0.15)
    raw = float(tr.t[int(np.argmax(np.abs(tr.E)))])
    ours = ph.peak_position(tr.t, tr.E)[0]
    log(u"011_a90.txt: DC = %+.4f, амплитуда импульса %.4f — смещение БОЛЬШЕ импульса"
        % (dc, np.abs(tr.E - dc).max()))
    log(u"максимум сырого |E| даёт %.3f пс, по трассе с вычтенным DC — %.3f пс"
        % (raw, ours))
    if abs(raw - ours) < 0.3:
        raise AssertionError(u"ожидалось расхождение — тест потерял смысл, проверить данные")
    if abs(ours - 14.583) > 0.05:
        raise AssertionError(u"пик %.3f, ожидалось около 14.583 пс" % ours)


# ------------------------------------------------------------- §3.3 окно
@test(u"§3.3 окно FWHM 1000 пс равносильно выключенному, узкое окно меняет спектр")
def t_window(log):
    sc = Scan(GRID)
    tr = [t for t in sc.signals if t.name == "003_a0.txt"][0]
    off = ph.compute_point(sc, tr, ph.Settings(window_on=False))
    wide = ph.compute_point(sc, tr, ph.Settings(window_on=True, window_fwhm_ps=1000.0))
    narrow = ph.compute_point(sc, tr, ph.Settings(window_on=True, window_fwhm_ps=5.0))

    n = min(off.T_nu.size, wide.T_nu.size)
    rel = float(np.max(np.abs(wide.T_nu[:n] - off.T_nu[:n]) / off.T_nu[:n]))
    if rel > 1e-9:
        raise AssertionError(u"FWHM 1000 пс должно выключать окно, расхождение %.2e" % rel)
    log(u"FWHM 1000 пс: спектр не изменился (макс. отклонение %.1e)" % rel)

    n2 = min(off.T_nu.size, narrow.T_nu.size)
    dev = np.abs(narrow.T_nu[:n2] - off.T_nu[:n2]) / off.T_nu[:n2]
    if float(np.max(dev)) < 1e-3:
        raise AssertionError(u"узкое окно обязано менять спектр")
    # Медиана, а не максимум: максимум набирается там, где T мало и относительная
    # разница взлетает до сотен процентов, ничего не говоря о типичном эффекте.
    log(u"FWHM 5 пс: спектр изменился на %.1f %% по медиане (макс. %.0f %%) — "
        u"эхо и водяные линии срезаны" % (100 * float(np.median(dev)),
                                          100 * float(np.max(dev))))
    log(u"режим окна виден в панели: «%s»" % narrow.settings.describe())


# ------------------------------------------- пропускание по полосе (Парсеваль)
@test(u"Парсеваль: энергия по полной полосе совпадает с суммой квадратов точно")
def t_parseval_exact(log):
    sc = Scan(GRID)
    tr = [t for t in sc.signals if t.name == "003_a0.txt"][0]
    s = ph.Settings(dc_mode=ph.DC_NONE)
    full = ph.band_energy(tr.t, tr.E, s, 0.0, 0.0, float("inf"))
    direct = float(np.sum(np.asarray(tr.E, dtype=float) ** 2))
    rel = abs(full - direct) / direct
    if rel > 1e-12:
        raise AssertionError(u"Σwₖ|FFT|²/N разошлась с ΣEₙ² на %.2e" % rel)
    log(u"Σ Eₙ² = %.12f, Σwₖ|FFTₖ|²/N = %.12f, расхождение %.1e"
        % (direct, full, rel))

    # Веса — не украшение: без них равенство ломается примерно вдвое.
    freqs, P = ph.power_spectrum(tr.t, tr.E, s, 0.0)
    naive = float(np.sum(P) / len(tr.E))
    log(u"без весов вышло бы %.6f вместо %.6f — ошибка в %.2f раза"
        % (naive, direct, direct / naive))


@test(u"Парсеваль на полной полосе схлопывает Δ; трапеция даёт остаток на концах")
def t_parseval_collapse(log):
    sc = Scan(GRID)
    tr = [t for t in sc.signals if t.name == "003_a0.txt"][0]
    wide = ph.Settings(parseval_on=True, int_lo=0.0, int_hi=1e9)
    r = ph.compute_point(sc, tr, wide)
    if abs(r.delta) > 4e-6:
        raise AssertionError(u"на полной полосе Δ должна быть ~10⁻⁶, вышло %.2e" % r.delta)
    log(u"полная полоса: T время %.10f, T полоса %.10f, Δ = %.2e"
        % (r.T, r.T_freq, r.delta))
    log(u"остаток — это разница трапеции и прямоугольной суммы, вклад дают концы трассы")

    # С окном концы трассы погашены, и способы суммирования сходятся.
    wide_w = ph.Settings(parseval_on=True, int_lo=0.0, int_hi=1e9,
                         window_on=True, window_fwhm_ps=20.0, window_in_time=True)
    rw = ph.compute_point(sc, tr, wide_w)
    if abs(rw.delta) >= abs(r.delta):
        raise AssertionError(u"окно обязано уменьшать остаток, стало %.2e" % rw.delta)
    log(u"с окном FWHM 20 пс остаток падает до %.2e — окно гасит концы" % rw.delta)


@test(u"Δ: в ярком положении ноль, в зоне гашения растёт при сужении полосы")
def t_delta_band(log):
    sc = Scan(GRID)
    bright = [t for t in sc.signals if t.name == "019_a10.txt"][0]
    dark = [t for t in sc.signals if t.name == "011_a90.txt"][0]

    r = ph.compute_point(sc, bright, ph.Settings(parseval_on=True,
                                                 int_lo=0.2, int_hi=3.0))
    if abs(r.delta_db) > 0.01:
        raise AssertionError(u"в ярком положении Δ обязана быть нулевой, вышло %+.3f дБ"
                             % r.delta_db)
    log(u"яркое 019_a10 (T = %.1f %%): Δ = %+.3f дБ — вся энергия в полосе"
        % (100 * r.T, r.delta_db))

    # Сужение полосы обязано УВЕЛИЧИВАТЬ расхождение по модулю: чем у́же полоса,
    # тем больше «пропускания» остаётся за её пределами. Проверяется монотонность,
    # а не совпадение с одним числом.
    seq = []
    for lo, hi in ((0.2, 3.0), (0.3, 1.2), (0.4, 1.0)):
        d = ph.compute_point(sc, dark, ph.Settings(parseval_on=True,
                                                   int_lo=lo, int_hi=hi))
        seq.append((lo, hi, d.delta_db, d.frac_in_band_s, d.frac_in_band_r))
    for (lo, hi, db, fs, fr) in seq:
        log(u"гашение 011_a90, полоса %.1f…%.1f: Δ = %+.3f дБ, "
            u"доля энергии в полосе образец %.1f %% референс %.1f %%"
            % (lo, hi, db, 100 * fs, 100 * fr))
    mags = [abs(x[2]) for x in seq]
    if not (mags[0] < mags[1] < mags[2]):
        raise AssertionError(u"сужение полосы обязано увеличивать |Δ|, вышло %s"
                             % [u"%.3f" % m for m in mags])
    if not (-0.9 < seq[0][2] < -0.3):
        raise AssertionError(u"на 0.2…3.0 ТГц ожидалось −0.3…−0.9 дБ, вышло %+.3f"
                             % seq[0][2])
    log(u"|Δ| растёт монотонно при сужении полосы — вклад набирается вне неё")


@test(u"учёт окна во временно́м интеграле меняет T и виден в описании режима")
def t_window_in_time(log):
    sc = Scan(GRID)
    dark = [t for t in sc.signals if t.name == "011_a90.txt"][0]
    base = ph.Settings(window_on=True, window_fwhm_ps=20.0, window_in_time=False)
    with_w = ph.Settings(window_on=True, window_fwhm_ps=20.0, window_in_time=True)
    a = ph.compute_point(sc, dark, base)
    b = ph.compute_point(sc, dark, with_w)
    if abs(a.T - b.T) < 1e-9:
        raise AssertionError(u"галочка обязана менять временно́е T")
    log(u"011_a90: без учёта окна T = %.5f %%, с учётом %.5f %% (%+.3f дБ)"
        % (100 * a.T, 100 * b.T, 10 * math.log10(b.T / a.T)))
    if u"окно учтено" not in b.settings.describe():
        raise AssertionError(u"режим обязан быть виден в панели данных")
    log(u"панель показывает: «%s»" % b.settings.describe())

    # Выключенное окно не должно давать эффекта, даже если галочка стоит.
    off = ph.compute_point(sc, dark, ph.Settings(window_on=False, window_in_time=True))
    plain = ph.compute_point(sc, dark, ph.Settings(window_on=False))
    if abs(off.T - plain.T) > 1e-15:
        raise AssertionError(u"без окна галочка не должна ни на что влиять")
    log(u"при выключенном окне галочка ничего не меняет — T тот же до 1e-15")


# ------------------------------------------- шумовой пол и перекрытый пучок
@test(u"три источника шумового пола; предымпульсная оценка измеренно завышена")
def t_noise_sources(log):
    sc = Scan(GRID)
    tr = [t for t in sc.signals if t.name == "003_a0.txt"][0]
    bg = [t for t in sc.backgrounds if t.name == "002_bg.txt"][0]
    s = ph.Settings()

    centre = ph.peak_position(bg.t, bg.E, s.dc_fraction)[0]
    pre = ph.noise_power(bg.t, bg.E, s, centre, s.dc_fraction)
    hf = ph.noise_power_hf(bg.t, bg.E, s, centre, 3.0)
    if not (pre > 0 and hf > 0):
        raise AssertionError(u"обе оценки обязаны быть положительными")
    if hf >= pre:
        raise AssertionError(u"ВЧ-хвост обязан быть НИЖЕ предымпульсной оценки")
    log(u"002_bg: предымпульс %.3g, ВЧ-хвост выше 3.0 ТГц %.3g — "
        u"разница %.1f дБ" % (pre, hf, 10 * math.log10(pre / hf)))
    log(u"предымпульсный участок содержит низкочастотный дрейф, а не только")
    log(u"белый шум, поэтому линия пола по нему нарисована заведомо выше реальной")

    # Каждый источник обязан дать конечный положительный пол и быть виден в панели.
    for src in (ph.NOISE_HF_TAIL, ph.NOISE_PRE_PULSE):
        r = ph.compute_point(sc, tr, ph.Settings(noise_source=src))
        if r.T_noise is None or not np.all(np.isfinite(r.T_noise)):
            raise AssertionError(u"источник %s дал неконечный пол" % src)
        if u"шумовой пол" not in r.settings.describe():
            raise AssertionError(u"источник обязан быть виден в панели данных")
    log(u"панель показывает источник: «%s»"
        % ph.Settings(noise_source=ph.NOISE_HF_TAIL).describe().split(u"; ")[-1])


@test(u"dark: парсер видит трассу перекрытого пучка и не пускает её в слайдер")
def t_dark_parse(log):
    from .core import road
    d = tempfile.mkdtemp(prefix="tv_dark_")
    try:
        for src, dst in (("002_bg.txt", "001_dark.txt"),
                         ("002_bg.txt", "002_bg.txt"),
                         ("003_a0.txt", "003_a0.txt"),
                         ("011_a90.txt", "004_a90.txt")):
            shutil.copy(os.path.join(GRID, src), os.path.join(d, dst))
        sc = Scan(d)
        if len(sc.darks) != 1 or sc.darks[0].name != "001_dark.txt":
            raise AssertionError(u"трасса перекрытого пучка не опознана")
        if any(t.name == "001_dark.txt" for t in sc.signals):
            raise AssertionError(u"dark не должен попадать в слайдер")
        if any(t.name == "001_dark.txt" for t in sc.backgrounds):
            raise AssertionError(u"dark не должен считаться референсом")
        c = sc.counts()
        log(u"каталог: %d сигнальных, %d фонов, %d перекрытых пучков"
            % (c["sig_done"], c["bg_done"], c["dark_done"]))
        log(u"ближайший dark к трассе №4: %s" % sc.nearest_dark(4).name)

        r = ph.compute_point(sc, sc.signals[0], ph.Settings(noise_source=ph.NOISE_DARK))
        if r.T_noise is None or not np.all(np.isfinite(r.T_noise)):
            raise AssertionError(u"пол по dark не посчитался")
        log(u"пол по dark посчитан, T = %.4f %% не изменилось — dark в трассы "
            u"не вмешивается" % (100 * r.T))
    finally:
        shutil.rmtree(d, ignore_errors=True)

    # Прежние каталоги обязаны дать прежние числа: расширение правила отбора
    # не должно ничего добавить и ничего не потерять.
    sc = Scan(GRID)
    c = sc.counts()
    if (c["sig_done"], c["bg_done"], c["dark_total"]) != (51, 17, 0):
        raise AssertionError(u"регрессия на test_grid_33_11: %s" % c)
    c2 = Scan(os.path.join(DATA, "purewave")).counts()
    if (c2["sig_done"], c2["dark_total"]) != (30, 0):
        raise AssertionError(u"регрессия на purewave: %s" % c2)
    log(u"регрессии нет: test_grid_33_11 те же 51/17, purewave те же 30, "
        u"файлов dark ни в одном каталоге нет")


@test(u"dark: выбор источника без dark-файлов откатывается на ВЧ-хвост с пометкой")
def t_dark_fallback(log):
    sc = Scan(GRID)
    tr = [t for t in sc.signals if t.name == "003_a0.txt"][0]
    r = ph.compute_point(sc, tr, ph.Settings(noise_source=ph.NOISE_DARK))
    notes = [w for w in r.warnings if u"перекрытого пучка" in w]
    if not notes:
        raise AssertionError(u"откат обязан быть помечен явно, а не молча")
    log(u"пометка: %s" % notes[0])

    ref = ph.compute_point(sc, tr, ph.Settings(noise_source=ph.NOISE_HF_TAIL))
    if float(np.max(np.abs(r.T_noise - ref.T_noise))) > 1e-12:
        raise AssertionError(u"откат обязан давать ровно тот же пол, что hf_tail")
    log(u"уровень совпал с ВЧ-хвостом до 1e-12 — откат, а не тихая замена")


@test(u"§5 генератор кладёт dark по настройке; гардрейл действует и на него")
def t_road_dark(log):
    from .core import road
    p = road.Plan(coarse_from=-30, coarse_to=30, coarse_step=15, fine_on=False,
                  mirror_fine=False, repeat_fine=False, anchors_rep2=(),
                  bg_mode=road.BG_START_END, dark_mode=road.DARK_START_END)
    recs = road.build_plan(p)
    names = [road.fname(r) for r in recs]
    if names[0] != "001_dark.txt" or names[-1] != "009_dark.txt":
        raise AssertionError(u"dark обязан открывать и замыкать прогон: %s" % names)
    # Порядок «сперва преграда, потом открытый пучок» — как это делается руками.
    if names[1] != "002_bg.txt":
        raise AssertionError(u"в начале dark идёт ПЕРЕД фоном: %s" % names[:3])
    if names[-2] != "008_bg.txt":
        raise AssertionError(u"в конце фон идёт ПЕРЕД dark: %s" % names[-3:])
    log(u"раскладка: %s … %s" % (u" ".join(names[:3]), u" ".join(names[-3:])))

    for mode, want in ((road.DARK_NONE, 0), (road.DARK_START, 1),
                       (road.DARK_START_END, 2)):
        p.dark_mode = mode
        got = road.summary(road.build_plan(p))["n_dark"]
        if got != want:
            raise AssertionError(u"режим %s дал %d dark вместо %d" % (mode, got, want))
    log(u"режимы «не снимать / один / начало-конец» дают 0, 1, 2 файла")

    p.dark_mode = road.DARK_EVERY_N
    p.dark_every = 2
    n_dark = road.summary(road.build_plan(p))["n_dark"]
    if n_dark < 3:
        raise AssertionError(u"каждые 2 трассы должно выйти больше двух dark")
    log(u"каждые 2 трассы: %d файлов перекрытого пучка" % n_dark)

    # Гардрейл §5.2 обязан действовать на новые имена без отдельного кода:
    # проверка идёт по списку плановых имён, а не по видам файлов.
    d = tempfile.mkdtemp(prefix="tv_dark_road_")
    try:
        p.dark_mode = road.DARK_START_END
        recs = road.build_plan(p)
        with io.open(os.path.join(d, "001_dark.txt"), "w") as fh:
            fh.write(u"0.0\t1.0\n0.1\t2.0\n")
        pv = road.preview(d, recs)
        if not pv.blocked:
            raise AssertionError(u"непустой dark обязан блокировать генерацию")
        log(u"непустой 001_dark.txt блокирует генерацию так же, как трасса: %d "
            u"конфликтов, не создано ничего" % len(pv.conflicts))
        if len(os.listdir(d)) != 1:
            raise AssertionError(u"каталог не должен измениться")

        # README прогона обязан объяснить, что dark не вычитается.
        text = road.readme_text(p, recs)
        if u"НЕ вычитается" not in text:
            raise AssertionError(u"README обязан предупредить о невычитании dark")
        log(u"README прогона содержит предупреждение о невычитании dark")
    finally:
        shutil.rmtree(d, ignore_errors=True)


# --------------------------------------------------- §4.1 предупреждение T > 1
@test(u"§4.1 перепутанные образец и референс дают T > 1 и явное предупреждение")
def t_over_unity(log):
    d = tempfile.mkdtemp(prefix="tv_swap_")
    try:
        # Тот самый сценарий, что дал test_grid_40_20 = 105.7 % и едва не ушёл в
        # отчёт: под именем образца лежит открытый пучок, под именем фона — образец.
        # Воспроизводим его в реалистичном масштабе (T чуть больше 100 %), а не
        # карикатурно: инструмент обязан ловить именно правдоподобный случай.
        shutil.copy(os.path.join(GRID, "002_bg.txt"), os.path.join(d, "003_a0.txt"))
        shutil.copy(os.path.join(GRID, "003_a0.txt"), os.path.join(d, "002_bg.txt"))
        sc = Scan(d)
        r = ph.compute_point(sc, sc.signals[0], ph.Settings())
        if r.T <= 1.0:
            raise AssertionError(u"ожидалось T > 1")
        if not any(u"> 100" in w for w in r.warnings):
            raise AssertionError(u"предупреждение о T > 100 %% не выдано")
        log(u"T = %.1f %%; предупреждение: %s" % (100 * r.T, r.warnings[0]))
    finally:
        shutil.rmtree(d, ignore_errors=True)


@test(u"§2.3 битый непустой файл — ошибка в панели, но трек не прерывается")
def t_broken_file(log):
    d = tempfile.mkdtemp(prefix="tv_broken_")
    try:
        shutil.copy(os.path.join(GRID, "002_bg.txt"), os.path.join(d, "002_bg.txt"))
        shutil.copy(os.path.join(GRID, "003_a0.txt"), os.path.join(d, "003_a0.txt"))
        with io.open(os.path.join(d, "004_a10.txt"), "w", encoding="utf-8") as fh:
            fh.write(u"это не таблица чисел\nсовсем не таблица\n")
        open(os.path.join(d, "005_a20.txt"), "wb").close()      # пустышка

        sc = Scan(d)
        if len(sc.errors) != 1 or sc.errors[0].name != "004_a10.txt":
            raise AssertionError(u"ожидалась ровно одна ошибка на 004_a10.txt")
        if len(sc.signals) != 1:
            raise AssertionError(u"на график должна попасть только одна трасса")
        pts = ph.compute_track(sc, ph.Settings())
        if len(pts) != 1:
            raise AssertionError(u"трек должен посчитаться, несмотря на битый файл")
        log(u"битый файл: %s" % sc.errors[0].error)
        log(u"пустышка 005_a20.txt пропущена молча; трек посчитан на 1 точке")
    finally:
        shutil.rmtree(d, ignore_errors=True)


@test(u"каталог без трасс даёт понятный отказ, а не пустое окно")
def t_empty_dir(log):
    d = tempfile.mkdtemp(prefix="tv_empty_")
    try:
        with io.open(os.path.join(d, "README.md"), "w", encoding="utf-8") as fh:
            fh.write(u"# пусто\n")
        sc = Scan(d)
        if sc.traces:
            raise AssertionError(u"трасс быть не должно")
        log(u"0 трасс — интерфейс обязан сказать об этом словами (см. gui.py)")
        log(u"так же выглядят каталоги старой схемы: test_grid_40_20 и specac")
        log(u"состоят только из нормализованных имён, журнальной нумерации в них нет")
    finally:
        shutil.rmtree(d, ignore_errors=True)


# -------------------------------------------------------- §5 генератор трека
@test(u"§5 генератор по умолчанию повторяет последовательность test_grid_33_11")
def t_road_reproduces(log):
    import re
    from .core import road

    # Историческая сверка идёт БЕЗ перекрытого пучка: в 2026 году его не снимали,
    # и файлов `NNN_dark.txt` в каталоге нет. Само значение по умолчанию с тех пор
    # изменилось (dark в начале и в конце, решение владельца 2026-08-12) — это
    # проверяется отдельно ниже, чтобы сдвиг нумерации не выдавался за регрессию.
    recs = road.build_plan(road.Plan(sample=u"test_grid_33_11",
                                     dark_mode=road.DARK_NONE))
    mine = [("bg", None) if r["kind"] == "bg" else ("sig", r["angle"]) for r in recs]

    names = sorted([n for n in os.listdir(GRID)
                    if re.match(r"^\d{3}_(a-?\d+|bg)\.txt$", n)])
    disk = []
    for n in names:
        m = re.match(r"^\d{3}_(a(-?\d+)|bg)\.txt$", n)
        disk.append(("bg", None) if m.group(2) is None
                    else ("sig", int(m.group(2))))

    if mine[:-1] != disk:
        for i, (a, b) in enumerate(zip(mine, disk)):
            if a != b:
                raise AssertionError(u"расхождение на позиции %d: %s против %s"
                                     % (i, a, b))
        raise AssertionError(u"длины планов не совпали: %d против %d"
                             % (len(mine) - 1, len(disk)))
    if mine[-1] != ("bg", None):
        raise AssertionError(u"последним обязан идти замыкающий фон")

    s = road.summary(recs)
    log(u"план: %d файлов = %d трасс + %d фонов" % (s["n_files"], s["n_sig"], s["n_bg"]))
    log(u"последовательность углов и фонов совпала с диском полностью")
    log(u"два намеренных отличия: нумерация сдвинута на 1 (META выведен из")
    log(u"сквозной нумерации, ТЗ §2.1) и добавлен замыкающий фон (ТЗ §5.1)")

    # Нынешнее умолчание: те же трассы плюс ровно два перекрытых пучка.
    now = road.build_plan(road.Plan(sample=u"test_grid_33_11"))
    s2 = road.summary(now)
    if s2["n_dark"] != 2:
        raise AssertionError(u"по умолчанию ожидались 2 перекрытых пучка, вышло %d"
                             % s2["n_dark"])
    if s2["n_sig"] != s["n_sig"] or s2["n_bg"] != s["n_bg"]:
        raise AssertionError(u"dark не должен менять состав трасс и фонов")
    if now[0]["kind"] != "dark" or now[-1]["kind"] != "dark":
        raise AssertionError(u"перекрытый пучок обязан открывать и замыкать прогон")
    log(u"нынешнее умолчание: те же %d трасс и %d фонов + 2 перекрытых пучка "
        u"(%s первым, %s последним)"
        % (s2["n_sig"], s2["n_bg"], road.fname(now[0]), road.fname(now[-1])))


@test(u"§5 параметры генератора работают: диапазоны, референсы, порядок, повторы")
def t_road_params(log):
    from .core import road

    p = road.Plan(coarse_from=-30, coarse_to=30, coarse_step=15,
                  fine_on=False, mirror_fine=False, repeat_fine=False,
                  anchors_rep2=(), bg_mode=road.BG_START,
                  dark_mode=road.DARK_NONE)
    recs = road.build_plan(p)
    angles = [r["angle"] for r in recs if r["kind"] == "sig"]
    if sorted(set(angles)) != [-30, -15, 0, 15, 30]:
        raise AssertionError(u"основной диапазон разложен неверно: %s" % angles)
    if len([r for r in recs if r["kind"] == "bg"]) != 1:
        raise AssertionError(u"в режиме «один в начале» должен быть ровно 1 фон")
    log(u"один референс в начале: %s" % [road.fname(r) for r in recs[:3]])

    p.bg_mode = road.BG_START_END
    recs = road.build_plan(p)
    if len([r for r in recs if r["kind"] == "bg"]) != 2:
        raise AssertionError(u"в режиме «начало и конец» должно быть 2 фона")
    if recs[0]["kind"] != "bg" or recs[-1]["kind"] != "bg":
        raise AssertionError(u"фоны должны стоять первым и последним")
    log(u"начало и конец: %s … %s" % (road.fname(recs[0]), road.fname(recs[-1])))

    p.order = road.ORDER_ASCENDING
    recs = road.build_plan(p)
    angles = [r["angle"] for r in recs if r["kind"] == "sig"]
    if angles != sorted(angles):
        raise AssertionError(u"обход по возрастанию не отсортирован: %s" % angles)
    log(u"обход по возрастанию: %s" % angles)

    p.order = road.ORDER_BLOCK
    p.repeats = 2
    recs = road.build_plan(p)
    n1 = len([r for r in road.build_plan(road.Plan(
        coarse_from=-30, coarse_to=30, coarse_step=15, fine_on=False,
        mirror_fine=False, repeat_fine=False, anchors_rep2=(),
        bg_mode=road.BG_START_END,
        dark_mode=road.DARK_NONE)) if r["kind"] == "sig"])
    n2 = len([r for r in recs if r["kind"] == "sig"])
    if n2 != 2 * n1:
        raise AssertionError(u"два повтора должны удвоить число трасс: %d против %d"
                             % (n2, 2 * n1))
    log(u"повторы: %d трасс при repeats=1, %d при repeats=2" % (n1, n2))

    bad = road.Plan(coarse_step=0).validate()
    if not bad:
        raise AssertionError(u"нулевой шаг обязан отвергаться")
    log(u"негодный план отвергается словами: %s" % bad[0])


@test(u"§5.2 ГАРДРЕЙЛ: непустой файл блокирует генерацию целиком")
def t_road_refuses(log):
    from .core import road

    d = tempfile.mkdtemp(prefix="tv_road_")
    try:
        plan = road.Plan(sample=u"проверка")
        recs = road.build_plan(plan)

        # Каталог с уже снятыми данными: кладём одну НЕПУСТУЮ трассу с плановым
        # именем — ровно тот случай, ради которого правило и написано.
        victim = road.fname([r for r in recs if r["kind"] == "sig"][0])
        shutil.copy(os.path.join(GRID, "003_a0.txt"), os.path.join(d, victim))
        before = sorted(os.listdir(d))

        pv = road.preview(d, recs)
        if not pv.blocked or victim not in pv.conflicts:
            raise AssertionError(u"конфликт не обнаружен")
        log(u"предпросмотр: %s" % pv.lines()[-2].strip())

        try:
            road.apply(d, recs, plan, pv)
        except ValueError as exc:
            log(u"создание отклонено: %s" % exc)
        else:
            raise AssertionError(u"apply обязан был отказать")

        after = sorted(os.listdir(d))
        if before != after:
            raise AssertionError(u"каталог изменился при отказе: %s" % after)
        if os.path.getsize(os.path.join(d, victim)) == 0:
            raise AssertionError(u"файл усечён — недопустимо")
        log(u"каталог не изменился: %d файл, %d байт — данные целы"
            % (len(after), os.path.getsize(os.path.join(d, victim))))
    finally:
        shutil.rmtree(d, ignore_errors=True)


@test(u"§5.2 создание кладёт только пустышки и не трогает то, что уже лежит")
def t_road_creates(log):
    from .core import road

    d = tempfile.mkdtemp(prefix="tv_road2_")
    try:
        plan = road.Plan(sample=u"проверка", coarse_from=-20, coarse_to=20,
                         coarse_step=10, fine_on=False, mirror_fine=False,
                         repeat_fine=False, anchors_rep2=())
        recs = road.build_plan(plan)

        created, skipped = road.apply(d, recs, plan)
        log(u"создано %d, пропущено %d" % (created, skipped))
        sizes = set()
        for rec in recs:
            path = os.path.join(d, road.fname(rec))
            if not os.path.exists(path):
                raise AssertionError(u"не создан %s" % road.fname(rec))
            sizes.add(os.path.getsize(path))
        if sizes != set([0]):
            raise AssertionError(u"созданы непустые файлы: размеры %s" % sizes)
        for extra in ("META.txt", "README.md"):
            if not os.path.exists(os.path.join(d, extra)):
                raise AssertionError(u"нет %s" % extra)
        log(u"все %d файлов нулевого размера; есть META.txt и README.md" % len(recs))

        # Свежая дорога должна открываться просмотрщиком как «ничего не снято».
        sc = Scan(d)
        c = sc.counts()
        if c["sig_done"] or c["bg_done"] or c["empty"] != len(recs):
            raise AssertionError(u"свежая дорога распознана неверно: %s" % c)
        log(u"просмотрщик видит её как %d пустышек, снято 0 — как и должно быть"
            % c["empty"])

        # Повторный запуск: ничего не создаёт, ничего не портит.
        with io.open(os.path.join(d, road.fname(recs[1])), "w", encoding="utf-8") as fh:
            fh.write(u"0.0\t1.0\n0.1\t2.0\n")          # «сняли» одну трассу
        pv = road.preview(d, recs)
        if not pv.blocked:
            raise AssertionError(u"снятая трасса обязана блокировать повтор")
        log(u"после «съёмки» одной трассы повторная генерация отказывает — "
            u"затереть измерение нельзя")
    finally:
        shutil.rmtree(d, ignore_errors=True)


# ------------------------------------------------- линеаризованный фит Малюса
def _synth_malus(theta_deg, t_perp, t_par, theta0_deg):
    """Когерентная модель: складываются ПОЛЯ, и лишь потом берётся квадрат."""
    t = (np.asarray(theta_deg, dtype=float) - theta0_deg) * (np.pi / 180.0)
    E = t_perp * np.cos(t) ** 2 + t_par * np.sin(t) ** 2
    return np.abs(E) ** 2


def _fit_points(**kw):
    """Трек test_grid_33_11 и его подгонка при заданных настройках."""
    sc = Scan(GRID)
    st = ph.Settings(**kw)
    pts = ph.compute_track(sc, st)
    return pts, fm.fit_track(pts, st)


@test(u"фит: order=4 точен для когерентной модели, order=2 — нет")
def t_fit_exact(log):
    th = np.arange(-90.0, 91.0, 10.0)
    t_perp, t_par, t0 = 0.93 + 0.11j, 0.031 * np.exp(1j * 0.7), 1.35
    U = _synth_malus(th, t_perp, t_par, t0)

    f4 = fm.fit_angular(th, U, fm.WEIGHT_UNIFORM, 4)
    if f4.info["resid_rms"] > 1e-14:
        raise AssertionError(u"order=4 обязан быть тождеством, rms = %.2e"
                             % f4.info["resid_rms"])
    pr = f4.params
    close(pr["theta0_deg_from_h2"], t0, 1e-9, u"theta0")
    close(pr["eta_amplitude"], abs(t_par) / abs(t_perp), 1e-12, u"утечка eta")
    close(pr["cos_dphi"], float(np.cos(np.angle(t_perp) - np.angle(t_par))),
          1e-10, u"cos относительной фазы")

    f2 = fm.fit_angular(th, U, fm.WEIGHT_UNIFORM, 2)
    if not f2.info["resid_rms"] > 1e3 * max(f4.info["resid_rms"], 1e-16):
        raise AssertionError(u"order=2 обязан оставить структурную невязку")
    log(u"rms: order=4 = %.2e (машинный нуль), order=2 = %.3e"
        % (f4.info["resid_rms"], f2.info["resid_rms"]))
    log(u"параметры восстановлены точно: eta = %.10f, cos(dphi) = %+.10f"
        % (pr["eta_amplitude"], pr["cos_dphi"]))


@test(u"фит: частотная зависимость не требует модели — сумма по частотам "
      u"остаётся многочленом 4-го порядка")
def t_fit_frequency_sum(log):
    th = np.arange(-90.0, 91.0, 10.0)
    Usum = np.zeros_like(th)
    for f in np.linspace(0.2, 1.5, 40):
        # eta(nu) и фаза(nu) взяты произвольными: суть в том, что их вид неважен
        tp = (0.95 - 0.1 * f) * np.exp(1j * 0.3 * f)
        tl = 0.02 * (1 + 0.5 * f) * np.exp(1j * (0.4 + 2 * np.pi * f * 0.19))
        Usum += _synth_malus(th, tp, tl, 1.35)
    fit = fm.fit_angular(th, Usum, fm.WEIGHT_UNIFORM, 4)
    rel = fit.info["resid_rms"] / Usum.max()
    if rel > 1e-14:
        raise AssertionError(u"сумма по частотам сломала базис: относит. rms = %.2e"
                             % rel)
    log(u"40 частот с произвольными eta(nu) и фазой(nu): относит. rms = %.2e" % rel)


@test(u"фит: широкополосный == спектрально-взвешенное среднее побинных")
def t_fit_identity(log):
    pts, _ = _fit_points()
    idx = [i for i, p in enumerate(pts) if p.freqs is not None and len(p.freqs)]
    grid = fm._common_grid(pts, idx)
    theta = np.array([pts[i].trace.angle for i in idx], dtype=float)
    Y = np.array([np.asarray(pts[i].T_nu, dtype=float)[
        np.searchsorted(pts[i].freqs, grid)] for i in idx])

    W = np.ones(len(grid)) / float(len(grid))       # любой вес; берём равномерный
    T_band = np.dot(Y, W)                            # та же величина во «времени»
    sigma = np.maximum(T_band, T_band.max() * 1e-6)  # ОДИН вектор весов на все бины

    wide = fm.fit_angular(theta, T_band, fm.WEIGHT_RELATIVE, 4)
    fine = fm.fit_spectral(theta, Y, grid, common_sigma=sigma)
    avg = np.dot(fine.coef, W)
    rel = np.max(np.abs(wide.p - avg) / np.maximum(np.abs(wide.p), 1e-30))
    if rel > 1e-9:
        raise AssertionError(u"тождество нарушено: относит. расхождение %.2e" % rel)
    log(u"5 коэффициентов совпали до %.1e — модель частотной зависимости не нужна"
        % rel)
    log(u"бинов %d, углов %d" % (len(grid), len(theta)))


@test(u"фит: два независимых угла столика согласованы, A4/A2 близко к 0.25")
def t_fit_two_offsets(log):
    _, fit = _fit_points(fit_on=True)
    pr = fit.time.params
    d = abs(pr["theta0_h4_minus_h2_deg"])
    if d > 0.05:
        raise AssertionError(u"офсеты из 2-й и 4-й гармоник разошлись на %.4f град" % d)
    ratio = pr["harmonic4_over_harmonic2"]
    if not (0.245 <= ratio <= 0.255):
        raise AssertionError(u"A4/A2 = %.5f вне 0.245…0.255" % ratio)
    log(u"theta0: %+.4f (2-я гармоника) и %+.4f (4-я), расхождение %+.4f град"
        % (pr["theta0_deg_from_h2"], pr["theta0_deg_from_h4_mod45"],
           pr["theta0_h4_minus_h2_deg"]))
    log(u"A4/A2 = %.5f при идеальных 0.25 для cos^4 — 4-я гармоника реальна"
        % ratio)


@test(u"фит: числа сходятся с сессией A (A13, модельно-независимый разбор)")
def t_fit_vs_session_a(log):
    _, fit = _fit_points(fit_on=True)
    pr = fit.time.params
    # research/results/two_wgp/a13_test_grid_33_11.json, блок harmonics:
    # theta0 = -1.2770 +- 0.0354 град, A4/A2 = 0.25183. Разбор гармонический,
    # то есть не затронут отзывом вывода A13 по закону H5 (см. A14).
    A13_THETA0, A13_RATIO = -1.2770, 0.25183
    for got, want, tol, what in (
            (pr["theta0_deg_from_h2"], A13_THETA0, 0.05, u"theta0 из 2-й гармоники"),
            (pr["theta0_deg_from_h4_mod45"], A13_THETA0, 0.05, u"theta0 из 4-й"),
            (pr["harmonic4_over_harmonic2"], A13_RATIO, 0.01, u"A4/A2")):
        close(got, want, tol, what)
    log(u"theta0 = %+.4f против %+.4f у A13; A4/A2 = %.5f против %.5f"
        % (pr["theta0_deg_from_h2"], A13_THETA0,
           pr["harmonic4_over_harmonic2"], A13_RATIO))
    log(u"расчёт независимый: у нас своя предобработка DC и свой выбор референса")


@test(u"фит: выбор весов меняет результат — равные веса слепы к зоне гашения")
def t_fit_weights(log):
    out = {}
    for mode in (fm.WEIGHT_UNIFORM, fm.WEIGHT_RELATIVE, fm.WEIGHT_NOISE):
        _, fit = _fit_points(fit_on=True, fit_weights=mode)
        out[mode] = fit.time.params
    rel = out[fm.WEIGHT_RELATIVE]["eta_amplitude"]
    uni = out[fm.WEIGHT_UNIFORM]["eta_amplitude"]
    if not abs(uni - rel) / rel > 0.01:
        raise AssertionError(u"веса обязаны влиять на утечку, вышло %.6f против %.6f"
                             % (uni, rel))
    for mode in out:
        if not np.isfinite(out[mode]["eta_amplitude"]):
            raise AssertionError(u"режим весов %s не дал утечки" % mode)
    log(u"eta: равные %.5f, относительные %.5f, по шумовому полу %.5f"
        % (uni, rel, out[fm.WEIGHT_NOISE]["eta_amplitude"]))
    log(u"cos(dphi): равные %+.4f против относительных %+.4f — разница не косметическая"
        % (out[fm.WEIGHT_UNIFORM]["cos_dphi"], out[fm.WEIGHT_RELATIVE]["cos_dphi"]))
    log(u"у сессии A равные веса на series1 дали |t_par|^2 < 0 — физически невозможное")


@test(u"фит: на неполном треке отказывается считать, а не выдаёт мусор")
def t_fit_guardrail(log):
    pts, _ = _fit_points()
    st = ph.Settings(fit_on=True)
    short = pts[:4]
    fit = fm.fit_track(short, st)
    if fit.ok:
        raise AssertionError(u"фит по %d точкам обязан быть отклонён" % len(short))
    if not fit.notes:
        raise AssertionError(u"отказ обязан быть объяснён, а не молчалив")
    log(u"%d точек: %s" % (len(short), fit.notes[0]))
    full = fm.fit_track(pts, st)
    if not full.ok:
        raise AssertionError(u"полный трек обязан фитироваться")
    log(u"%d точек: фит выполнен, cond = %.2f"
        % (len(pts), full.time.diag["cond_XtX"]))


@test(u"фит: обрезка углового диапазона портит обусловленность плана")
def t_fit_conditioning(log):
    sc = Scan(GRID)
    ang = np.array([tr.angle for tr in sc.signals], dtype=float)
    full = fm.design_diagnostics(ang, 4)
    half = fm.design_diagnostics(ang[(ang >= 0) & (ang <= 90)], 4)
    if not (full["cond_XtX"] < half["cond_XtX"]):
        raise AssertionError(u"обрезка обязана ухудшать обусловленность: %.2f и %.2f"
                             % (full["cond_XtX"], half["cond_XtX"]))
    gain = half["sigma_theta0_deg"] / full["sigma_theta0_deg"]
    log(u"cond: %.2f на %d углах против %.2f на %d — хуже в %.1f раза"
        % (full["cond_XtX"], full["n_angles"], half["cond_XtX"], half["n_angles"],
           half["cond_XtX"] / full["cond_XtX"]))
    log(u"предсказанная sigma(theta0) хуже в %.2f раза только из-за плана углов"
        % gain)


@test(u"фит: гармоники 6-8 считаются в метрике весов и ловят выход за модель")
def t_fit_high_harmonics(log):
    th = np.arange(-90.0, 91.0, 5.0)
    U = _synth_malus(th, 0.95 + 0.0j, 0.03 * np.exp(1j * 0.7), 1.35)

    # Без чужих гармоник доля обязана держаться на случайном уровне.
    f0 = fm.fit_angular(th, U, fm.WEIGHT_RELATIVE, 4)
    if not (f0.hh_share <= f0.hh_chance * 1.5):
        raise AssertionError(u"на чистых данных доля %.4f превысила случайный "
                             u"уровень %.4f" % (f0.hh_share, f0.hh_chance))

    # С подмешанной 6-й гармоникой — обязана вырасти на порядок.
    U6 = U + 0.01 * U.max() * np.cos(6 * np.radians(th - 1.35))
    f6 = fm.fit_angular(th, U6, fm.WEIGHT_RELATIVE, 4)
    if not (f6.hh_share > 10 * max(f0.hh_share, 1e-6)):
        raise AssertionError(u"подмешанная 6-я гармоника не замечена: %.4f против %.4f"
                             % (f6.hh_share, f0.hh_share))
    log(u"чистые данные %.2f %% при случайных %.2f %%; с 6-й гармоникой %.2f %%"
        % (100 * f0.hh_share, 100 * f0.hh_chance, 100 * f6.hh_share))

    # Метрика обязана быть той же, в которой минимизировался фит. Невзвешенный
    # счёт на реальных данных завышает долю вчетверо — так и была допущена
    # ошибка, из-за которой в отчёт едва не ушёл несуществующий выход за модель.
    _, fit = _fit_points(fit_on=True)
    f = fit.time
    naive = fm.high_harmonic_share(f.theta, f.info["resid"])[0]
    if not (naive > 2 * f.hh_share):
        raise AssertionError(u"тест метрики выродился: %.4f против %.4f"
                             % (naive, f.hh_share))
    log(u"на реальных данных: в метрике весов %.1f %% при случайных %.1f %%, "
        u"невзвешенно было бы %.1f %%"
        % (100 * f.hh_share, 100 * f.hh_chance, 100 * naive))
    log(u"выхода за модель на test_grid_33_11 НЕТ — доля ниже случайного уровня")


@test(u"фит: побинный разбор даёт спектры каналов и считает нарушения "
      u"неравенства Коши-Буняковского")
def t_fit_spectral(log):
    _, fit = _fit_points(fit_on=True, fit_spectral_on=True, band_hi=2.0)
    sp = fit.spectral
    if sp is None:
        raise AssertionError(u"побинный фит не выполнен")
    ch = sp.ch
    if sp.n_bins < 50:
        raise AssertionError(u"бинов слишком мало: %d" % sp.n_bins)
    if sp.n_neg > 0.1 * sp.n_bins:
        raise AssertionError(u"|t_par|^2 < 0 в %d бинах из %d — фит не физичен"
                             % (sp.n_neg, sp.n_bins))
    # восстановленная кривая обязана лечь на измеренную в ярком положении
    Tfit = sp.transmission(0.0)
    if not np.all(np.isfinite(Tfit)):
        raise AssertionError(u"восстановленная T(nu) содержит не-числа")
    lo = sp.freqs < 1.0
    log(u"бинов %d, |t_par|^2 < 0 в %d, нарушений |cos| <= 1 в %d"
        % (sp.n_bins, sp.n_neg, sp.n_cs_violations))
    log(u"утечка растёт с частотой: eta = %.4f ниже 1 ТГц и %.4f выше"
        % (float(np.nanmedian(ch["eta_amplitude"][lo])),
           float(np.nanmedian(ch["eta_amplitude"][~lo]))))
    log(u"восстановленная T(nu) при 0 град: от %.4f до %.4f"
        % (Tfit.min(), Tfit.max()))


# ------------------------------------------------------------- инварианты кода
@test(u"инвариант: в core/ нет ни tkinter, ни matplotlib")
def t_core_headless(log):
    core = os.path.join(os.path.dirname(os.path.abspath(__file__)), "core")
    bad = []
    for name in sorted(os.listdir(core)):
        if not name.endswith(".py"):
            continue
        with io.open(os.path.join(core, name), "r", encoding="utf-8") as fh:
            for i, line in enumerate(fh, 1):
                s = line.strip()
                if s.startswith(("import ", "from ")) and (
                        "tkinter" in s or "matplotlib" in s):
                    bad.append(u"%s:%d %s" % (name, i, s))
    if bad:
        raise AssertionError(u"ядро потеряло независимость от GUI: %s" % bad)
    log(u"проверено файлов: %d — приёмка идёт без графического окружения"
        % len([n for n in os.listdir(core) if n.endswith(".py")]))


def main(argv=None):
    argv = sys.argv[1:] if argv is None else argv
    log_path = None
    if "--log" in argv:
        i = argv.index("--log")
        if i + 1 < len(argv):
            log_path = argv[i + 1]

    setup_console()
    log_fh = None
    if log_path:
        # Лог — UTF-8 С BOM (`utf-8-sig`): его открывают Блокнотом, а Блокнот
        # Windows 7 распознаёт UTF-8 без BOM ненадёжно и показал бы ровно те
        # кракозябры, ради избавления от которых лог и заведён.
        # Если файл не создался — это не повод не прогонять приёмку.
        try:
            log_fh = io.open(log_path, "w", encoding="utf-8-sig")
            sys.stdout = Tee(sys.stdout, log_fh)
        except (IOError, OSError) as exc:
            print(u"не удалось открыть лог %s: %s" % (log_path, exc))

    console = sys.stdout
    try:
        code = _run()
    finally:
        if log_fh is not None:
            # Сначала вернуть консоль, потом закрыть файл: иначе последние
            # строки ушли бы в закрытый поток.
            sys.stdout = console._s if isinstance(console, Tee) else console
            log_fh.close()
            print(u"")
            print(u"полный вывод записан в %s" % os.path.abspath(log_path))
            print(u"этот файл открывается Блокнотом без искажений — его и присылать")
    return code


def _run():
    print(u"track_viewer — приёмка по §7 ТЗ")
    print(u"данные: %s" % DATA)
    print(u"python %s, numpy %s" % (sys.version.split()[0], np.__version__))
    print(u"")

    passed = failed = 0
    for title, fn in _TESTS:
        lines = []
        try:
            fn(lines.append)
        except Exception as exc:                                # noqa: BLE001
            failed += 1
            print(u"[ПРОВАЛ] %s" % title)
            for ln in lines:
                print(u"          %s" % ln)
            print(u"          %s: %s" % (type(exc).__name__, exc))
            if os.environ.get("TV_TRACEBACK"):
                traceback.print_exc()
        else:
            passed += 1
            print(u"[ ok  ] %s" % title)
            for ln in lines:
                print(u"          %s" % ln)
        print(u"")

    total = passed + failed
    print(u"=== пройдено %d из %d ===" % (passed, total))
    return 0 if failed == 0 else 1


if __name__ == "__main__":
    sys.exit(main())
