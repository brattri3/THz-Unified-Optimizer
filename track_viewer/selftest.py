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
import os
import shutil
import sys
import tempfile
import traceback

import numpy as np

from .core import physics as ph
from .core.compat import trapezoid
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
    # (069_a-40.txt). Прогон остановлен на нём, план предусматривал больше.
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


def main():
    if sys.stdout.encoding and sys.stdout.encoding.lower() not in ("utf-8", "utf8"):
        try:
            sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8")
        except Exception:                                       # noqa: BLE001
            pass

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
