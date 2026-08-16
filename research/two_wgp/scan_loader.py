#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Загрузчик прогонов в ЖУРНАЛЬНОМ формате прибора (`NNN_aXXX.txt`) — В-24.

Зачем
-----
Прибор пишет `NNN_a<угол>.txt` / `NNN_bg.txt` / `NNN_dark.txt`, а `DataManager`
разбирает только канонический шаблон `{dataset}_{angle}deg_rep{N}_{sig|bg}.txt`
(`unified_optimizer/data_manager.py:17`). Второй угловой прогон `specac`
(`data_pool/specac/exp_data_2/`, 68 трасс от 15.08) под шаблон не подходит и
молча пропускается — датасета для него не возникает.

Переименовывать нельзя: `data_pool/**` read-only для всех агентов. Прецедент
`normalize_test_grid_33_11.py` решал это КОПИЯМИ канонических имён рядом с
оригиналами, но копии — это запись в `data_pool`, требующая санкции владельца.
Здесь сопоставление делается **на стороне загрузчика**, без единого нового
файла на диске (путь, прямо рекомендованный в `HANDOFFS.md`, В-24).

Что делает
----------
Читает каталог, разбирает имена, сопоставляет каждой сигнальной трассе фон и
строит ту же пятёрку, что `exp_builder.load`:

    (angles, freqs, exp, valid, meta),  exp[угол, частота] = FFT(sig)/FFT(bg)

Правило фона — как в прецеденте и в маршрутном листе прибора: каждый `NNN_bg`
покрывает идущие ЗА ним сигнальные трассы до следующего фона. Политика
`nearest` (ближайший по журнальному номеру в обе стороны) доступна флагом:
для КОМПЛЕКСНОГО отношения ближе по времени — меньше фазовый дрейф, и разница
между политиками сама по себе есть мера этого дрейфа.

Повторы НЕ усредняются по умолчанию: у `specac/exp_data_2` на каждый угол их
2–4, и именно они дают прямую оценку шума УРОВНЯ УГЛА, которой в проекте не
было (A6/HN13 выносили вердикт по консервативной верхней границе, «на угол
один повтор»). Усреднение включается флагом — оно нужно лишь для сверки с
`DataManager`, который усредняет повторы (`data_manager.py:84-85`).

Зона A. `data_pool/**` и `unified_optimizer/**` — только чтение.

Запуск:
    .venv\\Scripts\\python.exe research/two_wgp/scan_loader.py --verify
    .venv\\Scripts\\python.exe research/two_wgp/scan_loader.py --dir data_pool/specac/exp_data_2
"""
from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[2]
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))

from unified_optimizer import config                                # noqa: E402
from unified_optimizer.optimizer_2d import get_transmission_spectra  # noqa: E402
from unified_optimizer.utils import find_auto_water_mask, normalize_angle  # noqa: E402

#: тот же контракт имён, что у инструмента прибора (track_viewer/core/scan.py:45)
TRACE_RE = re.compile(r"^(\d{3})_(?:a(-?\d+)|(bg|dark))\.txt$")

BG_PRECEDING = "preceding"
BG_NEAREST = "nearest"
#: практика владельца: против дрейфа усредняются ДВА ближайших фона, но только
#: по амплитуде спектра; фаза берётся у ближайшего. Комплексное усреднение
#: референсов в проекте отвергнуто по измерению (track_viewer/core/physics.py:506
#: — завышало T в 8.6 раза на 2.76 ТГц), а усреднение по мощности для
#: комплексного отношения не годится: фазу оно не даёт вовсе.
BG_TWO_ABS = "two_abs"


class Trace(object):
    """Одна трасса журнала."""

    __slots__ = ("seq", "kind", "angle", "path", "rep")

    def __init__(self, seq, kind, angle, path):
        self.seq = seq              # журнальный номер
        self.kind = kind            # 'sig' | 'bg' | 'dark'
        self.angle = angle          # float для 'sig', иначе None
        self.path = path
        self.rep = None             # номер повтора для данного угла

    def __repr__(self):
        return "<%03d %s %s>" % (self.seq, self.kind,
                                 "" if self.angle is None else "%g" % self.angle)


def scan_dir(directory):
    """Разбирает каталог -> список Trace, отсортированный по журнальному номеру.

    Файлы вне контракта имён (META, README, шаблоны) пропускаются молча — это
    та же политика, что у `DataManager`, и она делает недоснятый прогон
    обрабатываемым штатно.
    """
    directory = Path(directory)
    if not directory.is_dir():
        raise ValueError("нет такого каталога: %s" % directory)

    traces = []
    for p in sorted(directory.iterdir()):
        m = TRACE_RE.match(p.name)
        if not m:
            continue
        if p.stat().st_size == 0:               # недоснятая трасса
            continue
        seq = int(m.group(1))
        if m.group(2) is not None:
            traces.append(Trace(seq, "sig", normalize_angle(float(m.group(2))), p))
        else:
            traces.append(Trace(seq, m.group(3), None, p))
    traces.sort(key=lambda tr: tr.seq)

    seen = {}
    for tr in traces:                            # нумерация повторов по порядку
        if tr.kind == "sig":
            seen[tr.angle] = seen.get(tr.angle, 0) + 1
            tr.rep = seen[tr.angle]
    return traces


def assign_bg(traces, policy=BG_PRECEDING):
    """Сопоставляет каждой сигнальной трассе фоновую. -> {seq_sig: seq_bg}.

    Для `two_abs` первым идёт БЛИЖАЙШИЙ (его фаза и берётся), вторым — второй
    по близости; список из двух элементов.
    """
    bgs = [tr.seq for tr in traces if tr.kind == "bg"]
    if not bgs:
        raise ValueError("в каталоге нет ни одной фоновой трассы")
    out = {}
    for tr in traces:
        if tr.kind != "sig":
            continue
        if policy == BG_NEAREST:
            out[tr.seq] = min(bgs, key=lambda s: (abs(s - tr.seq), s))
        elif policy == BG_TWO_ABS:
            near = sorted(bgs, key=lambda s: (abs(s - tr.seq), s))
            out[tr.seq] = near[:2] if len(near) >= 2 else [near[0], near[0]]
        else:
            before = [s for s in bgs if s < tr.seq]
            out[tr.seq] = before[-1] if before else bgs[0]
    return out


def _transmission_two_abs(t_s, E_s, t_b1, E_b1, t_b2, E_b2):
    """t(nu) = FFT(sig) / [ср.модуль двух фонов * exp(i*arg ближайшего фона)].

    Приём владельца против дрейфа: усредняются АМПЛИТУДЫ двух ближайших фонов,
    фаза берётся у ближайшего. Амплитуда так усредняется честно (это оценка
    медленно плывущего уровня), а фаза — нет: её усреднение по двум трассам
    смешало бы разные задержки и испортило бы ровно ту величину, ради которой
    здесь всё и делается.
    """
    n = min(len(E_s), len(E_b1), len(E_b2))
    dt = float(t_s[1] - t_s[0])
    fs = np.fft.rfft(E_s[:n])
    f1 = np.fft.rfft(E_b1[:n])
    f2 = np.fft.rfft(E_b2[:n])
    freq = np.fft.rfftfreq(n, dt)
    amp = 0.5 * (np.abs(f1) + np.abs(f2))
    ref = amp * np.exp(1j * np.angle(f1))
    trans = np.zeros_like(fs, dtype=complex)
    good = amp > 1e-10
    trans[good] = fs[good] / ref[good]
    return freq, np.abs(fs), amp, trans


def _read(path):
    """(t, E) с вычтенной постоянной составляющей — как `data_manager.load_tds`.

    Третья колонка (позиция линии задержки) читается и игнорируется — решение
    владельца O-5, зафиксировано в `track_viewer/core/scan.py:99`.
    """
    raw = np.loadtxt(path, dtype=np.float64, ndmin=2)
    if raw.shape[1] < 2:
        raise ValueError("ожидались минимум 2 колонки: %s" % path)
    t, E = raw[:, 0], raw[:, 1].copy()
    E -= np.mean(E[:50])
    return t, E


def load_scan(directory, policy=BG_PRECEDING, average_reps=False,
              mask_angles=(0.0, 90.0), angles_limit=None):
    """Каталог журнального формата -> (angles, freqs, exp, valid, meta).

    Совместимо по смыслу с `exp_builder.load`: `exp` комплексная (углы x
    частоты), `valid` — водяная маска, растиражированная по углам.

    При `average_reps=False` строки соответствуют ОТДЕЛЬНЫМ трассам, и один и
    тот же угол может встречаться несколько раз — это и нужно для оценки шума
    уровня угла. `meta['rep']` хранит номер повтора каждой строки.
    """
    traces = scan_dir(directory)
    sigs = [tr for tr in traces if tr.kind == "sig"]
    if not sigs:
        raise ValueError("в каталоге нет сигнальных трасс: %s" % directory)
    bg_of = assign_bg(traces, policy)
    by_seq = {tr.seq: tr for tr in traces}

    if angles_limit is not None:
        lo, hi = angles_limit
        sigs = [tr for tr in sigs if lo <= tr.angle <= hi]
        if not sigs:
            raise ValueError("после фильтра angles_limit не осталось трасс")

    two_abs = (policy == BG_TWO_ABS)
    fields = []
    for tr in sigs:
        t_s, E_s = _read(tr.path)
        if two_abs:
            s1, s2 = bg_of[tr.seq]
            t_b, E_b = _read(by_seq[s1].path)
            t_b2, E_b2 = _read(by_seq[s2].path)
            fields.append((tr, t_s, E_s, t_b, E_b, t_b2, E_b2))
            continue
        t_b, E_b = _read(by_seq[bg_of[tr.seq]].path)
        if len(E_s) != len(E_b):
            raise ValueError("сигнал и фон разной длины: %s" % tr.path.name)
        fields.append((tr, t_s, E_s, t_b, E_b))

    if average_reps:
        # Усреднять ПОЛЯ во времени, а не готовые отношения: именно так делает
        # DataManager (data_manager.py:84-85), и это когерентное усреднение.
        # Среднее отношений не равно отношению средних — на этом сверка и
        # разошлась в первый раз (3.5e-2).
        grouped, out = {}, []
        for rec in fields:
            grouped.setdefault(rec[0].angle, []).append(rec)
        for a in sorted(grouped):
            g = grouped[a]
            n = min(len(r[2]) for r in g)
            E_s = np.mean([r[2][:n] for r in g], axis=0)
            E_b = np.mean([r[4][:n] for r in g], axis=0)
            rec = [g[0][0], g[0][1][:n], E_s, g[0][3][:n], E_b]
            if two_abs:
                rec += [g[0][5][:n], np.mean([r[6][:n] for r in g], axis=0)]
            out.append(tuple(rec))
        fields = out

    rows, freqs_common = [], None
    for rec in fields:
        tr, t_s, E_s, t_b, E_b = rec[:5]
        if two_abs:
            f, _, _, trans = _transmission_two_abs(t_s, E_s, t_b, E_b,
                                                   rec[5], rec[6])
        else:
            f, _, _, trans = get_transmission_spectra(t_s, E_s, t_b, E_b)
        if freqs_common is None:
            freqs_common = f
        elif len(f) != len(freqs_common):
            raise ValueError("трассы разной длины: %s (%d против %d точек); "
                             "окно сканирования менялось в течение прогона"
                             % (tr.path.name, len(f), len(freqs_common)))
        rows.append((tr, trans))

    f_max = getattr(config, "F_MAX_2D", config.F_MAX)
    idx = np.where((freqs_common >= config.F_MIN) & (freqs_common <= f_max))[0]
    analysis_freqs = freqs_common[idx]

    angles = np.array([tr.angle for tr, _ in rows], dtype=np.float64)
    reps = np.array([(1 if average_reps else tr.rep) for tr, _ in rows], dtype=int)
    exp = np.array([tr_[idx] for _, tr_ in rows], dtype=complex)

    order = np.argsort(angles, kind="stable")
    angles, exp, reps = angles[order], exp[order], reps[order]

    # --- водяная маска: по среднему ПО УГЛАМ, набор фиксируется явно
    if mask_angles is None:
        mask_rows = np.arange(len(angles))
    else:
        lo, hi = mask_angles
        mask_rows = np.array([i for i, a in enumerate(angles) if lo <= a <= hi])
        if mask_rows.size == 0:
            raise ValueError("mask_angles не покрывает ни одного угла")
    Ymean_db = np.mean(20 * np.log10(np.maximum(np.abs(exp[mask_rows, :]), 1e-12)),
                       axis=0)
    water_mask, intervals = find_auto_water_mask(analysis_freqs, Ymean_db)
    valid = np.tile(water_mask, (len(angles), 1))

    darks = [tr for tr in traces if tr.kind == "dark"]
    meta = {
        "directory": str(directory),
        "bg_policy": policy,
        "n_traces": int(len(rows)),
        "n_angles": int(len(angles)),
        "n_distinct_angles": int(len(np.unique(angles))),
        "angles": [float(a) for a in angles],
        "rep": [int(r) for r in reps],
        "n_freqs": int(len(analysis_freqs)),
        "f_min_thz": float(analysis_freqs[0]),
        "f_max_thz": float(analysis_freqs[-1]),
        "n_valid_freqs": int(np.sum(water_mask)),
        "water_mask": water_mask,
        "masked_intervals": [(float(a), float(b)) for a, b in intervals],
        "n_bg": int(sum(1 for tr in traces if tr.kind == "bg")),
        "n_dark": len(darks),
        "dark_paths": [str(tr.path) for tr in darks],
        "averaged_reps": bool(average_reps),
    }
    return angles, analysis_freqs, exp, valid, meta


def noise_floor_from_dark(directory, freqs_ref=None):
    """Шумовой пол по перекрытому пучку -> (freqs, средний |FFT|^2 трасс dark).

    Прямое измерение шума закрытого тракта, независимое от метода Нафтали по
    ВЧ-хвосту. Из трасс НЕ вычитается ничего: шум некогерентен, и вычитание
    добавило бы второй экземпляр вместо снятия первого.
    """
    darks = [tr for tr in scan_dir(directory) if tr.kind == "dark"]
    if not darks:
        return None, None
    spec, f = [], None
    for tr in darks:
        t, E = _read(tr.path)
        dt = float(t[1] - t[0])
        S = np.abs(np.fft.rfft(E)) ** 2
        f = np.fft.rfftfreq(len(E), dt) if f is None else f
        spec.append(S)
    return f, np.mean(spec, axis=0)


# =========================================================== сверка
def verify(verbose=True):
    """Сверка с `exp_builder.load` на датасете, доступном ОБОИМИ путями.

    `test_grid_33_11` лежит в журнальном формате И имеет канонические копии
    (сделаны A12), поэтому один и тот же прогон читается двумя независимыми
    путями. Расхождение обязано быть машинным нулём: физика одна.
    """
    import exp_builder                                             # noqa: E402

    ok = True
    ds, d = "test_grid_33_11", REPO / "data_pool" / "test_grid_33_11"

    a_ref, f_ref, e_ref, v_ref, m_ref = exp_builder.load(
        ds, angles_limit=None, mask_angles=(0.0, 90.0))
    a_my, f_my, e_my, v_my, m_my = load_scan(
        d, policy=BG_PRECEDING, average_reps=True, mask_angles=(0.0, 90.0))

    same_ang = len(a_ref) == len(a_my) and np.allclose(a_ref, a_my)
    if verbose:
        print("  углов: %d (эталон) против %d (загрузчик) — %s"
              % (len(a_ref), len(a_my), "совпало" if same_ang else "РАСХОЖДЕНИЕ"))
    ok &= same_ang

    same_f = len(f_ref) == len(f_my) and np.allclose(f_ref, f_my)
    if verbose:
        print("  частотная сетка: %d бинов — %s"
              % (len(f_my), "совпала" if same_f else "РАСХОЖДЕНИЕ"))
    ok &= same_f

    if same_ang and same_f:
        d_exp = float(np.max(np.abs(e_ref - e_my)))
        if verbose:
            print("  максимальное расхождение exp: %.3e" % d_exp)
        ok &= d_exp < 1e-12
        same_v = bool(np.array_equal(v_ref, v_my))
        if verbose:
            print("  водяная маска: %s (валидных %d против %d)"
                  % ("совпала" if same_v else "РАСХОЖДЕНИЕ",
                     m_ref["n_valid_freqs"], m_my["n_valid_freqs"]))
        ok &= same_v

    if verbose:
        print("  ИТОГ: %s" % ("сверка пройдена" if ok else "СВЕРКА ПРОВАЛЕНА"))
    return ok


def main():
    ap = argparse.ArgumentParser(description="Загрузчик журнального формата (В-24)")
    ap.add_argument("--verify", action="store_true",
                    help="сверка с exp_builder на test_grid_33_11")
    ap.add_argument("--dir", type=str, default=None, help="каталог прогона")
    ap.add_argument("--policy", choices=(BG_PRECEDING, BG_NEAREST, BG_TWO_ABS),
                    default=BG_PRECEDING)
    ap.add_argument("--average-reps", action="store_true")
    args = ap.parse_args()

    if args.verify:
        print("Сверка загрузчика с exp_builder (test_grid_33_11)")
        return 0 if verify() else 1

    if not args.dir:
        ap.print_help()
        return 0

    ang, freqs, exp, valid, meta = load_scan(
        args.dir, policy=args.policy, average_reps=args.average_reps)
    print("каталог: %s" % meta["directory"])
    print("  трасс %d, различных углов %d (%.0f..%.0f), фонов %d, dark %d"
          % (meta["n_traces"], meta["n_distinct_angles"], min(meta["angles"]),
             max(meta["angles"]), meta["n_bg"], meta["n_dark"]))
    print("  бинов %d (%.2f..%.2f ТГц), после водяной маски %d"
          % (meta["n_freqs"], meta["f_min_thz"], meta["f_max_thz"],
             meta["n_valid_freqs"]))
    uniq, cnt = np.unique(ang, return_counts=True)
    print("  повторов на угол: %s" % dict(zip(uniq.astype(int).tolist(),
                                              cnt.tolist())))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
