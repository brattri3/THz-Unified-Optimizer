# -*- coding: utf-8 -*-
"""Консольный режим track_viewer: печать трека и экспорт CSV.

То же ядро, что и у графического интерфейса, но без окна. Нужен для трёх вещей:
сверки чисел с ядром без запуска GUI, работы по RDP на машине без графики и
выгрузки CSV для вёрстки в Origin Pro (ТЗ §4.5).

Запуск из корня репозитория:
    .venv\\Scripts\\python.exe -m track_viewer.cli data_pool/test_grid_33_11
    .venv\\Scripts\\python.exe -m track_viewer.cli data_pool/purewave --dc none
    .venv\\Scripts\\python.exe -m track_viewer.cli <каталог> --csv трек.csv
"""
from __future__ import annotations

import argparse
import io
import os
import sys

import numpy as np

from .core import physics as ph
from .core.scan import Scan

CSV_HEADER = (u"seq;file;angle_deg;rep;reference;energy_sig;energy_ref;"
              u"T;T_dB;T_freq;T_freq_dB;delta;delta_dB;"
              u"frac_in_band_sig;frac_in_band_ref;peak_sig_ps;peak_ref_ps;warnings")


def export_csv(points, path):
    """Выгрузка трека. Разделитель — `;`: Origin и Excel в русской локали ждут его.

    Колонки полосового пропускания присутствуют всегда; при выключенном
    Парсевале они пусты — а не заполнены нулями, которые Origin отложил бы
    на график как настоящие точки.
    """
    with io.open(path, "w", encoding="utf-8") as fh:
        fh.write(CSV_HEADER + u"\n")
        for p in points:
            tr = p.trace
            band = (u"%.10g;%.6f;%.10g;%.6f;%.8f;%.8f" % (
                p.T_freq, p.T_freq_db, p.delta, p.delta_db,
                p.frac_in_band_s, p.frac_in_band_r)
                if p.settings.parseval_on else u";;;;;")
            fh.write(u"%d;%s;%d;%d;%s;%.10g;%.10g;%.10g;%.6f;%s;%.6f;%s;%s\n" % (
                tr.seq, tr.name, tr.angle, tr.rep,
                u"+".join([r.name for r in p.refs]),
                p.energy_s,
                sum(p.energy_r) / float(len(p.energy_r)) if p.energy_r else float("nan"),
                p.T, p.T_db, band, p.peak_s,
                u"+".join([u"%.6f" % x for x in p.peak_r]),
                u" | ".join(p.warnings),
            ))
    return path


def build_settings(args):
    return ph.Settings(
        dc_mode={"pre": ph.DC_PRE_PULSE, "full": ph.DC_FULL_MEAN,
                 "none": ph.DC_NONE}[args.dc],
        dc_fraction=args.dc_fraction,
        window_on=args.window is not None,
        window_fwhm_ps=args.window if args.window is not None else 20.0,
        window_center=ph.CENTER_REF if args.common_centre else ph.CENTER_OWN,
        ref_mode={"earlier": ph.REF_EARLIER, "later": ph.REF_LATER,
                  "average": ph.REF_AVERAGE}[args.ref],
        band_lo=args.band[0], band_hi=args.band[1],
        parseval_on=bool(args.parseval or args.int_band is not None),
        int_lo=args.int_band[0] if args.int_band else 0.2,
        int_hi=args.int_band[1] if args.int_band else 3.0,
        window_in_time=args.window_in_time,
    )


def main(argv=None):
    ap = argparse.ArgumentParser(
        description=u"track_viewer: угловой трек THz-TDS в консоли")
    ap.add_argument("directory", help=u"каталог прогона")
    ap.add_argument("--dc", choices=["pre", "full", "none"], default="pre",
                    help=u"вычитание DC: pre — предымпульс (по умолчанию), "
                         u"full — полное среднее, none — как в ядре")
    ap.add_argument("--dc-fraction", type=float, default=0.15,
                    help=u"доля предымпульсного участка (по умолчанию 0.15)")
    ap.add_argument("--ref", choices=["earlier", "later", "average"],
                    default="earlier", help=u"выбор референса")
    ap.add_argument("--window", type=float, default=None, metavar=u"FWHM_ПС",
                    help=u"наложить гауссово окно указанной ширины (FWHM), пс")
    ap.add_argument("--common-centre", action="store_true",
                    help=u"общий центр окна по пику референса")
    ap.add_argument("--band", type=float, nargs=2, default=[0.2, 3.0],
                    metavar=("НИЗ", "ВЕРХ"), help=u"полоса спектра, ТГц")
    ap.add_argument("--parseval", action="store_true",
                    help=u"считать пропускание и по полосе частот, показать Δ")
    ap.add_argument("--int-band", type=float, nargs=2, default=None,
                    metavar=("НИЗ", "ВЕРХ"),
                    help=u"полоса ИНТЕГРИРОВАНИЯ, ТГц (включает --parseval); "
                         u"отдельна от --band, которая задаёт полосу показа")
    ap.add_argument("--window-in-time", action="store_true",
                    help=u"учитывать гауссово окно и во временно́м интеграле — "
                         u"без этого сравнение с полосой несимметрично")
    ap.add_argument("--csv", default=None, help=u"выгрузить трек в CSV")
    ap.add_argument("--seq", type=int, default=None,
                    help=u"подробности по одной трассе (сквозной номер)")
    args = ap.parse_args(argv)

    if sys.stdout.encoding and sys.stdout.encoding.lower() not in ("utf-8", "utf8"):
        try:
            sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8")
        except Exception:                                       # noqa: BLE001
            pass

    try:
        scan = Scan(args.directory)
    except ValueError as exc:
        print(u"ОШИБКА: %s" % exc)
        return 2

    c = scan.counts()
    print(u"каталог: %s" % os.path.abspath(args.directory))
    print(u"снято %d из %d сигнальных, %d из %d референсов, пустых %d, ошибок %d"
          % (c["sig_done"], c["sig_total"], c["bg_done"], c["bg_total"],
             c["empty"], c["errors"]))
    for tr in scan.errors:
        print(u"  ОШИБКА %s: %s" % (tr.name, tr.error))
    nxt = scan.next_unmeasured()
    if nxt is not None:
        print(u"следующий незаполненный файл: %s" % nxt.name)

    if not scan.signals:
        print(u"")
        print(u"Сигнальных трасс не найдено. Инструмент опознаёт только журнальную")
        print(u"нумерацию NNN_a<угол>.txt и NNN_bg.txt (см. docs/track_viewer/01_FORMATS.md).")
        return 1

    settings = build_settings(args)
    print(u"режим: %s" % settings.describe())
    print(u"")

    points = ph.compute_track(scan, settings)

    if args.seq is not None:
        sel = [p for p in points if p.trace.seq == args.seq]
        if not sel:
            print(u"трассы с номером %03d нет среди снятых" % args.seq)
            return 1
        for line in sel[0].summary_lines():
            print(line)
        return 0

    extra = u" %11s %9s" % (u"T полоса, %", u"Δ, дБ") if settings.parseval_on else u""
    print(u"%-4s %-16s %7s %4s  %-22s %9s %9s%s" %
          (u"№", u"файл", u"угол", u"rep", u"референс", u"T, %", u"T, дБ", extra))
    for p in points:
        cols = (u" %11.4f %+9.3f" % (100 * p.T_freq, p.delta_db)
                if settings.parseval_on else u"")
        print(u"%-4d %-16s %+7d %4d  %-22s %9.4f %9.3f%s%s" % (
            p.trace.seq, p.trace.name, p.trace.angle, p.trace.rep,
            p.ref_names, 100 * p.T, p.T_db, cols,
            u"  ⚠" if p.warnings else u""))

    warned = [p for p in points if p.warnings]
    if warned:
        print(u"")
        print(u"предупреждений: %d" % len(warned))
        for p in warned:
            for w in p.warnings:
                print(u"  %s: %s" % (p.trace.name, w))

    Ts = np.array([p.T for p in points])
    print(u"")
    print(u"T: минимум %.4f %% (%s), максимум %.4f %% (%s)"
          % (100 * Ts.min(), points[int(np.argmin(Ts))].trace.name,
             100 * Ts.max(), points[int(np.argmax(Ts))].trace.name))

    if args.csv:
        export_csv(points, args.csv)
        print(u"CSV: %s" % os.path.abspath(args.csv))
    return 0


if __name__ == "__main__":
    sys.exit(main())
