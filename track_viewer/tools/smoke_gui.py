# -*- coding: utf-8 -*-
"""Дымовой прогон интерфейса без участия человека.

Строит окно, прячет его, открывает каталог, дёргает слайдер и все переключатели
и сохраняет снимок каждой конфигурации в PNG. Ловит то, чего не видит
`selftest.py`: опечатки в именах виджетов, развалившуюся раскладку, падения при
переключении режимов, пустой каталог. Полноценную ручную проверку не заменяет —
цвета, читаемость и удобство слайдера смотрит человек.

Запуск из корня репозитория:
    .venv\\Scripts\\python.exe track_viewer/tools/smoke_gui.py
    .venv\\Scripts\\python.exe track_viewer/tools/smoke_gui.py --out <каталог>
"""
from __future__ import annotations

import argparse
import io
import os
import sys
import tempfile
import tkinter as tk

# Консоль Windows отдаёт cp1251, а в отчёте есть «Δ» и «⚠» — без перевода потока
# в UTF-8 прогон падает на печати, хотя сам интерфейс исправен.
if sys.stdout.encoding and sys.stdout.encoding.lower() not in ("utf-8", "utf8"):
    try:
        sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8")
    except Exception:                                            # noqa: BLE001
        pass

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
REPO = os.path.dirname(HERE)
if REPO not in sys.path:
    sys.path.insert(0, REPO)

from track_viewer.gui import TrackViewer                     # noqa: E402

GRID = os.path.join(REPO, "data_pool", "test_grid_33_11")
PUREWAVE = os.path.join(REPO, "data_pool", "purewave")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default=os.path.join(tempfile.gettempdir(), "tv_smoke"))
    args = ap.parse_args()
    if not os.path.isdir(args.out):
        os.makedirs(args.out)

    root = tk.Tk()
    root.geometry("1320x900")
    root.withdraw()                     # окно не мешает, отрисовка идёт в фигуру
    app = TrackViewer(root, GRID)
    root.update()

    shots = []

    def shot(name):
        path = os.path.join(args.out, name + ".png")
        app.fig.savefig(path, dpi=90)
        shots.append(path)
        print(u"  %-34s %s" % (name, path))

    print(u"каталог %s: %d точек" % (GRID, len(app.points)))
    shot("01_track_default")

    app.select(len(app.points) // 2)
    root.update()
    shot("02_slider_middle")

    print(u"проход слайдером по всем точкам…")
    for i in range(len(app.points)):
        app.select(i)
    root.update()
    print(u"  без исключений, %d шагов" % len(app.points))
    shot("03_slider_last")

    app.var_dc.set(u"нет (как в ядре)")
    app.recompute()
    root.update()
    shot("04_dc_none")

    app.var_dc.set(u"предымпульс")
    app.var_ref.set(u"среднее двух")
    app.var_window.set(True)
    app.var_fwhm.set("12")
    app.var_centre.set(u"общий по пику референса")
    app.recompute()
    root.update()
    shot("05_window_avg_ref")

    app.var_window.set(False)
    app.var_ref.set(u"ближайший ранний")
    app.var_band_hi.set("1.5")
    app.recompute()
    root.update()
    shot("06_band_1p5")

    # Некорректный ввод не должен ронять пересчёт — должен откатываться.
    app.var_fwhm.set(u"абв")
    app.var_dc_frac.set("-5")
    app.recompute()
    root.update()
    print(u"мусор в полях: FWHM -> %s, доля -> %s (откат сработал)"
          % (app.var_fwhm.get(), app.var_dc_frac.get()))

    app.var_band_hi.set("3.0")

    # --- пропускание по полосе: раскладка 3×2 со средним рядом дельты
    app.var_parseval.set(True)
    app.recompute()
    root.update()
    print(u"полоса счёта 0.2…3.0: панелей на фигуре %d" % len(app.fig.axes))
    shot("09_parseval_0p2_3p0")

    app.var_int_lo.set("0.3")
    app.var_int_hi.set("1.2")
    app.recompute()
    root.update()
    bright = [p for p in app.points if p.trace.angle == 10]
    dark = [p for p in app.points if p.trace.angle == 90]
    if bright and dark:
        print(u"  Δ на 0.3…1.2 ТГц: яркое %+0.3f дБ, гашение %+0.3f дБ"
              % (bright[0].delta_db, dark[0].delta_db))
    shot("10_parseval_0p3_1p2")

    app.var_window.set(True)
    app.var_fwhm.set("20")
    app.var_win_time.set(True)
    app.recompute()
    root.update()
    shot("11_window_in_time")
    app.var_window.set(False)
    app.var_win_time.set(False)

    app.full_band()
    root.update()
    if dark:
        d2 = [p for p in app.points if p.trace.angle == 90]
        print(u"  полная полоса 0…%s: Δ гашения %+.2e дБ — Парсеваль на экране"
              % (app.var_int_hi.get(), d2[0].delta_db))
    shot("12_full_band_parseval")

    app.var_parseval.set(False)
    app.var_int_lo.set("0.2")
    app.var_int_hi.set("3.0")
    app.recompute()
    root.update()
    print(u"выключение Парсеваля: панелей снова %d" % len(app.fig.axes))

    app.open_directory(PUREWAVE)
    root.update()
    print(u"каталог %s: %d точек" % (PUREWAVE, len(app.points)))
    shot("07_purewave")

    empty = tempfile.mkdtemp(prefix="tv_empty_")
    app.open_directory(empty)
    root.update()
    print(u"пустой каталог: %d точек, показано сообщение" % len(app.points))
    shot("08_empty_dir")

    # --- окно генерации дороги: параметры → предпросмотр → создание
    print(u"")
    print(u"окно генерации трека:")
    from track_viewer.gui import RoadDialog
    dlg = RoadDialog(root, app)
    root.update()

    target = tempfile.mkdtemp(prefix="tv_road_gui_")
    dlg.v["dir"].set(target)
    dlg.v["sample"].set(u"проверка_gui")
    dlg.do_preview()
    root.update()
    print(u"  предпросмотр на пустом каталоге: кнопка создания %s"
          % (u"активна" if str(dlg.btn_create["state"]) == "normal" else u"ЗАБЛОКИРОВАНА"))

    dlg.v["coarse_step"].set("20")
    print(u"  после правки параметра кнопка создания %s"
          % (u"активна" if str(dlg.btn_create["state"]) == "normal" else u"погасла"))

    dlg.do_preview()
    dlg.pv = dlg.pv
    from track_viewer.core import road
    created, skipped = road.apply(target, dlg.records, dlg.plan, dlg.pv)
    print(u"  создано %d пустышек, README.md и META.txt на месте: %s"
          % (created, os.path.isfile(os.path.join(target, "META.txt"))))

    dlg.do_preview()
    root.update()
    print(u"  повторный предпросмотр по тем же файлам: блокировка %s"
          % (u"есть" if dlg.pv.blocked else u"нет (все пустые — так и должно быть)"))

    with open(os.path.join(target, sorted(os.listdir(target))[1]), "w") as fh:
        fh.write("0.0\t1.0\n0.1\t2.0\n")
    dlg.do_preview()
    root.update()
    print(u"  после «съёмки» одной трассы: блокировка %s, кнопка %s"
          % (u"есть" if dlg.pv.blocked else u"НЕТ",
             u"погашена" if str(dlg.btn_create["state"]) == "disabled" else u"АКТИВНА"))
    dlg.win.destroy()

    root.destroy()
    print(u"")
    print(u"снимков: %d, каталог %s" % (len(shots), args.out))
    return 0


if __name__ == "__main__":
    sys.exit(main())
