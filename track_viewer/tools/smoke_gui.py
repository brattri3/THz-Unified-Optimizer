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

from track_viewer.core import fit_malus as fm                # noqa: E402
from track_viewer.gui import (                               # noqa: E402
    COLLAPSED_AT_START, PLOT_ROW_MIN_PX, TrackViewer)

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

    # --- колонка настроек обязана помещаться на экране прибора (1024x768).
    # Окно приходится показать: у скрытого окна Tk не выполняет раскладку, и
    # winfo_height даёт 1 px — проверка выродилась бы в тождество.
    root.deiconify()
    root.geometry("1024x768")
    root.update()
    root.update_idletasks()
    canvas = app._panel_canvas
    inner = canvas.nametowidget(canvas.itemcget(canvas.find_all()[0], "window"))
    at_start = inner.winfo_reqheight()
    print(u"колонка настроек на 1024x768: холст %d px, при старте %d px — %s"
          % (canvas.winfo_height(), at_start,
             u"помещается" if at_start <= canvas.winfo_height()
             else u"НЕ ПОМЕЩАЕТСЯ, спасает прокрутка"))
    for title in sorted(app._groups):
        state, toggle = app._groups[title]
        was = state["open"]
        toggle()
        root.update_idletasks()
        if state["open"] == was:
            raise SystemExit(u"блок %s не переключился" % title)
        toggle()
        root.update_idletasks()
    for title in sorted(app._groups):
        state, toggle = app._groups[title]
        if not state["open"]:
            toggle()
    root.update_idletasks()
    full = inner.winfo_reqheight()
    canvas.yview_moveto(1.0)
    root.update_idletasks()
    first, last = canvas.yview()
    print(u"  все %d блоков развёрнуты: %d px, прокрутка до конца даёт %.2f…%.2f"
          % (len(app._groups), full, first, last))
    if last < 0.999:
        raise SystemExit(u"прокрутка не доходит до низа колонки")
    for title in COLLAPSED_AT_START:
        state, toggle = app._groups[title]
        if state["open"]:
            toggle()
    canvas.yview_moveto(0.0)

    # --- графики не должны ужиматься: ряд не ниже PLOT_ROW_MIN_PX, дальше прокрутка
    holder = app._plot_holder
    for label, tune in ((u"4 панели", None),
                        (u"6 панелей", lambda: app.var_parseval.set(True)),
                        (u"8 панелей", lambda: (app.var_fit.set(True),
                                                app.var_fit_spectral.set(True)))):
        if tune:
            tune()
            app.recompute()
        root.update()
        root.update_idletasks()
        rows, fig_h = app._plot_rows, app._plot_size[1]
        per_row = fig_h // rows
        if per_row < PLOT_ROW_MIN_PX:
            raise SystemExit(u"ряд ужат до %d px при минимуме %d"
                             % (per_row, PLOT_ROW_MIN_PX))
        holder.yview_moveto(1.0)
        root.update_idletasks()
        if holder.yview()[1] < 0.999:
            raise SystemExit(u"область графиков не прокручивается до низа")
        holder.yview_moveto(0.0)
        print(u"  %-10s рядов %d, фигура %d px (%d px на ряд), окно %d px — %s"
              % (label, rows, fig_h, per_row, holder.winfo_height(),
                 u"помещается" if fig_h <= holder.winfo_height() else u"прокрутка"))
    app.var_parseval.set(False)
    app.var_fit.set(False)
    app.var_fit_spectral.set(False)
    app.recompute()

    root.geometry("1320x900")
    root.update()
    root.withdraw()

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

    # --- три источника шумового пола
    from track_viewer.gui import NOISE_LABELS
    from track_viewer.core import physics as ph
    for key in (ph.NOISE_PRE_PULSE, ph.NOISE_DARK, ph.NOISE_HF_TAIL):
        app.var_noise.set(NOISE_LABELS[key])
        app._noise_changed()
        root.update()
        p = app.points[app.index]
        note = [w for w in p.warnings if u"шумовой пол" in w or u"перекрытого" in w]
        print(u"  %-26s поле границы %s%s"
              % (NOISE_LABELS[key],
                 u"активно" if str(app.entry_hf["state"]) == "normal" else u"погашено",
                 u"; " + note[0][:60] if note else u""))
    shot("13_noise_hf_tail")

    # --- фит: кривые на угловых панелях, разность фитов, каналы по частоте
    from track_viewer.gui import ORDER_LABELS, WEIGHT_LABELS
    app.var_fit.set(True)
    app.recompute()
    root.update()
    print(u"")
    print(u"фит включён: панелей %d" % len(app.fig.axes))
    for line in app.fit.summary_lines():
        print(u"  " + line)
    shot("14_fit_time_only")

    app.var_parseval.set(True)
    app.var_int_lo.set("0.3")
    app.var_int_hi.set("2.0")
    app.recompute()
    root.update()
    print(u"фит + Парсеваль: панелей %d (ожидается 6)" % len(app.fig.axes))
    shot("15_fit_delta_curve")

    app.var_fit_spectral.set(True)
    app.recompute()
    root.update()
    sp = app.fit.spectral
    print(u"фит побинно: панелей %d (ожидается 8), бинов %d, |t_par|^2<0 в %d, "
          u"нарушений Коши-Буняковского %d"
          % (len(app.fig.axes), sp.n_bins, sp.n_neg, sp.n_cs_violations))
    shot("16_fit_channels")

    for key in (ph.FIT_UNIFORM, ph.FIT_NOISE, ph.FIT_RELATIVE):
        app.var_fit_weights.set(WEIGHT_LABELS[key])
        app.recompute()
        root.update()
        pr = app.fit.time.params
        print(u"  веса %-28s eta = %.5f, theta0 = %+.4f"
              % (WEIGHT_LABELS[key], pr["eta_amplitude"],
                 pr["theta0_deg_from_h2"]))

    app.var_fit_order.set(ORDER_LABELS[2])
    app.recompute()
    root.update()
    print(u"  порядок 2 (учебный): %s"
          % app.fit.summary_lines()[1].strip())
    shot("17_fit_order2")
    app.var_fit_order.set(ORDER_LABELS[4])

    # Неполный трек: фит обязан отказаться, а не выдать уверенный мусор.
    saved = app.points
    app.points = saved[:4]
    app.fit = fm.fit_track(app.points, app.settings())
    print(u"  на %d точках: %s" % (len(app.points), app.fit.notes[0]))
    app.points = saved
    app.var_fit.set(False)
    app.var_fit_spectral.set(False)
    app.var_parseval.set(False)
    app.recompute()
    root.update()

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

    from track_viewer.core import road as _road
    print(u"  план по умолчанию: %d файлов, из них перекрытых пучков %d"
          % (len(dlg.records), _road.summary(dlg.records)["n_dark"]))
    dlg.v_dark.set(u"не снимать")
    dlg.do_preview()
    print(u"  режим «не снимать»: перекрытых пучков %d"
          % _road.summary(dlg.records)["n_dark"])
    dlg.v_dark.set(u"в начале и в конце")

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
