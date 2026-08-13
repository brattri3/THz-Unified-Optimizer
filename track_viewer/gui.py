# -*- coding: utf-8 -*-
"""Графический интерфейс track_viewer: Tkinter + matplotlib.

Компоновка по ТЗ §4: слева панель управления, справа контент. Четыре панели
графиков сеткой 2×2 — верхняя пара угловая зависимость, нижняя спектр выбранного
измерения; в каждой паре слева линейная шкала 0…1, справа та же величина в дБ.

**Двойных осей Y нет нигде.** Линейная и логарифмическая — это две отдельные
панели, а не две шкалы на одном поле: совмещённые оси делают неотличимыми
«кривая изменилась» и «сменилась шкала».

Запуск из корня репозитория:
    .venv\\Scripts\\python.exe -m track_viewer.gui
    .venv\\Scripts\\python.exe -m track_viewer.gui data_pool/test_grid_33_11

Физики здесь нет: все числа приходят из `core.physics`, ровно те же, что печатает
`cli.py`. Это не стилистика, а условие того, чтобы панель данных и выгрузка CSV
не разошлись между собой.
"""
from __future__ import annotations

import os
import sys
import tkinter as tk
from tkinter import filedialog, messagebox, ttk

import numpy as np

import matplotlib
matplotlib.use("TkAgg")
from matplotlib.backends.backend_tkagg import (         # noqa: E402
    FigureCanvasTkAgg, NavigationToolbar2Tk)
from matplotlib.figure import Figure                    # noqa: E402
from matplotlib.gridspec import GridSpec                # noqa: E402

from .core import fit_malus as fm                        # noqa: E402
from .core import physics as ph                         # noqa: E402
from .core.scan import Scan                             # noqa: E402

# Палитра: повторы одного угла — одинаковый цвет, РАЗНАЯ форма маркера (O-7).
# Идентичность серии не кодируется одним лишь цветом (ТЗ §4.1).
COLOR_TRACK = "#1f77b4"
COLOR_WARN = "#d62728"
COLOR_SEL = "#ff7f0e"
COLOR_SPEC = "#2ca02c"
COLOR_FREQ = "#9467bd"      # пропускание по полосе — отдельный цвет И форма маркера
REP_MARKERS = ["o", "s", "^", "D", "v"]

FLOOR_DB = -60.0        # нижняя граница дБ-панелей: ниже уже шум установки

# Порядок в списке — порядок предпочтения: первым стоит метод по умолчанию.
NOISE_ORDER = [ph.NOISE_HF_TAIL, ph.NOISE_PRE_PULSE, ph.NOISE_DARK]
NOISE_LABELS = {
    ph.NOISE_HF_TAIL: u"ВЧ-хвост (Нафтали)",
    ph.NOISE_PRE_PULSE: u"предымпульс (завышает)",
    ph.NOISE_DARK: u"перекрытый пучок (dark)",
}
NOISE_BY_LABEL = dict((v, k) for k, v in NOISE_LABELS.items())

COLOR_FIT = "#8c564b"        # кривая фита T во времени
COLOR_FIT_BAND = "#e377c2"   # кривая фита T по полосе
COLOR_PAR = "#7f7f7f"        # тёмный канал |t_par|^2 на панелях каналов

WEIGHT_ORDER = [ph.FIT_RELATIVE, ph.FIT_NOISE, ph.FIT_UNIFORM]
WEIGHT_LABELS = {
    ph.FIT_RELATIVE: u"относительные (sigma ~ U)",
    ph.FIT_NOISE: u"по шумовому полу",
    ph.FIT_UNIFORM: u"равные (слепы к тени)",
}
WEIGHT_BY_LABEL = dict((v, k) for k, v in WEIGHT_LABELS.items())

ORDER_LABELS = {4: u"4 — точный (когерентный)", 2: u"2 — учебный, теряет фазу"}
ORDER_BY_LABEL = dict((v, k) for k, v in ORDER_LABELS.items())


class TrackViewer(object):

    def __init__(self, root, directory=None):
        self.root = root
        self.root.title(u"tds_track_viewer — контроль углового прогона THz-TDS")
        self.scan = None
        self.points = []
        self.index = 0
        self.fit = None

        self.var_dir = tk.StringVar(value=u"каталог не выбран")
        self.var_dc = tk.StringVar(value=u"предымпульс")
        self.var_dc_frac = tk.StringVar(value="15")
        self.var_window = tk.BooleanVar(value=False)
        self.var_fwhm = tk.StringVar(value="20")
        self.var_centre = tk.StringVar(value=u"свой пик")
        self.var_ref = tk.StringVar(value=u"ближайший ранний")
        self.var_band_lo = tk.StringVar(value="0.2")
        self.var_band_hi = tk.StringVar(value="3.0")
        self.var_parseval = tk.BooleanVar(value=False)
        self.var_int_lo = tk.StringVar(value="0.2")
        self.var_int_hi = tk.StringVar(value="3.0")
        self.var_win_time = tk.BooleanVar(value=False)
        self.var_noise = tk.StringVar(value=NOISE_LABELS[ph.NOISE_HF_TAIL])
        self.var_noise_hf = tk.StringVar(value="3.0")
        self.var_fit = tk.BooleanVar(value=False)
        self.var_fit_order = tk.StringVar(value=u"4 — точный (когерентный)")
        self.var_fit_weights = tk.StringVar(value=WEIGHT_LABELS[ph.FIT_RELATIVE])
        self.var_fit_spectral = tk.BooleanVar(value=False)
        self.var_fit_drop = tk.BooleanVar(value=False)
        self.var_status = tk.StringVar(value=u"")
        self.var_pos = tk.StringVar(value="0")

        # Опциональные ряды панелей (дельта и каналы фита) создаются по
        # необходимости; ключ помнит текущую раскладку, чтобы не пересоздавать
        # оси зря — пересборка сетки стоит дороже перерисовки.
        self._layout_key = None

        self._build()
        if directory:
            self.open_directory(directory)

    # ------------------------------------------------------------- построение
    def _build(self):
        outer = ttk.Frame(self.root, padding=6)
        outer.pack(fill=tk.BOTH, expand=True)

        left = ttk.Frame(outer, width=260)
        left.pack(side=tk.LEFT, fill=tk.Y, padx=(0, 8))
        left.pack_propagate(False)
        right = ttk.Frame(outer)
        right.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)

        self._build_controls(left)
        self._build_plots(right)
        self._build_databar(right)
        self._bind_keys()

    def _group(self, parent, title):
        f = ttk.LabelFrame(parent, text=title, padding=6)
        f.pack(fill=tk.X, pady=(0, 8))
        return f

    def _build_controls(self, parent):
        g = self._group(parent, u"ПАПКА")
        ttk.Button(g, text=u"выбрать…", command=self.choose_directory).pack(fill=tk.X)
        ttk.Label(g, textvariable=self.var_dir, wraplength=230,
                  foreground="#444").pack(fill=tk.X, pady=(4, 4))
        # Обе кнопки пересчитывают физику целиком. Разница ТОЛЬКО в том, берётся
        # ли содержимое файлов из кэша (ключ: имя, размер, mtime) или читается с
        # диска заново. Прежнее название «полный пересчёт» вводило в заблуждение:
        # звучало так, будто первая кнопка считает не всё.
        ttk.Button(g, text=u"ОБНОВИТЬ", command=self.refresh).pack(fill=tk.X)
        ttk.Button(g, text=u"перечитать все файлы с диска",
                   command=lambda: self.refresh(drop_cache=True)).pack(fill=tk.X, pady=(3, 0))
        ttk.Label(g, textvariable=self.var_status, wraplength=230).pack(fill=tk.X, pady=(4, 0))

        g = self._group(parent, u"РАСЧЁТ")
        ttk.Label(g, text=u"вычитание DC:").pack(anchor=tk.W)
        cb = ttk.Combobox(g, textvariable=self.var_dc, state="readonly", width=24,
                          values=[u"предымпульс", u"полное среднее", u"нет (как в ядре)"])
        cb.pack(fill=tk.X)
        cb.bind("<<ComboboxSelected>>", lambda e: self.recompute())

        row = ttk.Frame(g)
        row.pack(fill=tk.X, pady=(3, 6))
        ttk.Label(row, text=u"предымпульс, %:").pack(side=tk.LEFT)
        e = ttk.Entry(row, textvariable=self.var_dc_frac, width=6)
        e.pack(side=tk.RIGHT)
        e.bind("<Return>", lambda ev: self.recompute())

        ttk.Checkbutton(g, text=u"наложить гауссово окно", variable=self.var_window,
                        command=self.recompute).pack(anchor=tk.W)
        row = ttk.Frame(g)
        row.pack(fill=tk.X, pady=(2, 2))
        # Подпись «FWHM» стоит в интерфейсе явно: ширина гаусса задаётся
        # по-разному (σ, FWHM, 1/e), и без указания её читают наугад (ТЗ §3.3).
        ttk.Label(row, text=u"ширина FWHM, пс:").pack(side=tk.LEFT)
        e = ttk.Entry(row, textvariable=self.var_fwhm, width=6)
        e.pack(side=tk.RIGHT)
        e.bind("<Return>", lambda ev: self.recompute())
        ttk.Label(g, text=u"1000 пс = окно выключено", foreground="#777").pack(anchor=tk.W)

        ttk.Label(g, text=u"центр окна:").pack(anchor=tk.W, pady=(4, 0))
        cb = ttk.Combobox(g, textvariable=self.var_centre, state="readonly", width=24,
                          values=[u"свой пик", u"общий по пику референса"])
        cb.pack(fill=tk.X)
        cb.bind("<<ComboboxSelected>>", lambda e: self.recompute())

        ttk.Label(g, text=u"референс:").pack(anchor=tk.W, pady=(6, 0))
        for text in (u"ближайший ранний", u"ближайший поздний", u"среднее двух"):
            ttk.Radiobutton(g, text=text, value=text, variable=self.var_ref,
                            command=self.recompute).pack(anchor=tk.W)

        row = ttk.Frame(g)
        row.pack(fill=tk.X, pady=(6, 0))
        ttk.Label(row, text=u"полоса показа, ТГц:").pack(side=tk.LEFT)
        e2 = ttk.Entry(row, textvariable=self.var_band_hi, width=5)
        e2.pack(side=tk.RIGHT)
        e1 = ttk.Entry(row, textvariable=self.var_band_lo, width=5)
        e1.pack(side=tk.RIGHT, padx=(0, 3))
        for w in (e1, e2):
            w.bind("<Return>", lambda ev: self.recompute())

        ttk.Label(g, text=u"шумовой пол:").pack(anchor=tk.W, pady=(6, 0))
        cb = ttk.Combobox(g, textvariable=self.var_noise, state="readonly", width=24,
                          values=[NOISE_LABELS[k] for k in NOISE_ORDER])
        cb.pack(fill=tk.X)
        cb.bind("<<ComboboxSelected>>", lambda e: self._noise_changed())
        row = ttk.Frame(g)
        row.pack(fill=tk.X, pady=(2, 0))
        ttk.Label(row, text=u"шум по полосе выше, ТГц:").pack(side=tk.LEFT)
        self.entry_hf = ttk.Entry(row, textvariable=self.var_noise_hf, width=5)
        self.entry_hf.pack(side=tk.RIGHT)
        self.entry_hf.bind("<Return>", lambda ev: self.recompute())

        g = self._group(parent, u"ПРОПУСКАНИЕ ПО ПОЛОСЕ")
        ttk.Checkbutton(g, text=u"считать T по полосе (Парсеваль)",
                        variable=self.var_parseval,
                        command=self.recompute).pack(anchor=tk.W)
        row = ttk.Frame(g)
        row.pack(fill=tk.X, pady=(2, 0))
        # Полоса ИНТЕГРИРОВАНИЯ намеренно отдельна от полосы показа (решение
        # владельца): менять границы счёта, не трогая масштаб спектра, — это и
        # есть рабочий приём для отслеживания вклада шумовой части.
        ttk.Label(row, text=u"полоса счёта, ТГц:").pack(side=tk.LEFT)
        e2 = ttk.Entry(row, textvariable=self.var_int_hi, width=5)
        e2.pack(side=tk.RIGHT)
        e1 = ttk.Entry(row, textvariable=self.var_int_lo, width=5)
        e1.pack(side=tk.RIGHT, padx=(0, 3))
        for w in (e1, e2):
            w.bind("<Return>", lambda ev: self.recompute())
        ttk.Button(g, text=u"полная полоса — проверка Парсеваля",
                   command=self.full_band).pack(fill=tk.X, pady=(4, 0))
        ttk.Checkbutton(g, text=u"учитывать окно во временно́м T",
                        variable=self.var_win_time,
                        command=self.recompute).pack(anchor=tk.W, pady=(4, 0))
        ttk.Label(g, text=u"без этой галочки слева трасса не взвешена\n"
                          u"окном, а справа взвешена — Δ наберётся не\n"
                          u"только за счёт полосы",
                  foreground="#777", justify=tk.LEFT).pack(anchor=tk.W)

        g = self._group(parent, u"ФИТИРОВАНИЕ")
        ttk.Checkbutton(g, text=u"подгонять закон Малюса",
                        variable=self.var_fit,
                        command=self.recompute).pack(anchor=tk.W)
        ttk.Label(g, text=u"порядок базиса:").pack(anchor=tk.W, pady=(4, 0))
        cb = ttk.Combobox(g, textvariable=self.var_fit_order, state="readonly",
                          width=24, values=[ORDER_LABELS[4], ORDER_LABELS[2]])
        cb.pack(fill=tk.X)
        cb.bind("<<ComboboxSelected>>", lambda e: self.recompute())
        ttk.Label(g, text=u"веса МНК:").pack(anchor=tk.W, pady=(4, 0))
        cb = ttk.Combobox(g, textvariable=self.var_fit_weights, state="readonly",
                          width=24, values=[WEIGHT_LABELS[k] for k in WEIGHT_ORDER])
        cb.pack(fill=tk.X)
        cb.bind("<<ComboboxSelected>>", lambda e: self.recompute())
        # Предупреждение вынесено в интерфейс, а не только в гайд: равные веса
        # взвешивают точку как U², и зона гашения — единственное место, где
        # видна утечка, — просто не влияет на результат.
        ttk.Label(g, text=u"равные веса делают фит слепым к зоне\n"
                          u"гашения: у сессии A они дали |t_прод|²<0",
                  foreground="#777", justify=tk.LEFT).pack(anchor=tk.W)
        ttk.Checkbutton(g, text=u"фитировать каждую частоту",
                        variable=self.var_fit_spectral,
                        command=self.recompute).pack(anchor=tk.W, pady=(4, 0))
        ttk.Checkbutton(g, text=u"не брать точки с предупреждениями",
                        variable=self.var_fit_drop,
                        command=self.recompute).pack(anchor=tk.W)

        g = self._group(parent, u"ВЫГРУЗКА")
        ttk.Button(g, text=u"сохранить трек в CSV",
                   command=self.export_csv).pack(fill=tk.X)
        ttk.Label(g, text=u"для вёрстки в Origin Pro", foreground="#777").pack(anchor=tk.W)

        g = self._group(parent, u"ГЕНЕРАЦИЯ ТРЕКА")
        ttk.Button(g, text=u"разложить дорогу…",
                   command=self.open_road_dialog).pack(fill=tk.X)
        ttk.Label(g, text=u"пустые файлы под будущий прогон",
                  foreground="#777").pack(anchor=tk.W)

    def open_road_dialog(self):
        RoadDialog(self.root, self)

    def _build_plots(self, parent):
        self.fig = Figure(figsize=(10.5, 6.6), dpi=100)
        self._layout_axes(False)

        self.canvas = FigureCanvasTkAgg(self.fig, master=parent)
        self.canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)
        NavigationToolbar2Tk(self.canvas, parent).update()
        self.canvas.mpl_connect("button_press_event", self._on_click)

        bar = ttk.Frame(parent)
        bar.pack(fill=tk.X, pady=(4, 0))
        ttk.Label(bar, text=u"измерение №").pack(side=tk.LEFT)
        self.slider = ttk.Scale(bar, from_=0, to=0, orient=tk.HORIZONTAL,
                                command=self._on_slider)
        self.slider.pack(side=tk.LEFT, fill=tk.X, expand=True, padx=6)
        ttk.Button(bar, text=u"◀", width=3,
                   command=lambda: self.select(self.index - 1)).pack(side=tk.LEFT)
        ttk.Button(bar, text=u"▶", width=3,
                   command=lambda: self.select(self.index + 1)).pack(side=tk.LEFT)
        entry = ttk.Entry(bar, textvariable=self.var_pos, width=6)
        entry.pack(side=tk.LEFT, padx=(6, 0))
        entry.bind("<Return>", self._on_entry)

    def _layout_axes(self, with_delta, with_channels=False):
        """Пересобрать сетку панелей: угловые, дельта, спектр, каналы фита.

        Опциональные ряды не прячутся, а именно не создаются: пустая полоса
        съедала бы высоту на машине у прибора, где экран невелик. По той же
        причине оба ряда включаются независимо — при всём включённом на экране
        восемь панелей, и это уже предел читаемости.
        """
        key = (bool(with_delta), bool(with_channels))
        if self._layout_key == key:
            return
        self._layout_key = key
        self.fig.clear()

        rows = [1.0]                                   # угловые панели
        d_row = ch_row = None
        if with_delta:
            d_row = len(rows)
            rows.append(0.55)
        sp_row = len(rows)
        rows.append(1.0)
        if with_channels:
            ch_row = len(rows)
            rows.append(0.85)

        tight = len(rows) >= 4
        gs = GridSpec(len(rows), 2, figure=self.fig, height_ratios=rows,
                      left=0.07, right=0.985, top=0.95 if len(rows) > 2 else 0.94,
                      bottom=0.06 if tight else (0.07 if len(rows) > 2 else 0.09),
                      hspace=0.60 if tight else (0.55 if len(rows) > 2 else 0.42),
                      wspace=0.20)
        self.ax_tr_lin = self.fig.add_subplot(gs[0, 0])
        self.ax_tr_db = self.fig.add_subplot(gs[0, 1])
        self.ax_d_lin = self.fig.add_subplot(gs[d_row, 0]) if d_row else None
        self.ax_d_db = self.fig.add_subplot(gs[d_row, 1]) if d_row else None
        self.ax_sp_lin = self.fig.add_subplot(gs[sp_row, 0])
        self.ax_sp_db = self.fig.add_subplot(gs[sp_row, 1])
        self.ax_ch_lin = self.fig.add_subplot(gs[ch_row, 0]) if ch_row else None
        self.ax_ch_db = self.fig.add_subplot(gs[ch_row, 1]) if ch_row else None

    def _all_axes(self):
        axes = [self.ax_tr_lin, self.ax_tr_db, self.ax_sp_lin, self.ax_sp_db]
        return axes + [ax for ax in (self.ax_d_lin, self.ax_d_db,
                                     self.ax_ch_lin, self.ax_ch_db)
                       if ax is not None]

    def _build_databar(self, parent):
        f = ttk.LabelFrame(parent, text=u"измерение", padding=6)
        f.pack(fill=tk.X, pady=(6, 0))
        # Высота выросла с 11 строк: блок фита добавляет до семи. Полосы
        # прокрутки хватает, растягивать панель дальше нельзя — она отъедала бы
        # высоту у графиков.
        self.data_text = tk.Text(f, height=15, wrap=tk.WORD, relief=tk.FLAT,
                                 background="#f7f7f7")
        sb = ttk.Scrollbar(f, orient=tk.VERTICAL, command=self.data_text.yview)
        self.data_text.configure(yscrollcommand=sb.set)
        sb.pack(side=tk.RIGHT, fill=tk.Y)
        self.data_text.pack(fill=tk.X)
        self.data_text.configure(state=tk.DISABLED)
        self.data_text.tag_configure("warn", foreground=COLOR_WARN)

    def _bind_keys(self):
        self.root.bind("<Left>", lambda e: self.select(self.index - 1))
        self.root.bind("<Right>", lambda e: self.select(self.index + 1))
        self.root.bind("<Prior>", lambda e: self.select(self.index - 10))
        self.root.bind("<Next>", lambda e: self.select(self.index + 10))

    # ----------------------------------------------------------------- данные
    def settings(self):
        dc = {u"предымпульс": ph.DC_PRE_PULSE,
              u"полное среднее": ph.DC_FULL_MEAN,
              u"нет (как в ядре)": ph.DC_NONE}[self.var_dc.get()]
        ref = {u"ближайший ранний": ph.REF_EARLIER,
               u"ближайший поздний": ph.REF_LATER,
               u"среднее двух": ph.REF_AVERAGE}[self.var_ref.get()]
        centre = (ph.CENTER_OWN if self.var_centre.get() == u"свой пик"
                  else ph.CENTER_REF)
        return ph.Settings(
            dc_mode=dc,
            dc_fraction=self._num(self.var_dc_frac, 15.0, 0.5, 90.0) / 100.0,
            window_on=self.var_window.get(),
            window_fwhm_ps=self._num(self.var_fwhm, 20.0, 0.05, 1e6),
            window_center=centre,
            ref_mode=ref,
            band_lo=self._num(self.var_band_lo, 0.2, 0.0, 100.0),
            band_hi=self._num(self.var_band_hi, 3.0, 0.01, 100.0),
            parseval_on=self.var_parseval.get(),
            int_lo=self._num(self.var_int_lo, 0.2, 0.0, 1000.0),
            int_hi=self._num(self.var_int_hi, 3.0, 0.01, 1000.0),
            window_in_time=self.var_win_time.get(),
            noise_source=NOISE_BY_LABEL[self.var_noise.get()],
            noise_hf_from=self._num(self.var_noise_hf, 3.0, 0.05, 1000.0),
            fit_on=self.var_fit.get(),
            fit_order=ORDER_BY_LABEL[self.var_fit_order.get()],
            fit_weights=WEIGHT_BY_LABEL[self.var_fit_weights.get()],
            fit_spectral_on=self.var_fit_spectral.get(),
            fit_drop_warned=self.var_fit_drop.get(),
        )

    def _noise_changed(self):
        """Поле границы ВЧ-хвоста имеет смысл только для метода Нафтали."""
        active = NOISE_BY_LABEL[self.var_noise.get()] == ph.NOISE_HF_TAIL
        self.entry_hf.configure(state=(tk.NORMAL if active else tk.DISABLED))
        self.recompute()

    def nyquist(self):
        """Частота Найквиста текущего прогона: 1/(2·Δt), Δt читается из файла."""
        for p in self.points:
            t = p.trace.t
            if t is not None and len(t) > 1:
                return 0.5 / float(t[1] - t[0])
        return 100.0

    def full_band(self):
        """Полоса счёта → 0…Найквист. Тогда Δ обязана схлопнуться к ~10⁻⁶.

        Это не удобство, а живая проверка: равенство Парсеваля превращается из
        утверждения в наблюдаемое на экране число. Расхождение остаётся ровно
        потому, что временной интеграл берётся трапецией, а спектральный —
        прямоугольной суммой; вклад дают концы трассы.
        """
        self.var_parseval.set(True)
        self.var_int_lo.set("0")
        self.var_int_hi.set("%.4g" % (self.nyquist() * 1.001))
        self.recompute()

    def _num(self, var, default, lo, hi):
        """Числовое поле с откатом: неверный ввод не должен ронять пересчёт."""
        try:
            v = float(str(var.get()).replace(",", "."))
        except (TypeError, ValueError):
            var.set(str(default))
            return default
        if not (lo <= v <= hi):
            var.set(str(default))
            return default
        return v

    def choose_directory(self):
        d = filedialog.askdirectory(title=u"каталог углового прогона")
        if d:
            self.open_directory(d)

    def open_directory(self, directory):
        try:
            self.scan = Scan(directory)
        except ValueError as exc:
            messagebox.showerror(u"Каталог не открылся", u"%s" % exc)
            return
        self.var_dir.set(os.path.abspath(directory))
        self.index = 0
        self.recompute()

    def refresh(self, drop_cache=False):
        """Перечитать каталог: новые файлы добавляются, прежние берутся из кэша."""
        if self.scan is None:
            return
        keep = self.points[self.index].trace.seq if self.points else None
        self.scan.refresh(drop_cache=drop_cache)
        self.recompute()
        if keep is not None:
            for i, p in enumerate(self.points):
                if p.trace.seq == keep:
                    self.select(i)
                    break

    def recompute(self):
        if self.scan is None:
            return
        s = self.settings()
        self.points = ph.compute_track(self.scan, s)
        # Фит считается ЗДЕСЬ, один раз на пересчёт, а не в отрисовке: `select`
        # вызывается на каждый шаг слайдера, и подгонка по треку там была бы
        # лишней работой на каждое нажатие стрелки.
        self.fit = (fm.fit_track(self.points, s)
                    if s.fit_on and self.points else None)

        c = self.scan.counts()
        status = u"снято %d из %d трасс, референсов %d из %d" % (
            c["sig_done"], c["sig_total"], c["bg_done"], c["bg_total"])
        if c["dark_total"]:
            status += u", перекрытый пучок %d из %d" % (c["dark_done"],
                                                        c["dark_total"])
        if c["errors"]:
            status += u"\n⚠ файлов с ошибкой: %d" % c["errors"]
        nxt = self.scan.next_unmeasured()
        if nxt is not None:
            status += u"\nследующий: %s" % nxt.name
        self.var_status.set(status)

        n = len(self.points)
        self.slider.configure(to=max(0, n - 1))
        if n == 0:
            self._draw_empty()
            return
        self.index = min(self.index, n - 1)
        self.select(self.index)

    def select(self, index):
        if not self.points:
            return
        self.index = max(0, min(index, len(self.points) - 1))
        self.var_pos.set(str(self.points[self.index].trace.seq))
        try:
            self.slider.set(self.index)
        except tk.TclError:
            pass
        self._draw()
        self._fill_databar()

    def _on_slider(self, value):
        try:
            i = int(round(float(value)))
        except (TypeError, ValueError):
            return
        if self.points and i != self.index:
            self.select(i)

    def _on_entry(self, _event):
        """Поле ввода принимает СКВОЗНОЙ номер файла, а не позицию в списке.

        Оператор читает номер с имени файла (`047_a-94.txt`), а не считает,
        какой он по счёту среди снятых — пустышки нумерацию не занимают.
        """
        try:
            seq = int(str(self.var_pos.get()).strip())
        except ValueError:
            return
        for i, p in enumerate(self.points):
            if p.trace.seq == seq:
                self.select(i)
                return
        self.var_pos.set(str(self.points[self.index].trace.seq))

    def _on_click(self, event):
        """Клик по угловой панели выбирает ближайшую точку."""
        if event.inaxes not in (self.ax_tr_lin, self.ax_tr_db) or not self.points:
            return
        if event.xdata is None:
            return
        angles = np.array([p.trace.angle for p in self.points], dtype=float)
        self.select(int(np.argmin(np.abs(angles - event.xdata))))

    # ---------------------------------------------------------------- графики
    def _draw_empty(self):
        self._layout_axes(False)
        for ax in self._all_axes():
            ax.clear()
        self.ax_tr_lin.text(
            0.5, 0.5,
            u"Сигнальных трасс не найдено.\n\n"
            u"Инструмент опознаёт только журнальную нумерацию\n"
            u"NNN_a<угол>.txt и NNN_bg.txt — например 047_a-94.txt.\n"
            u"Каталоги старой схемы (только имена вида\n"
            u"<набор>_<угол>deg_repN_sig.txt) он не открывает.",
            ha="center", va="center", fontsize=10, color="#666",
            transform=self.ax_tr_lin.transAxes)
        for ax in self._all_axes():
            ax.set_xticks([])
            ax.set_yticks([])
        self.canvas.draw_idle()
        self._set_databar([u"каталог открыт, но трасс в нём нет"], warn_from=0)

    def _draw(self):
        sp = (self.fit is not None and self.fit.spectral is not None
              and self.fit.spectral.ch is not None)
        self._layout_axes(self.points[0].settings.parseval_on, sp)
        self._draw_track()
        if self.ax_d_lin is not None:
            self._draw_delta()
        self._draw_spectrum()
        if self.ax_ch_lin is not None:
            self._draw_channels()
        self.canvas.draw_idle()

    def _fit_grid(self):
        """Плотная сетка углов для непрерывной кривой фита."""
        ang = [p.trace.angle for p in self.points]
        return np.linspace(min(ang), max(ang), 721)

    def _draw_track(self):
        """Верхняя пара: угловая зависимость в линейной шкале и в дБ."""
        ax1, ax2 = self.ax_tr_lin, self.ax_tr_db
        ax1.clear()
        ax2.clear()

        sel = self.points[self.index]
        max_rep = max([p.trace.rep for p in self.points])

        for rep in range(1, max_rep + 1):
            grp = [p for p in self.points if p.trace.rep == rep]
            if not grp:
                continue
            marker = REP_MARKERS[(rep - 1) % len(REP_MARKERS)]
            ang = np.array([p.trace.angle for p in grp], dtype=float)
            T = np.array([p.T for p in grp], dtype=float)
            # Точка выше единицы прижимается к границе и красится в
            # предупреждающий цвет — но НЕ исчезает (ТЗ §4.1). Прецедент:
            # test_grid_40_20 дал 105.7 % из-за перепутанных sig/bg.
            over = T > 1.0
            label = u"rep%d" % rep if max_rep > 1 else u"измерения"
            ax1.plot(ang[~over], T[~over], marker, color=COLOR_TRACK,
                     markersize=5, linestyle="none", label=label)
            if np.any(over):
                ax1.plot(ang[over], np.ones(int(over.sum())), marker,
                         color=COLOR_WARN, markersize=7, linestyle="none",
                         markeredgecolor="black",
                         label=u"T > 100 %% (%d)" % int(over.sum()))
            with np.errstate(divide="ignore", invalid="ignore"):
                db = 10.0 * np.log10(np.where(T > 0, T, np.nan))
            ax2.plot(ang, db, marker, color=COLOR_TRACK, markersize=5,
                     linestyle="none", label=label)

            if sel.settings.parseval_on:
                # Второй ряд точек: то же пропускание, но набранное только в
                # полосе счёта. Форма маркера та же (это тот же повтор), цвет
                # другой — различается величина, а не серия.
                Tf = np.array([p.T_freq for p in grp], dtype=float)
                lab_f = (u"по полосе %.3g…%.3g ТГц"
                         % (sel.settings.int_lo, sel.settings.int_hi)) if rep == 1 else None
                ax1.plot(ang, np.minimum(Tf, 1.0), marker, color=COLOR_FREQ,
                         markersize=5, linestyle="none", markerfacecolor="none",
                         label=lab_f)
                with np.errstate(divide="ignore", invalid="ignore"):
                    dbf = 10.0 * np.log10(np.where(Tf > 0, Tf, np.nan))
                ax2.plot(ang, dbf, marker, color=COLOR_FREQ, markersize=5,
                         linestyle="none", markerfacecolor="none", label=lab_f)

        # Кривые фита: линией, а не маркером — различается величина, а не серия.
        if self.fit is not None and self.fit.ok:
            grid = self._fit_grid()
            curves = [(self.fit.time, COLOR_FIT, u"фит T во времени")]
            if self.fit.band is not None:
                curves.append((self.fit.band, COLOR_FIT_BAND, u"фит T по полосе"))
            for f, color, label in curves:
                y = f.curve(grid)
                ax1.plot(grid, np.minimum(y, 1.0), "-", color=color, linewidth=1.6,
                         label=label)
                with np.errstate(divide="ignore", invalid="ignore"):
                    ax2.plot(grid, 10.0 * np.log10(np.where(y > 0, y, np.nan)),
                             "-", color=color, linewidth=1.6, label=label)

        ax1.plot([sel.trace.angle], [min(sel.T, 1.0)], "o", markersize=12,
                 markerfacecolor="none", markeredgecolor=COLOR_SEL,
                 markeredgewidth=2, label=u"выбрано")
        if sel.T > 0:
            ax2.plot([sel.trace.angle], [10 * np.log10(sel.T)], "o", markersize=12,
                     markerfacecolor="none", markeredgecolor=COLOR_SEL,
                     markeredgewidth=2, label=u"выбрано")

        ax1.set_ylim(0.0, 1.0)          # шкала зафиксирована по требованию владельца
        ax1.set_xlabel(u"угол, °")
        ax1.set_ylabel(u"T, отн. ед.")
        ax1.set_title(u"пропускание по углу")
        ax2.set_xlabel(u"угол, °")
        ax2.set_ylabel(u"T, дБ")
        ax2.set_title(u"то же в логарифмической шкале")
        ax2.set_ylim(FLOOR_DB, 3.0)
        for ax in (ax1, ax2):
            ax.grid(True, alpha=0.3)
            ax.legend(loc="best", fontsize=8, framealpha=0.9)

    def _draw_delta(self):
        """Средний ряд: расхождение временно́го пропускания с полосовым.

        Обе панели идут с **автомасштабом**, а не с фиксированной шкалой 0…1,
        как угловые. Это осознанное отступление: разность на шкале пропускания
        легла бы плоской линией у нуля и не показывала бы ровно ничего, ради
        чего панель и заводилась. Линия нуля нарисована явно, чтобы автомасштаб
        не выдавал шум за сигнал.
        """
        ax1, ax2 = self.ax_d_lin, self.ax_d_db
        ax1.clear()
        ax2.clear()
        sel = self.points[self.index]
        max_rep = max([p.trace.rep for p in self.points])

        for rep in range(1, max_rep + 1):
            grp = [p for p in self.points if p.trace.rep == rep]
            if not grp:
                continue
            marker = REP_MARKERS[(rep - 1) % len(REP_MARKERS)]
            ang = np.array([p.trace.angle for p in grp], dtype=float)
            ax1.plot(ang, [100.0 * p.delta for p in grp], marker,
                     color=COLOR_FREQ, markersize=4, linestyle="none")
            ax2.plot(ang, [p.delta_db for p in grp], marker,
                     color=COLOR_FREQ, markersize=4, linestyle="none")

        # Непрерывная Δ — именно РАЗНОСТЬ ДВУХ ПОДОГНАННЫХ ФУНКЦИЙ, а не
        # подгонка к точкам Δ (решение владельца 2026-08-13). Разница
        # существенна: фит каждой величины делается по своим весам, и их
        # разность честнее, чем подгонка к уже вычтенной разности.
        if self.fit is not None and self.fit.ok and self.fit.band is not None:
            grid = self._fit_grid()
            d = self.fit.delta_curve(grid)
            ax1.plot(grid, 100.0 * d, "-", color=COLOR_FIT, linewidth=1.4,
                     label=u"разность фитов")
            ft, fb = self.fit.time.curve(grid), self.fit.band.curve(grid)
            with np.errstate(divide="ignore", invalid="ignore"):
                ratio = np.where((ft > 0) & (fb > 0), ft / fb, np.nan)
                ax2.plot(grid, 10.0 * np.log10(ratio), "-", color=COLOR_FIT,
                         linewidth=1.4, label=u"разность фитов")
            for ax in (ax1, ax2):
                ax.legend(loc="best", fontsize=7, framealpha=0.9)

        for ax in (ax1, ax2):
            ax.axhline(0.0, color="#888", linewidth=0.8)
            ax.axvline(sel.trace.angle, color=COLOR_SEL, linewidth=1.0, alpha=0.6)
            ax.grid(True, alpha=0.3)
            ax.set_xlabel(u"угол, °")
        ax1.set_ylabel(u"Δ, п.п.")
        ax2.set_ylabel(u"Δ, дБ")
        ax1.set_title(u"Δ = T во времени − T по полосе", fontsize=9)
        ax2.set_title(u"то же в дБ: 10·lg(T время / T полоса)", fontsize=9)

    def _draw_spectrum(self):
        """Нижняя пара: спектральное пропускание выбранного измерения."""
        ax1, ax2 = self.ax_sp_lin, self.ax_sp_db
        ax1.clear()
        ax2.clear()
        p = self.points[self.index]

        if p.freqs is None or p.freqs.size == 0:
            for ax in (ax1, ax2):
                ax.text(0.5, 0.5, u"спектр не построен:\nнет точек выше порога",
                        ha="center", va="center", color="#666",
                        transform=ax.transAxes)
        else:
            ax1.plot(p.freqs, p.T_nu, "-", color=COLOR_SPEC, linewidth=1.2,
                     label=u"T(ν)")
            with np.errstate(divide="ignore", invalid="ignore"):
                db = 10.0 * np.log10(np.where(p.T_nu > 0, p.T_nu, np.nan))
            ax2.plot(p.freqs, db, "-", color=COLOR_SPEC, linewidth=1.2,
                     label=u"T(ν)")
            ax2.set_ylim(FLOOR_DB, 3.0)
            if p.T_noise is not None:
                # Где зелёная кривая ложится на серую — там уже нечего измерять,
                # и «структура» спектра выше этой линии читаться не должна.
                with np.errstate(divide="ignore", invalid="ignore"):
                    ndb = 10.0 * np.log10(np.where(p.T_noise > 0, p.T_noise, np.nan))
                ax1.plot(p.freqs, p.T_noise, "--", color="#888", linewidth=1.0,
                         label=u"шумовой пол")
                ax2.plot(p.freqs, ndb, "--", color="#888", linewidth=1.0,
                         label=u"шумовой пол")
            sp = self.fit.spectral if self.fit is not None else None
            if sp is not None:
                # Восстановленная фитом кривая СГЛАЖЕНА по построению: каждая её
                # точка получена по всем трассам трека сразу, а не по одной паре
                # образец/референс. Расхождение с измеренной — мера шума, а не
                # признак ошибки.
                Tfit = sp.transmission(p.trace.angle)
                ax1.plot(sp.freqs, Tfit, "-", color=COLOR_FIT, linewidth=1.4,
                         label=u"восстановлено фитом")
                with np.errstate(divide="ignore", invalid="ignore"):
                    ax2.plot(sp.freqs,
                             10.0 * np.log10(np.where(Tfit > 0, Tfit, np.nan)),
                             "-", color=COLOR_FIT, linewidth=1.4,
                             label=u"восстановлено фитом")
            if p.settings.parseval_on:
                # Полоса счёта показана заливкой прямо на спектре: иначе её
                # границы существуют только в полях ввода, и связь между
                # числом Δ и участком кривой приходится держать в голове.
                for ax in (ax1, ax2):
                    ax.axvspan(p.settings.int_lo, p.settings.int_hi,
                               color=COLOR_FREQ, alpha=0.10, zorder=0,
                               label=u"полоса счёта")
            for ax in (ax1, ax2):
                ax.legend(loc="best", fontsize=8, framealpha=0.9)

        ax1.set_ylim(0.0, 1.0)
        ax1.set_xlabel(u"частота, ТГц")
        ax1.set_ylabel(u"T(ν), отн. ед.")
        ax1.set_title(u"спектр: %s" % p.trace.name)
        ax2.set_xlabel(u"частота, ТГц")
        ax2.set_ylabel(u"T(ν), дБ")
        ax2.set_title(u"то же в логарифмической шкале")
        for ax in (ax1, ax2):
            ax.grid(True, alpha=0.3)

    def _draw_channels(self):
        """Свой ряд: спектры каналов |t⊥|²(ν) и |t∥|²(ν) из побинного фита.

        Это и есть тот график пропускания, который из одной трассы получить
        нельзя: каждая точка кривой — результат подгонки по всем углам на
        данной частоте. Никакой модели частотной зависимости при этом не
        используется, бины независимы. Расстояние между кривыми — экстинкция.
        """
        ax1, ax2 = self.ax_ch_lin, self.ax_ch_db
        ax1.clear()
        ax2.clear()
        sp = self.fit.spectral
        ch = sp.ch
        perp, par = ch["I_perp"], ch["I_par"]

        ax1.plot(sp.freqs, perp, "-", color=COLOR_SPEC, linewidth=1.4,
                 label=u"|t⊥|² — яркий канал")
        # Отрицательные значения тёмного канала физически невозможны и означают,
        # что в тени фит упёрся в шум. Их не прячем — иначе кривая соврёт.
        neg = par <= 0
        ax1.plot(sp.freqs[~neg], par[~neg], "-", color=COLOR_PAR, linewidth=1.4,
                 label=u"|t∥|² — утечка")
        if np.any(neg):
            ax1.plot(sp.freqs[neg], np.zeros(int(neg.sum())), "x",
                     color=COLOR_WARN, markersize=5, linestyle="none",
                     label=u"|t∥|² < 0 (%d)" % int(neg.sum()))
        with np.errstate(divide="ignore", invalid="ignore"):
            for y, color, label in ((perp, COLOR_SPEC, u"|t⊥|²"),
                                    (par, COLOR_PAR, u"|t∥|²")):
                ax2.plot(sp.freqs, 10.0 * np.log10(np.where(y > 0, y, np.nan)),
                         "-", color=color, linewidth=1.4, label=label)

        ax1.set_ylim(0.0, 1.0)
        ax2.set_ylim(FLOOR_DB, 3.0)
        ax1.set_ylabel(u"каналы, отн. ед.")
        ax2.set_ylabel(u"каналы, дБ")
        ax1.set_title(u"каналы по частоте (фит по всем углам)", fontsize=9)
        ax2.set_title(u"расстояние между кривыми = экстинкция(ν)", fontsize=9)
        for ax in (ax1, ax2):
            ax.set_xlabel(u"частота, ТГц")
            ax.grid(True, alpha=0.3)
            ax.legend(loc="best", fontsize=7, framealpha=0.9)

    # ------------------------------------------------------------ панель данных
    def _fill_databar(self):
        p = self.points[self.index]
        lines = list(p.summary_lines())
        lines.insert(1, u"измерение %d из %d в этом треке"
                     % (self.index + 1, len(self.points)))
        if self.fit is not None:
            # Строки фита вставляются ДО расчёта warn_from: иначе они попали бы
            # в красный хвост предупреждений, которым не являются.
            extra = list(self.fit.summary_lines())
            if self.fit.ok:
                line = self.fit.point_line(self.index)
                if line:
                    extra.insert(1, u"  " + line)
            lines[len(lines) - len(p.warnings):len(lines) - len(p.warnings)] = extra
        warn_from = len(lines) - len(p.warnings)
        for tr in self.scan.errors:
            lines.append(u"⚠ файл %s не разобран: %s" % (tr.name, tr.error))
        self._set_databar(lines, warn_from)

    def _set_databar(self, lines, warn_from):
        self.data_text.configure(state=tk.NORMAL)
        self.data_text.delete("1.0", tk.END)
        for i, line in enumerate(lines):
            self.data_text.insert(tk.END, line + u"\n",
                                  ("warn",) if i >= warn_from else ())
        self.data_text.configure(state=tk.DISABLED)

    # ---------------------------------------------------------------- выгрузка
    def export_csv(self):
        if not self.points:
            messagebox.showinfo(u"Нечего выгружать", u"Трек пуст.")
            return
        from .cli import export_csv
        default = os.path.basename(os.path.normpath(self.scan.directory)) + u"_track.csv"
        path = filedialog.asksaveasfilename(
            title=u"сохранить трек", defaultextension=".csv",
            initialfile=default, filetypes=[(u"CSV", "*.csv")])
        if not path:
            return
        try:
            export_csv(self.points, path, self.fit)
        except (IOError, OSError) as exc:
            messagebox.showerror(u"Не сохранилось", u"%s" % exc)
            return
        messagebox.showinfo(u"Сохранено", u"%d точек:\n%s" % (len(self.points), path))


class RoadDialog(object):
    """Окно генерации дороги: параметры → предпросмотр → создание.

    Создание физически невозможно до предпросмотра: кнопка «СОЗДАТЬ ФАЙЛЫ»
    включается только после успешного просмотра и гаснет обратно при любом
    изменении параметров. Это пункт 4 гардрейлов ТЗ §5.2, вынесенный в
    состояние виджета, а не оставленный на дисциплину оператора.
    """

    def __init__(self, parent, app):
        self.app = app
        self.records = None
        self.plan = None
        self.pv = None

        self.win = tk.Toplevel(parent)
        self.win.title(u"Генерация трека — пустые файлы будущего прогона")
        self.win.transient(parent)
        self.win.geometry("880x680")

        self.v = {}
        for name, value in (
                ("sample", u"образец"), ("dir", u""),
                ("coarse_from", "-70"), ("coarse_to", "70"), ("coarse_step", "10"),
                ("fine_from", "80"), ("fine_to", "100"), ("fine_step", "2"),
                ("anchor", "0"), ("anchors_rep2", "40, 0, -40"),
                ("bg_every", "3"), ("repeats", "1"), ("dark_every", "20")):
            self.v[name] = tk.StringVar(value=value)
        self.v_fine = tk.BooleanVar(value=True)
        self.v_mirror = tk.BooleanVar(value=True)
        self.v_repeat = tk.BooleanVar(value=True)
        self.v_bg = tk.StringVar(value=u"в начале, каждые n и в конце")
        self.v_order = tk.StringVar(value=u"блочный")
        self.v_dark = tk.StringVar(value=u"в начале и в конце")

        self._build()

        # Гардрейл «сначала предпросмотр» держится на изменении САМИХ переменных,
        # а не на событии клавиатуры: вставка мышью, выбор из списка и правка из
        # кода тоже обязаны гасить кнопку создания. Привязка к <KeyRelease> это
        # пропускала.
        for var in list(self.v.values()) + [self.v_fine, self.v_mirror,
                                            self.v_repeat, self.v_bg, self.v_order,
                                            self.v_dark]:
            var.trace_add("write", lambda *a: self._invalidate())

    def _build(self):
        top = ttk.Frame(self.win, padding=8)
        top.pack(fill=tk.BOTH, expand=True)

        left = ttk.Frame(top, width=330)
        left.pack(side=tk.LEFT, fill=tk.Y, padx=(0, 8))
        left.pack_propagate(False)
        right = ttk.Frame(top)
        right.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)

        def row(parent, label, key, width=7):
            f = ttk.Frame(parent)
            f.pack(fill=tk.X, pady=1)
            ttk.Label(f, text=label).pack(side=tk.LEFT)
            e = ttk.Entry(f, textvariable=self.v[key], width=width)
            e.pack(side=tk.RIGHT)
            return e

        g = ttk.LabelFrame(left, text=u"КУДА", padding=6)
        g.pack(fill=tk.X, pady=(0, 6))
        ttk.Button(g, text=u"выбрать каталог…",
                   command=self._choose).pack(fill=tk.X)
        ttk.Label(g, textvariable=self.v["dir"], wraplength=300,
                  foreground="#444").pack(fill=tk.X, pady=(3, 0))
        row(g, u"имя образца:", "sample", 18)

        g = ttk.LabelFrame(left, text=u"ОСНОВНОЙ ДИАПАЗОН", padding=6)
        g.pack(fill=tk.X, pady=(0, 6))
        row(g, u"от, °:", "coarse_from")
        row(g, u"до, °:", "coarse_to")
        row(g, u"шаг, °:", "coarse_step")
        row(g, u"яркий якорь, °:", "anchor")

        g = ttk.LabelFrame(left, text=u"ВТОРОЙ ДИАПАЗОН (зона гашения)", padding=6)
        g.pack(fill=tk.X, pady=(0, 6))
        ttk.Checkbutton(g, text=u"включить", variable=self.v_fine,
                        command=self._invalidate).pack(anchor=tk.W)
        row(g, u"от, °:", "fine_from")
        row(g, u"до, °:", "fine_to")
        row(g, u"шаг, °:", "fine_step")
        ttk.Checkbutton(g, text=u"зеркальная тень (контроль π-периодичности)",
                        variable=self.v_mirror,
                        command=self._invalidate).pack(anchor=tk.W)
        ttk.Checkbutton(g, text=u"повтор обратным ходом (мера дрейфа)",
                        variable=self.v_repeat,
                        command=self._invalidate).pack(anchor=tk.W)

        g = ttk.LabelFrame(left, text=u"РЕФЕРЕНСЫ И ПОРЯДОК", padding=6)
        g.pack(fill=tk.X, pady=(0, 6))
        cb = ttk.Combobox(g, textvariable=self.v_bg, state="readonly", width=30,
                          values=[u"один в начале", u"в начале и в конце",
                                  u"в начале, каждые n и в конце"])
        cb.pack(fill=tk.X)
        cb.bind("<<ComboboxSelected>>", lambda e: self._invalidate())
        row(g, u"n =", "bg_every")
        cb = ttk.Combobox(g, textvariable=self.v_order, state="readonly", width=30,
                          values=[u"блочный", u"по возрастанию"])
        cb.pack(fill=tk.X, pady=(4, 0))
        cb.bind("<<ComboboxSelected>>", lambda e: self._invalidate())
        row(g, u"якоря rep2, °:", "anchors_rep2", 14)
        row(g, u"повторов плана:", "repeats")

        g = ttk.LabelFrame(left, text=u"ПЕРЕКРЫТЫЙ ПУЧОК (dark)", padding=6)
        g.pack(fill=tk.X, pady=(0, 6))
        cb = ttk.Combobox(g, textvariable=self.v_dark, state="readonly", width=30,
                          values=[u"не снимать", u"один в начале",
                                  u"в начале и в конце",
                                  u"в начале, каждые n и в конце"])
        cb.pack(fill=tk.X)
        cb.bind("<<ComboboxSelected>>", lambda e: self._invalidate())
        row(g, u"n =", "dark_every")
        ttk.Label(g, text=u"NNN_dark.txt — непрозрачная преграда в пучке.\n"
                          u"Шум закрытого тракта: даёт динамический\n"
                          u"диапазон независимо от метода Нафтали.\n"
                          u"Из трасс НЕ вычитается.",
                  foreground="#777", justify=tk.LEFT).pack(anchor=tk.W, pady=(3, 0))

        bar = ttk.Frame(left)
        bar.pack(fill=tk.X, pady=(4, 0))
        ttk.Button(bar, text=u"ПРЕДПРОСМОТР",
                   command=self.do_preview).pack(fill=tk.X)
        self.btn_create = ttk.Button(bar, text=u"СОЗДАТЬ ФАЙЛЫ",
                                     command=self.do_create, state=tk.DISABLED)
        self.btn_create.pack(fill=tk.X, pady=(4, 0))
        ttk.Label(bar, text=u"создаются только пустые файлы;\n"
                            u"существующие не перезаписываются никогда",
                  foreground="#777", justify=tk.LEFT).pack(anchor=tk.W, pady=(4, 0))

        self.text = tk.Text(right, wrap=tk.NONE, relief=tk.FLAT,
                            background="#f7f7f7")
        sb = ttk.Scrollbar(right, orient=tk.VERTICAL, command=self.text.yview)
        self.text.configure(yscrollcommand=sb.set)
        sb.pack(side=tk.RIGHT, fill=tk.Y)
        self.text.pack(fill=tk.BOTH, expand=True)
        self.text.tag_configure("warn", foreground=COLOR_WARN)
        self._show([u"Задайте параметры и нажмите ПРЕДПРОСМОТР.",
                    u"",
                    u"Значения по умолчанию повторяют проверенный план",
                    u"test_grid_33_11: яркий якорь, тень 80…100° с шагом 2°,",
                    u"яркая зона ±70° с шагом 10°, зеркальная тень и её повтор",
                    u"обратным ходом, референс каждые 3 трассы."], 99)

    def _choose(self):
        d = filedialog.askdirectory(title=u"каталог будущего прогона")
        if d:
            self.v["dir"].set(d)
            self._invalidate()

    def _invalidate(self):
        """Любая правка параметров гасит кнопку создания — предпросмотр устарел."""
        self.btn_create.configure(state=tk.DISABLED)
        self.pv = None

    def _show(self, lines, warn_from):
        self.text.configure(state=tk.NORMAL)
        self.text.delete("1.0", tk.END)
        for i, line in enumerate(lines):
            self.text.insert(tk.END, line + u"\n",
                             ("warn",) if i >= warn_from else ())
        self.text.configure(state=tk.DISABLED)

    def _int(self, key, default):
        try:
            return int(float(str(self.v[key].get()).replace(",", ".")))
        except (TypeError, ValueError):
            self.v[key].set(str(default))
            return default

    def build(self):
        from .core import road
        anchors = []
        for chunk in str(self.v["anchors_rep2"].get()).replace(";", ",").split(","):
            chunk = chunk.strip()
            if chunk:
                try:
                    anchors.append(int(float(chunk)))
                except ValueError:
                    pass
        plan = road.Plan(
            sample=str(self.v["sample"].get()).strip() or u"образец",
            coarse_from=self._int("coarse_from", -70),
            coarse_to=self._int("coarse_to", 70),
            coarse_step=self._int("coarse_step", 10),
            fine_on=self.v_fine.get(),
            fine_from=self._int("fine_from", 80),
            fine_to=self._int("fine_to", 100),
            fine_step=self._int("fine_step", 2),
            mirror_fine=self.v_mirror.get(),
            repeat_fine=self.v_repeat.get(),
            anchor=self._int("anchor", 0),
            anchors_rep2=tuple(anchors),
            bg_mode={u"один в начале": road.BG_START,
                     u"в начале и в конце": road.BG_START_END,
                     u"в начале, каждые n и в конце": road.BG_EVERY_N}[self.v_bg.get()],
            bg_every=self._int("bg_every", 3),
            order=(road.ORDER_BLOCK if self.v_order.get() == u"блочный"
                   else road.ORDER_ASCENDING),
            repeats=self._int("repeats", 1),
            dark_mode={u"не снимать": road.DARK_NONE,
                       u"один в начале": road.DARK_START,
                       u"в начале и в конце": road.DARK_START_END,
                       u"в начале, каждые n и в конце": road.DARK_EVERY_N}[
                           self.v_dark.get()],
            dark_every=self._int("dark_every", 20),
        )
        return plan, road.build_plan(plan)

    def do_preview(self):
        from .core import road
        directory = str(self.v["dir"].get()).strip()
        if not directory:
            messagebox.showinfo(u"Нужен каталог", u"Выберите каталог будущего прогона.")
            return
        try:
            self.plan, self.records = self.build()
        except ValueError as exc:
            self._show([u"План неисполним:", u"", u"%s" % exc], 0)
            self._invalidate()
            return

        self.pv = road.preview(directory, self.records)
        lines = list(self.pv.lines())
        warn_from = len(lines)
        if self.pv.blocked:
            warn_from = next(i for i, s in enumerate(lines) if s.startswith(u"ОТКАЗ"))
        else:
            lines.append(u"")
            lines.append(u"Порядок съёмки (следующий пустой файл = следующее действие):")
            lines.append(u"")
            for rec in self.records:
                if rec["kind"] == "bg":
                    what = u"фон (пустой пучок)"
                elif rec["kind"] == "dark":
                    what = u"ПЕРЕКРЫТЬ ПУЧОК непрозрачной преградой"
                else:
                    what = u"угол %+d°%s" % (rec["angle"],
                                             u"   повтор %d" % rec["rep"]
                                             if rec["rep"] and rec["rep"] > 1
                                             else u"")
                lines.append(u"  %-16s %s" % (road.fname(rec), what))
            warn_from = len(lines)
        self._show(lines, warn_from)
        self.btn_create.configure(
            state=tk.DISABLED if self.pv.blocked else tk.NORMAL)

    def do_create(self):
        from .core import road
        if self.pv is None or self.pv.blocked or not self.records:
            return
        directory = str(self.v["dir"].get()).strip()
        n = len(self.pv.to_create)
        if not messagebox.askyesno(
                u"Создать файлы?",
                u"Будет создано %d пустых файлов в\n%s\n\n"
                u"Существующие файлы не изменяются." % (n, directory)):
            return

        # Свежий предпросмотр прямо перед записью: между показом и нажатием
        # прибор мог дописать трассу в этот каталог. Правило «непустой файл
        # отменяет операцию целиком» обязано действовать и в этом окне.
        fresh = road.preview(directory, self.records)
        if fresh.blocked:
            self.pv = fresh
            self._show(fresh.lines(),
                       next(i for i, s in enumerate(fresh.lines())
                            if s.startswith(u"ОТКАЗ")))
            self.btn_create.configure(state=tk.DISABLED)
            messagebox.showerror(
                u"Не создано",
                u"Каталог изменился с момента предпросмотра: появилось %d непустых "
                u"файлов с плановыми именами. Не создано ничего."
                % len(fresh.conflicts))
            return

        try:
            created, skipped = road.apply(directory, self.records, self.plan, fresh)
        except (ValueError, IOError, OSError) as exc:
            messagebox.showerror(u"Не создано", u"%s" % exc)
            return
        messagebox.showinfo(
            u"Дорога разложена",
            u"создано пустых файлов: %d\nуже были: %d\n\nв каталоге также META.txt "
            u"и README.md с планом" % (created, skipped))
        self._invalidate()
        self.app.open_directory(directory)
        self.win.destroy()


def main(argv=None):
    argv = sys.argv[1:] if argv is None else argv
    directory = argv[0] if argv else None
    root = tk.Tk()
    root.geometry("1320x900")
    TrackViewer(root, directory)
    root.mainloop()
    return 0


if __name__ == "__main__":
    sys.exit(main())
