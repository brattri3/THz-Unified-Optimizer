"""tkinter GUI поверх `service_calc.py` -- обслуживание аттенюатора в THz-TDS
спектрометре (задача C9, санкция владельца 2026-08-19). Один сценарий:
двунаправленный калькулятор (дБ<->угол), относительный/абсолютный режим,
выбираемая метрика полосы, без доверительного интервала.

Соглашения (владелец 2026-08-18) -- подробности в docstring `service_calc`:
  * затухание -- ОТРИЦАТЕЛЬНЫЕ децибелы по мощности (10*log10 T);
  * SET OFFSET (theta0) -- физический офсет совмещения WGP1/WGP2, зашит в
    калибровку, подгоночный параметр модели; автокалибровка -- следующая версия;
  * SET ZERO -- рабочая точка отсчёта на кривой, от которой оператор считает
    добавочное затухание/усиление; по умолчанию = SET OFFSET;
  * P_0 -- мощность на входе аттенюатора (знаменатель T в режиме absolute).

График: две панели с ЖЁСТКИМИ шкалами -- пропускание 0-100 %, затухание
0...-45 дБ (автоскейл выключен; верх раздвигается только если сдвинутый SET ZERO
уводит кривую в усиление, иначе часть данных была бы просто не видна). Блёклая
серая кривая -- всегда полная мощность по всей полосе, яркая -- выбранная
метрика. Точки запроса отмечены отрезками до осей с подписями угла, пропускания
и затухания.

НЕ клиентский `attenuator_app.gui`/`cli` (v0.2/v0.3, отдельный трек, не
трогается) -- самостоятельный инструмент поверх модели C8
(`measured_curve.py`), обвязка не пересчитывает физику заново, только вызывает
`service_calc`.

Запуск из корня репозитория:
    .venv\\Scripts\\python.exe -m attenuator_app.tools.service_gui
"""
from __future__ import annotations

import sys
import tkinter as tk
from pathlib import Path
from tkinter import ttk, messagebox

import numpy as np
import matplotlib
matplotlib.use("TkAgg")
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg  # noqa: E402
from matplotlib.figure import Figure                              # noqa: E402

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))

from attenuator_app.tools.service_calc import (            # noqa: E402
    FULL, MODE_LABEL, MODES, Metric, REF_POWER_SHORT,
    angle_for_db, attenuation_db, attenuation_db_array, load_calibration,
    power_field_ratios)

#: (kind, подпись в списке, [(подпись поля, значение по умолчанию), ...])
METRIC_ITEMS: list[tuple[str, str, list[tuple[str, str]]]] = [
    ("full", "полная мощность (вся записанная полоса)", []),
    ("single", "одна частота", [("частота, ТГц", "0.800")]),
    ("band_cw", "полоса: центр + ширина", [("центр, ТГц", "0.800"), ("ширина, ТГц", "0.400")]),
    ("band_minmax", "полоса: мин + макс", [("f_min, ТГц", "0.300"), ("f_max, ТГц", "1.500")]),
]
METRIC_BY_LABEL = {label: (kind, fields) for kind, label, fields in METRIC_ITEMS}

PCT_LIM = (0.0, 100.0)      # жёсткая шкала пропускания, %
DB_LIM = (-45.0, 0.0)       # жёсткая шкала затухания, дБ

C_BRIGHT = "#2a78d6"        # выбранная метрика
C_BRIGHT_DB = "#eb6834"
C_DIM = "#c2c2c2"           # полная мощность, фоном
C_MARK = "#12996a"


def _annotate_point(ax, x: float, y: float, x_text: str, y_text: str,
                    color: str = C_MARK, draw_h: bool = True) -> None:
    """Отметить точку (x, y): вертикальный отрезок вниз до оси X и
    горизонтальный влево до оси Y (оба ЗАКАНЧИВАЮТСЯ на графике, а не уходят в
    бесконечность), подписи посередине каждого отрезка. Подписи отодвигаются
    внутрь при приближении к границе поля, значение вне шкалы помечается."""
    xlo, xhi = ax.get_xlim()
    ylo, yhi = ax.get_ylim()
    xspan, yspan = xhi - xlo, yhi - ylo
    off_scale = not (ylo <= y <= yhi)
    y_draw = float(np.clip(y, ylo, yhi))

    ax.plot([x, x], [ylo, y_draw], color=color, lw=1.1, alpha=0.9, zorder=4)
    if draw_h:
        ax.plot([xlo, x], [y_draw, y_draw], color=color, lw=1.1, alpha=0.9, zorder=4)
    if not off_scale:
        ax.plot([x], [y_draw], "o", color=color, ms=4.5, zorder=5)

    box = dict(boxstyle="round,pad=0.18", fc="white", ec=color, lw=0.6, alpha=0.88)
    # подпись угла -- посередине вертикального отрезка, уводим влево у правого края
    near_right = (x - xlo) / xspan > 0.80
    ax.text(x + (-0.012 if near_right else 0.012) * xspan, (ylo + y_draw) / 2.0,
            x_text, ha="right" if near_right else "left", va="center",
            fontsize=8, color=color, bbox=box, zorder=6)
    # подпись значения -- посередине горизонтального отрезка, уводим вниз у верха
    if draw_h:
        near_top = (y_draw - ylo) / yspan > 0.88
        label = y_text + (" (вне шкалы)" if off_scale else "")
        ax.text((xlo + x) / 2.0, y_draw + (-0.02 if near_top else 0.02) * yspan,
                label, ha="center", va="top" if near_top else "bottom",
                fontsize=8, color=color, bbox=box, zorder=6)


class ServiceGUI(tk.Tk):
    def __init__(self):
        super().__init__()
        self.title("Аттенюатор -- обслуживание THz-TDS спектрометра")
        self.geometry("1240x820")
        self.minsize(980, 640)

        self.cal = load_calibration()
        self.offset: float = self.cal.theta0_calibration_deg   # SET OFFSET (физика)
        self.zero: float = self.offset                          # SET ZERO (точка отсчёта)
        self._markers: list[float] = []
        self._field_vars: list[tk.StringVar] = []

        controls = ttk.Frame(self, padding=8)
        controls.pack(fill="x")

        info = ttk.LabelFrame(controls, text="устройство (зашитая калибровка)", padding=6)
        info.pack(fill="x")
        ttk.Label(info, justify="left", text=(
            f"{self.cal.device_id}  --  калибровка {self.cal.dataset} ({self.cal.generated})\n"
            f"P={self.cal.P_um:.2f} мкм, D={self.cal.D_um:.2f} мкм, "
            f"потери={self.cal.loss_db:.3f} дБ/ТГц^{self.cal.gamma:.2f}, "
            f"полоса {self.cal.band_thz[0]:.2f}-{self.cal.band_thz[1]:.2f} ТГц")
        ).pack(anchor="w")

        # --- два разных нуля -------------------------------------------
        zeros = ttk.Frame(controls, padding=(0, 4, 0, 0))
        zeros.pack(fill="x")

        offf = ttk.LabelFrame(zeros, text="SET OFFSET (theta0) -- физический офсет "
                                         "совмещения WGP1/WGP2, входит в модель", padding=6)
        offf.pack(side="left", fill="both", expand=True, padx=(0, 4))
        ttk.Label(offf, justify="left", font=("TkDefaultFont", 8), text=(
            "подобран по калибровочной серии и зашит; менять только если прибор\n"
            "переустановили в роторе (автокалибровка офсета -- следующая версия)")
        ).pack(anchor="w")
        r1 = ttk.Frame(offf)
        r1.pack(fill="x", pady=(4, 0))
        self.offset_var = tk.StringVar(value=f"{self.offset:.3f}")
        ttk.Entry(r1, textvariable=self.offset_var, width=10).pack(side="left")
        ttk.Button(r1, text="применить", command=self._apply_offset).pack(side="left", padx=(6, 0))
        ttk.Button(r1, text="из калибровки", command=self._restore_offset).pack(side="left", padx=(4, 0))
        self.offset_status = tk.StringVar()
        ttk.Label(r1, textvariable=self.offset_status, foreground="#357").pack(side="left", padx=(10, 0))

        zerf = ttk.LabelFrame(zeros, text="SET ZERO -- рабочая точка отсчёта на кривой",
                              padding=6)
        zerf.pack(side="left", fill="both", expand=True, padx=(4, 0))
        ttk.Label(zerf, justify="left", font=("TkDefaultFont", 8), text=(
            "от неё считается ДОБАВОЧНОЕ затухание (или усиление, если идти\n"
            "к совмещению); в режиме relative она же точка нормировки (0 дБ)")
        ).pack(anchor="w")
        r2 = ttk.Frame(zerf)
        r2.pack(fill="x", pady=(4, 0))
        self.zero_var = tk.StringVar(value=f"{self.zero:.3f}")
        ttk.Entry(r2, textvariable=self.zero_var, width=10).pack(side="left")
        ttk.Button(r2, text="SET ZERO", command=self._apply_zero).pack(side="left", padx=(6, 0))
        ttk.Button(r2, text="= SET OFFSET", command=self._reset_zero).pack(side="left", padx=(4, 0))
        self.zero_status = tk.StringVar()
        ttk.Label(r2, textvariable=self.zero_status, foreground="#357").pack(side="left", padx=(10, 0))

        # --- режим и метрика -------------------------------------------
        opts = ttk.Frame(controls, padding=(0, 6, 0, 0))
        opts.pack(fill="x")
        ttk.Label(opts, text="режим:").pack(side="left")
        self.mode_var = tk.StringVar(value="relative")
        for value in MODES:
            ttk.Radiobutton(opts, text=MODE_LABEL[value], value=value, variable=self.mode_var,
                            command=self._redraw_plot).pack(side="left", padx=(4, 10))

        met = ttk.Frame(controls, padding=(0, 4, 0, 0))
        met.pack(fill="x")
        ttk.Label(met, text="метрика:").pack(side="left")
        self.metric_var = tk.StringVar(value=METRIC_ITEMS[0][1])
        cb = ttk.Combobox(met, textvariable=self.metric_var, state="readonly", width=38,
                          values=[label for _, label, _ in METRIC_ITEMS])
        cb.pack(side="left", padx=(6, 10))
        cb.bind("<<ComboboxSelected>>", lambda _e: self._rebuild_metric_fields())
        self.metric_fields = ttk.Frame(met)
        self.metric_fields.pack(side="left")
        ttk.Button(met, text="обновить график", command=self._redraw_plot).pack(side="left", padx=(10, 0))

        # --- калькулятор ------------------------------------------------
        calc = ttk.Frame(controls, padding=(0, 6, 0, 0))
        calc.pack(fill="x")
        fwd = ttk.LabelFrame(calc, text="угол -> затухание", padding=6)
        fwd.pack(side="left", fill="x", expand=True, padx=(0, 4))
        ttk.Label(fwd, text="угол WGP1, град:").pack(side="left")
        self.angle_var = tk.StringVar(value="")
        e = ttk.Entry(fwd, textvariable=self.angle_var, width=10)
        e.pack(side="left", padx=(6, 10))
        e.bind("<Return>", lambda _e: self._forward())
        ttk.Button(fwd, text="Вычислить затухание", command=self._forward).pack(side="left")

        inv = ttk.LabelFrame(calc, text="затухание -> угол (дБ отрицательные)", padding=6)
        inv.pack(side="left", fill="x", expand=True, padx=(4, 0))
        ttk.Label(inv, text="желаемое затухание, дБ:").pack(side="left")
        self.db_var = tk.StringVar(value="")
        e2 = ttk.Entry(inv, textvariable=self.db_var, width=10)
        e2.pack(side="left", padx=(6, 10))
        e2.bind("<Return>", lambda _e: self._inverse())
        ttk.Button(inv, text="Вычислить угол", command=self._inverse).pack(side="left")

        body = ttk.Frame(self, padding=(8, 0, 8, 8))
        body.pack(fill="both", expand=True)
        self.plot_frame = ttk.Frame(body)
        self.plot_frame.pack(side="left", fill="both", expand=True)
        self.canvas = None

        outf = ttk.LabelFrame(body, text="результат", padding=6)
        outf.pack(side="right", fill="y", padx=(8, 0))
        self.text = tk.Text(outf, wrap="word", font=("Consolas", 9), width=48)
        self.text.pack(fill="both", expand=True)

        self._sync_status()
        self._rebuild_metric_fields()

    # -- вспомогательное -------------------------------------------------
    def _log(self, msg: str = "") -> None:
        self.text.insert("end", msg + "\n")
        self.text.see("end")

    def _mode(self) -> str:
        return self.mode_var.get()

    def _sync_status(self) -> None:
        same_off = np.isclose(self.offset, self.cal.theta0_calibration_deg, atol=1e-6)
        self.offset_status.set(f"{self.offset:+.3f} град "
                               f"({'из калибровки' if same_off else 'вручную'})")
        same_zero = np.isclose(self.zero, self.offset, atol=1e-6)
        self.zero_status.set(f"{self.zero:+.3f} град "
                             f"({'= SET OFFSET' if same_zero else 'сдвинут'})")

    def _rebuild_metric_fields(self) -> None:
        for w in self.metric_fields.winfo_children():
            w.destroy()
        kind, fields = METRIC_BY_LABEL[self.metric_var.get()]
        self._field_vars = []
        for label, default in fields:
            ttk.Label(self.metric_fields, text=label + ":").pack(side="left", padx=(0, 4))
            var = tk.StringVar(value=default)
            ent = ttk.Entry(self.metric_fields, textvariable=var, width=9)
            ent.pack(side="left", padx=(0, 10))
            ent.bind("<Return>", lambda _e: self._redraw_plot())
            self._field_vars.append(var)
        self._redraw_plot()

    def _metric(self) -> Metric:
        """Собрать метрику из выпадающего списка и полей. ValueError при мусоре."""
        kind, fields = METRIC_BY_LABEL[self.metric_var.get()]
        vals: list[float] = []
        for var, (label, _) in zip(self._field_vars, fields):
            s = var.get().strip()
            if not s:
                raise ValueError(f"не заполнено поле «{label}»")
            try:
                vals.append(float(s))
            except ValueError:
                raise ValueError(f"поле «{label}»: ожидается число, получено {s!r}") from None
        return Metric(kind, *(vals + [None, None])[:2])

    # -- график ------------------------------------------------------
    def _redraw_plot(self) -> None:
        mode = self._mode()
        try:
            metric = self._metric()
            metric_err = None
        except ValueError as e:
            metric, metric_err = FULL, str(e)

        grid = np.linspace(self.offset - 90.0, self.offset + 90.0, 721)
        try:
            db_sel = attenuation_db_array(grid, self.offset, self.cal, metric, mode, self.zero)
        except ValueError as e:                       # напр. пустая полоса
            metric, metric_err = FULL, str(e)
            db_sel = attenuation_db_array(grid, self.offset, self.cal, FULL, mode, self.zero)
        db_full = (db_sel if metric.kind == "full" else
                   attenuation_db_array(grid, self.offset, self.cal, FULL, mode, self.zero))
        pct_sel, pct_full = 10.0 ** (db_sel / 10.0) * 100.0, 10.0 ** (db_full / 10.0) * 100.0

        fig = Figure(figsize=(6.6, 6.4), dpi=100)
        ax1 = fig.add_subplot(211)
        ax2 = fig.add_subplot(212, sharex=ax1)

        for ax, y_full, y_sel, color in ((ax1, pct_full, pct_sel, C_BRIGHT),
                                         (ax2, db_full, db_sel, C_BRIGHT_DB)):
            if metric.kind != "full":
                ax.plot(grid, y_full, color=C_DIM, lw=1.4, zorder=1,
                        label="полная мощность")
            ax.plot(grid, y_sel, color=color, lw=1.8, zorder=2, label=metric.label)
            ax.axvline(self.offset, color="#9aa", ls="--", lw=1, zorder=0)
            if not np.isclose(self.zero, self.offset, atol=1e-6):
                ax.axvline(self.zero, color=C_MARK, ls=":", lw=1.2, zorder=0)
            ax.grid(alpha=0.25, lw=0.6)

        # ЖЁСТКИЕ шкалы; верх раздвигаем только при усилении (сдвинутый SET ZERO)
        ax1.set_xlim(grid[0], grid[-1])
        ax1.set_ylim(PCT_LIM[0], max(PCT_LIM[1], float(np.nanmax(pct_sel)) * 1.02))
        ax2.set_ylim(DB_LIM[0], max(DB_LIM[1], float(np.nanmax(db_sel)) + 0.5))

        ref = REF_POWER_SHORT[mode]
        ax1.set_ylabel(f"пропускание T = P/{ref}, %")
        ax1.set_title(f"{MODE_LABEL[mode]}\nSET OFFSET {self.offset:+.2f}°, "
                      f"SET ZERO {self.zero:+.2f}°", fontsize=8.5)
        ax1.legend(fontsize=7.5, loc="upper left")
        ax2.set_ylabel("затухание, дБ по мощности")
        ax2.set_xlabel("угол WGP1 (показание шкалы ротатора), град")

        # симметричные решения дают одинаковое значение -- горизонтальный
        # отрезок и его подпись рисуем один раз, иначе подписи наложатся
        drawn: list[float] = []
        for theta in self._markers:
            db = attenuation_db(theta, self.offset, self.cal, metric, mode, self.zero)
            pct = 10.0 ** (db / 10.0) * 100.0
            draw_h = not any(abs(db - d) < 1e-3 for d in drawn)
            drawn.append(db)
            _annotate_point(ax1, theta, pct, f"{theta:+.2f}°", f"T = {pct:.2f} %", draw_h=draw_h)
            _annotate_point(ax2, theta, db, f"{theta:+.2f}°", f"{db:+.2f} дБ", draw_h=draw_h)

        fig.tight_layout()
        if self.canvas is not None:
            self.canvas.get_tk_widget().destroy()
        self.canvas = FigureCanvasTkAgg(fig, master=self.plot_frame)
        self.canvas.draw()
        self.canvas.get_tk_widget().pack(fill="both", expand=True)

        if metric_err:
            self._log(f"[!] метрика: {metric_err} -- график по полной мощности")

    # -- действия ----------------------------------------------------
    def _apply_offset(self) -> None:
        try:
            v = float(self.offset_var.get())
        except ValueError:
            messagebox.showerror("Неверное значение", "SET OFFSET должен быть числом")
            return
        # SET ZERO -- точка на шкале ротатора: если он был привязан к офсету,
        # едет вместе с ним, иначе остаётся там, где его поставил оператор
        was_bound = bool(np.isclose(self.zero, self.offset, atol=1e-6))
        self.offset = v
        if was_bound:
            self.zero = v
            self.zero_var.set(f"{v:.3f}")
        self._markers.clear()
        self._sync_status()
        self._log(f"SET OFFSET = {self.offset:+.3f} град")
        self._log("")
        self._redraw_plot()

    def _restore_offset(self) -> None:
        self.offset_var.set(f"{self.cal.theta0_calibration_deg:.3f}")
        self._apply_offset()

    def _apply_zero(self) -> None:
        try:
            v = float(self.zero_var.get())
        except ValueError:
            messagebox.showerror("Неверное значение", "SET ZERO должен быть числом")
            return
        self.zero = v
        self._sync_status()
        db_zero_abs = attenuation_db(self.zero, self.offset, self.cal, FULL, "absolute")
        self._log(f"SET ZERO = {self.zero:+.3f} град "
                  f"(там абсолютное пропускание {10.0 ** (db_zero_abs / 10.0) * 100:.2f} % от P_0)")
        self._log("")
        self._redraw_plot()

    def _reset_zero(self) -> None:
        self.zero_var.set(f"{self.offset:.3f}")
        self._apply_zero()

    def _forward(self) -> None:
        try:
            theta = float(self.angle_var.get())
        except ValueError:
            messagebox.showerror("Неверное значение", "угол должен быть числом")
            return
        try:
            metric = self._metric()
        except ValueError as e:
            messagebox.showerror("Метрика", str(e))
            return
        mode = self._mode()
        w = metric.warning(self.cal)
        if w:
            self._log(w)

        db = attenuation_db(theta, self.offset, self.cal, metric, mode, self.zero)
        p_r, f_r = power_field_ratios(db)
        db_zero = attenuation_db(self.zero, self.offset, self.cal, metric, mode, self.zero)
        d_zero = db - db_zero

        self._log(f"[{mode}] {metric.label}")
        self._log(f"  угол WGP1 = {theta:+.3f} град "
                  f"(от OFFSET {theta - self.offset:+.3f}, от ZERO {theta - self.zero:+.3f})")
        self._log(f"  затухание = {db:+.2f} дБ по мощности")
        self._log(f"  T = P/{REF_POWER_SHORT[mode]} = {p_r * 100:.3f} %   "
                  f"(поле E/E_ref = {f_r:.4g})")
        self._log(f"  относительно SET ZERO: {d_zero:+.2f} дБ "
                  f"({'усиление' if d_zero > 0 else 'затухание'})")
        self._log("")

        self._markers = [theta]
        self._redraw_plot()

    def _inverse(self) -> None:
        try:
            target = float(self.db_var.get())
        except ValueError:
            messagebox.showerror("Неверное значение", "затухание должно быть числом")
            return
        try:
            metric = self._metric()
        except ValueError as e:
            messagebox.showerror("Метрика", str(e))
            return
        mode = self._mode()
        w = metric.warning(self.cal)
        if w:
            self._log(w)
        try:
            sol = angle_for_db(target, self.offset, self.cal, metric, mode, self.zero)
        except ValueError as e:
            messagebox.showerror("Недостижимо", str(e))
            return

        self._log(f"[{mode}] {metric.label}")
        self._log(f"  желаемое затухание {target:+.2f} дБ по мощности; диапазон "
                  f"[{sol['db_min']:.2f}, {sol['db_max']:.2f}] дБ")
        self._log(f"  угол WGP1 = {sol['theta_plus_deg']:+.3f} град")
        self._log(f"  или       = {sol['theta_minus_deg']:+.3f} град")
        self._log("  (два симметричных решения -- выбрать ближе к текущему положению ротатора)")
        self._log("")

        self._markers = [sol["theta_plus_deg"], sol["theta_minus_deg"]]
        self._redraw_plot()


def main() -> int:
    app = ServiceGUI()
    app.mainloop()
    return 0


if __name__ == "__main__":
    sys.exit(main())
