"""tkinter GUI поверх `service_calc.py` -- обслуживание аттенюатора в THz-TDS
спектрометре (задача C9, санкция владельца 2026-08-19). Один сценарий:
двунаправленный калькулятор (дБ<->угол), интегрально и/или на выбранной
частоте, относительный/абсолютный режим, без доверительного интервала.

theta0 (SET ZERO) по умолчанию берётся из калибровки прибора (подгонка по
серии `ATT-11-16-CA85_02721_s2` уже достаточно хороша -- владелец 2026-08-19)
и приложение готово к расчётам сразу при старте; поле и кнопка SET ZERO
остаются, чтобы оператор мог сдвинуть точку отсчёта вручную (например, если
прибор переустановлен в роторе, или чтобы задать рабочую точку и считать
затухание относительно неё) -- алгоритм автоматической калибровки офсета
запланирован на следующую версию.

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
    angle_for_db, attenuation_db, attenuation_db_array, band_warning,
    load_calibration, power_field_ratios)

MODE_LABELS = [("relative", "относительный (к theta0, T(theta0)=1)"),
              ("absolute", "абсолютный (доля входной мощности)")]


class ServiceGUI(tk.Tk):
    def __init__(self):
        super().__init__()
        self.title("Аттенюатор -- обслуживание THz-TDS спектрометра")
        self.geometry("1160x720")
        self.minsize(920, 600)

        self.cal = load_calibration()
        self.theta0: float = self.cal.theta0_calibration_deg
        self._markers: list[float] = []

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

        zero = ttk.LabelFrame(
            controls, text="theta0 (SET ZERO) -- офсет совмещения WGP1/WGP2; "
                          "алгоритм авто-калибровки -- в следующей версии", padding=6)
        zero.pack(fill="x", pady=(4, 0))
        ttk.Label(zero, justify="left", text=(
            "показание шкалы WGP1, град (по умолчанию -- значение из калибровки; "
            "переустановите, если прибор передвигали в роторе, "
            "или чтобы считать затухание от произвольной рабочей точки):")
        ).pack(anchor="w")
        row = ttk.Frame(zero)
        row.pack(fill="x", pady=(4, 0))
        self.zero_var = tk.StringVar(value=f"{self.theta0:.3f}")
        ttk.Entry(row, textvariable=self.zero_var, width=12).pack(side="left")
        ttk.Button(row, text="SET ZERO", command=self._set_zero).pack(side="left", padx=(8, 0))
        ttk.Button(row, text="из калибровки", command=self._restore_zero).pack(side="left", padx=(6, 0))
        self.zero_status = tk.StringVar(
            value=f"theta0 = {self.theta0:+.3f} град (из калибровки)")
        ttk.Label(row, textvariable=self.zero_status, foreground="#357").pack(side="left", padx=(12, 0))

        opts = ttk.Frame(controls, padding=(0, 6, 0, 0))
        opts.pack(fill="x")
        ttk.Label(opts, text="режим:").pack(side="left")
        self.mode_var = tk.StringVar(value="relative")
        for value, label in MODE_LABELS:
            ttk.Radiobutton(opts, text=label, value=value, variable=self.mode_var,
                           command=self._redraw_plot).pack(side="left", padx=(4, 10))
        ttk.Label(opts, text="частота, ТГц (пусто = интегрально):").pack(side="left", padx=(16, 0))
        self.freq_var = tk.StringVar(value="")
        ttk.Entry(opts, textvariable=self.freq_var, width=8).pack(side="left", padx=(6, 0))
        ttk.Button(opts, text="обновить график", command=self._redraw_plot).pack(side="left", padx=(8, 0))

        calc = ttk.Frame(controls, padding=(0, 6, 0, 0))
        calc.pack(fill="x")
        fwd = ttk.LabelFrame(calc, text="угол -> затухание", padding=6)
        fwd.pack(side="left", fill="x", expand=True, padx=(0, 4))
        ttk.Label(fwd, text="угол WGP1, град:").pack(side="left")
        self.angle_var = tk.StringVar(value="")
        ttk.Entry(fwd, textvariable=self.angle_var, width=10).pack(side="left", padx=(6, 10))
        ttk.Button(fwd, text="Вычислить затухание", command=self._forward).pack(side="left")

        inv = ttk.LabelFrame(calc, text="затухание -> угол", padding=6)
        inv.pack(side="left", fill="x", expand=True, padx=(4, 0))
        ttk.Label(inv, text="желаемое затухание, дБ:").pack(side="left")
        self.db_var = tk.StringVar(value="")
        ttk.Entry(inv, textvariable=self.db_var, width=10).pack(side="left", padx=(6, 10))
        ttk.Button(inv, text="Вычислить угол", command=self._inverse).pack(side="left")

        body = ttk.Frame(self, padding=(8, 0, 8, 8))
        body.pack(fill="both", expand=True)

        self.plot_frame = ttk.Frame(body)
        self.plot_frame.pack(side="left", fill="both", expand=True)
        self.canvas = None

        outf = ttk.LabelFrame(body, text="результат", padding=6)
        outf.pack(side="right", fill="y", padx=(8, 0))
        self.text = tk.Text(outf, wrap="word", font=("Consolas", 9), width=46)
        self.text.pack(fill="both", expand=True)

        self._redraw_plot()

    # -- вспомогательное -------------------------------------------------
    def _log(self, msg: str = "") -> None:
        self.text.insert("end", msg + "\n")
        self.text.see("end")

    def _parse_freq(self) -> float | None:
        s = self.freq_var.get().strip()
        return float(s) if s else None

    def _mode(self) -> str:
        return self.mode_var.get()

    # -- график ------------------------------------------------------
    def _redraw_plot(self) -> None:
        try:
            freq = self._parse_freq()
        except ValueError:
            freq = None
        mode = self._mode()
        theta_grid = np.linspace(self.theta0 - 90.0, self.theta0 + 90.0, 361)
        atten_db = attenuation_db_array(theta_grid, self.theta0, self.cal, freq, mode)
        T_pct = 10.0 ** (-atten_db / 10.0) * 100.0

        fig = Figure(figsize=(6.4, 6.0), dpi=100)
        ax1 = fig.add_subplot(211)
        ax2 = fig.add_subplot(212, sharex=ax1)

        for ax in (ax1, ax2):
            ax.axvline(self.theta0, color="#999", ls="--", lw=1,
                      label="theta0 (SET ZERO)" if ax is ax1 else None)
            for m in self._markers:
                ax.axvline(m, color="#1baf7a", lw=1.3, alpha=0.8)

        ax1.plot(theta_grid, T_pct, color="#2a78d6")
        ylabel_pct = ("T(theta)/T(theta0), %" if mode == "relative"
                     else "T = P/P_вход, %")
        ax1.set_ylabel(ylabel_pct)
        ax1.set_title(f"угловая кривая, {mode}, "
                     f"{f'{freq:.3f} ТГц' if freq is not None else 'интегрально'}", fontsize=9)
        ax1.legend(fontsize=7.5, loc="upper right")

        ax2.plot(theta_grid, atten_db, color="#eb6834")
        ax2.set_ylabel("затухание, дБ (по мощности)")
        ax2.set_xlabel("угол WGP1 (показание шкалы ротатора), град")

        fig.tight_layout()

        if self.canvas is not None:
            self.canvas.get_tk_widget().destroy()
        self.canvas = FigureCanvasTkAgg(fig, master=self.plot_frame)
        self.canvas.draw()
        self.canvas.get_tk_widget().pack(fill="both", expand=True)

    # -- действия ----------------------------------------------------
    def _set_zero(self) -> None:
        try:
            z = float(self.zero_var.get())
        except ValueError:
            messagebox.showerror("Неверное значение", "показание шкалы должно быть числом")
            return
        self.theta0 = z
        manual = not np.isclose(z, self.cal.theta0_calibration_deg, atol=1e-6)
        self.zero_status.set(f"theta0 = {z:+.3f} град ({'задан вручную' if manual else 'из калибровки'})")
        self._markers.clear()
        self.text.delete("1.0", "end")
        self._log(f"SET ZERO: theta0 = {z:+.3f} град")
        self._log("")
        self._redraw_plot()

    def _restore_zero(self) -> None:
        self.zero_var.set(f"{self.cal.theta0_calibration_deg:.3f}")
        self._set_zero()

    def _forward(self) -> None:
        try:
            theta = float(self.angle_var.get())
            freq = self._parse_freq()
        except ValueError:
            messagebox.showerror("Неверное значение", "угол и частота должны быть числами")
            return
        mode = self._mode()

        db_int = attenuation_db(theta, self.theta0, self.cal, None, mode)
        p_r, f_r = power_field_ratios(db_int)
        self._log(f"[{mode}] угол WGP1 = {theta:+.3f} град (delta от theta0 = {theta - self.theta0:+.3f})")
        self._log(f"  интегральное затухание = {db_int:.2f} дБ по мощности "
                 f"(P/P_ref={p_r:.4g}, поле E/E_ref={f_r:.4g})")
        if freq is not None:
            w = band_warning(freq, self.cal)
            if w:
                self._log(f"  {w}")
            db_f = attenuation_db(theta, self.theta0, self.cal, freq, mode)
            p_rf, f_rf = power_field_ratios(db_f)
            self._log(f"  затухание на {freq:.3f} ТГц = {db_f:.2f} дБ по мощности "
                     f"(P/P_ref={p_rf:.4g}, поле E/E_ref={f_rf:.4g})")
        self._log("")

        self._markers = [theta]
        self._redraw_plot()

    def _inverse(self) -> None:
        try:
            target = float(self.db_var.get())
            freq = self._parse_freq()
        except ValueError:
            messagebox.showerror("Неверное значение", "затухание и частота должны быть числами")
            return
        mode = self._mode()

        w = band_warning(freq, self.cal)
        if w:
            self._log(w)
        try:
            sol = angle_for_db(target, self.theta0, self.cal, freq, mode)
        except ValueError as e:
            messagebox.showerror("Недостижимо", str(e))
            return

        metric = f"на {freq:.3f} ТГц" if freq is not None else "интегрально"
        self._log(f"[{mode}] желаемое затухание {target:.2f} дБ по мощности ({metric}); "
                 f"диапазон на этой калибровке [{sol['floor_min_db']:.2f}, {sol['floor_max_db']:.2f}] дБ")
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
