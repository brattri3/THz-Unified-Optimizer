"""tkinter GUI поверх `service_calc.py` -- обслуживание аттенюатора в THz-TDS
спектрометре (задача C9, санкция владельца 2026-08-19). Один сценарий,
минимальный: SET ZERO вручную + двунаправленный калькулятор (дБ<->угол),
интегрально и/или на выбранной частоте, без доверительного интервала.

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

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))

from attenuator_app.tools.service_calc import (            # noqa: E402
    angle_for_db, attenuation_db, band_warning, load_calibration)


class ServiceGUI(tk.Tk):
    def __init__(self):
        super().__init__()
        self.title("Аттенюатор -- обслуживание THz-TDS спектрометра")
        self.geometry("640x580")
        self.minsize(560, 480)

        self.cal = load_calibration()
        self.theta0: float | None = None

        info = ttk.LabelFrame(self, text="устройство (зашитая калибровка)", padding=6)
        info.pack(fill="x", padx=8, pady=(8, 4))
        ttk.Label(info, justify="left", text=(
            f"{self.cal.device_id}  --  калибровка {self.cal.dataset} ({self.cal.generated})\n"
            f"P={self.cal.P_um:.2f} мкм, D={self.cal.D_um:.2f} мкм, "
            f"потери={self.cal.loss_db:.3f} дБ/ТГц^{self.cal.gamma:.2f}, "
            f"полоса {self.cal.band_thz[0]:.2f}-{self.cal.band_thz[1]:.2f} ТГц")
        ).pack(anchor="w")

        zero = ttk.LabelFrame(
            self, text="SET ZERO -- офсет ротатора (алгоритм авто-калибровки -- в следующей версии)",
            padding=6)
        zero.pack(fill="x", padx=8, pady=4)
        ttk.Label(zero, justify="left", text=(
            "показание шкалы WGP1 в СОВМЕЩЁННОМ положении\n"
            "(минимум затухания, WGP1 || WGP2/детектор), град:")
        ).pack(anchor="w")
        row = ttk.Frame(zero)
        row.pack(fill="x", pady=(4, 0))
        self.zero_var = tk.StringVar(value="0.000")
        ttk.Entry(row, textvariable=self.zero_var, width=12).pack(side="left")
        ttk.Button(row, text="SET ZERO", command=self._set_zero).pack(side="left", padx=(8, 0))
        self.zero_status = tk.StringVar(value="ноль НЕ установлен")
        ttk.Label(row, textvariable=self.zero_status, foreground="#a33").pack(side="left", padx=(12, 0))

        freqf = ttk.Frame(self, padding=(8, 4, 8, 0))
        freqf.pack(fill="x")
        ttk.Label(freqf, text="частота, ТГц (необязательно -- иначе интегрально по всей полосе):").pack(side="left")
        self.freq_var = tk.StringVar(value="")
        ttk.Entry(freqf, textvariable=self.freq_var, width=10).pack(side="left", padx=(6, 0))

        fwd = ttk.LabelFrame(self, text="угол -> затухание", padding=6)
        fwd.pack(fill="x", padx=8, pady=4)
        ttk.Label(fwd, text="угол WGP1, град:").pack(side="left")
        self.angle_var = tk.StringVar(value="")
        ttk.Entry(fwd, textvariable=self.angle_var, width=10).pack(side="left", padx=(6, 10))
        ttk.Button(fwd, text="Вычислить затухание", command=self._forward).pack(side="left")

        inv = ttk.LabelFrame(self, text="затухание -> угол", padding=6)
        inv.pack(fill="x", padx=8, pady=4)
        ttk.Label(inv, text="желаемое затухание, дБ:").pack(side="left")
        self.db_var = tk.StringVar(value="")
        ttk.Entry(inv, textvariable=self.db_var, width=10).pack(side="left", padx=(6, 10))
        ttk.Button(inv, text="Вычислить угол", command=self._inverse).pack(side="left")

        outf = ttk.LabelFrame(self, text="результат", padding=6)
        outf.pack(fill="both", expand=True, padx=8, pady=(4, 8))
        self.text = tk.Text(outf, wrap="word", font=("Consolas", 10), height=14)
        self.text.pack(fill="both", expand=True)

    # -- вспомогательное -------------------------------------------------
    def _log(self, msg: str = "") -> None:
        self.text.insert("end", msg + "\n")
        self.text.see("end")

    def _parse_freq(self) -> float | None:
        s = self.freq_var.get().strip()
        return float(s) if s else None

    # -- действия ----------------------------------------------------
    def _set_zero(self) -> None:
        try:
            z = float(self.zero_var.get())
        except ValueError:
            messagebox.showerror("Неверное значение", "показание шкалы должно быть числом")
            return
        self.theta0 = z
        self.zero_status.set(f"ноль установлен: theta0 = {z:+.3f} град")
        self.text.delete("1.0", "end")
        self._log(f"SET ZERO: theta0 = {z:+.3f} град")
        self._log("")

    def _forward(self) -> None:
        if self.theta0 is None:
            messagebox.showinfo("SET ZERO", "сначала выполните SET ZERO")
            return
        try:
            theta = float(self.angle_var.get())
            freq = self._parse_freq()
        except ValueError:
            messagebox.showerror("Неверное значение", "угол и частота должны быть числами")
            return

        db_int = attenuation_db(theta, self.theta0, self.cal, None)
        self._log(f"угол WGP1 = {theta:+.3f} град (delta от нуля = {theta - self.theta0:+.3f})")
        self._log(f"  интегральное затухание = {db_int:.2f} дБ")
        if freq is not None:
            w = band_warning(freq, self.cal)
            if w:
                self._log(f"  {w}")
            db_f = attenuation_db(theta, self.theta0, self.cal, freq)
            self._log(f"  затухание на {freq:.3f} ТГц = {db_f:.2f} дБ")
        self._log("")

    def _inverse(self) -> None:
        if self.theta0 is None:
            messagebox.showinfo("SET ZERO", "сначала выполните SET ZERO")
            return
        try:
            target = float(self.db_var.get())
            freq = self._parse_freq()
        except ValueError:
            messagebox.showerror("Неверное значение", "затухание и частота должны быть числами")
            return

        w = band_warning(freq, self.cal)
        if w:
            self._log(w)
        try:
            sol = angle_for_db(target, self.theta0, self.cal, freq)
        except ValueError as e:
            messagebox.showerror("Недостижимо", str(e))
            return

        metric = f"на {freq:.3f} ТГц" if freq is not None else "интегрально"
        self._log(f"желаемое затухание {target:.2f} дБ ({metric}); "
                 f"максимум на этой калибровке {sol['floor_db']:.2f} дБ")
        self._log(f"  угол WGP1 = {sol['theta_plus_deg']:+.3f} град")
        self._log(f"  или       = {sol['theta_minus_deg']:+.3f} град")
        self._log("  (два симметричных решения -- выбрать ближе к текущему положению ротатора)")
        self._log("")


def main() -> int:
    app = ServiceGUI()
    app.mainloop()
    return 0


if __name__ == "__main__":
    sys.exit(main())
