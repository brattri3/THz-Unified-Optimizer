"""GUI-обёртка над `measured_curve.py` — прототип для демонстрации в отчёте.

Выбор серии измерений (`att-11-16-s2`/`s3`), кнопка «Посчитать», график и
сводка прямо в окне. Расчёт полностью в `measured_curve.analyze()` — здесь
только tkinter-обвязка, ничего не пересчитывается заново.

Ручная подгонка Бланко (P, D, потери, offset ротатора) — поля ввода + кнопка
«Пересчитать» (не слайдеры: владелец 2026-08-17 попросил именно поля, точное
значение важнее плавности перетаскивания).

ВНИМАНИЕ: прототип для завтрашнего отчёта, без глубокого тестирования и без
сверки с зоной A (санкция владельца 2026-08-17). Читает `data_pool/`, которого
нет в клиентском Win7-комплекте -- это лабораторный инструмент, не часть
`attenuator_app.gui`/`cli`, которые уезжают покупателю.

Запуск из корня репозитория:
    .venv\\Scripts\\python.exe -m attenuator_app.tools.measured_gui
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

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))

from attenuator_app.tools.measured_curve import (          # noqa: E402
    DATA_ROOT, GAMMA_FIXED, analyze, analyze_road, blanco_channel_model,
    build_figure, channels_at_theta0, summary_text)

#: номинал паспорта (attenuator_app/tools/make_passport.py:SPEC) -- только для
#: подсказки в заголовке блока ручной подгонки, поля стартуют с результата
#: авто-фита, а не с номинала.
P_NOMINAL_UM = 16.0
D_NOMINAL_UM = 11.0
LOSS_NOMINAL_DB = 0.255

#: только серии, где вращается ПЕРВЫЙ WGP (`data_pool/two_wgp_attenuator/README.md`).
#: `att-11-16-356`/`s1` вращают ВТОРОЙ WGP -- другая схема, подпись "угол
#: первого WGP" на графике для них была бы неверной.
FLAT_DATASETS = ["att-11-16-s2", "att-11-16-s3"]
ROAD_PREFIX = "дорога: "


def find_road_folders():
    """Подкаталоги data_pool/two_wgp_attenuator с дорожной конвенцией (NNN_aXX.txt)."""
    if not DATA_ROOT.exists():
        return []
    return sorted(p.name for p in DATA_ROOT.iterdir()
                  if p.is_dir() and any(p.glob("[0-9][0-9][0-9]_a*.txt")))


class MeasuredCurveGUI(tk.Tk):
    def __init__(self):
        super().__init__()
        self.title("Аттенюатор -- угловая кривая по измерениям (прототип)")
        self.geometry("1200x820")
        self.minsize(900, 620)

        values = list(FLAT_DATASETS) + [ROAD_PREFIX + n for n in find_road_folders()]

        top = ttk.Frame(self, padding=6)
        top.pack(fill="x")
        ttk.Label(top, text="Серия:").pack(side="left")
        self.dataset_var = tk.StringVar(value=values[0] if values else "")
        cb = ttk.Combobox(top, textvariable=self.dataset_var, values=values,
                          state="readonly", width=32)
        cb.pack(side="left", padx=(4, 10))
        ttk.Label(top, text="полоса, ТГц:").pack(side="left")
        self.band_lo = tk.StringVar(value="0.3")
        self.band_hi = tk.StringVar(value="1.5")
        ttk.Entry(top, textvariable=self.band_lo, width=6).pack(side="left", padx=(3, 2))
        ttk.Entry(top, textvariable=self.band_hi, width=6).pack(side="left", padx=(2, 10))
        ttk.Button(top, text="Посчитать", command=self._run).pack(side="left")
        self.status_var = tk.StringVar(value="выберите серию и нажмите «Посчитать»")
        ttk.Label(top, textvariable=self.status_var, foreground="#666",
                  wraplength=760, justify="left").pack(side="left", padx=(12, 0))

        # -- ручная подгонка: поля ввода + кнопка «Пересчитать» --------------
        manual = ttk.LabelFrame(
            self, text=f"Ручная подгонка Бланко (номинал паспорта: "
                      f"P={P_NOMINAL_UM:g} D={D_NOMINAL_UM:g} мкм, "
                      f"потери={LOSS_NOMINAL_DB:g} дБ/ТГц^{GAMMA_FIXED:g})", padding=6)
        manual.pack(fill="x", padx=6, pady=(0, 6))

        self.p_var = tk.StringVar(value=f"{P_NOMINAL_UM:.3f}")
        self.d_var = tk.StringVar(value=f"{D_NOMINAL_UM:.3f}")
        self.loss_var = tk.StringVar(value=f"{LOSS_NOMINAL_DB:.4f}")
        self.gamma_var = tk.StringVar(value=f"{GAMMA_FIXED:.3f}")
        self.theta0_var = tk.StringVar(value="0.000")

        for label, var, width in (("P, мкм:", self.p_var, 10),
                                  ("D, мкм:", self.d_var, 10),
                                  ("потери, дБ/ТГц^gamma:", self.loss_var, 10),
                                  ("gamma:", self.gamma_var, 7),
                                  ("офсет ротатора theta0, град:", self.theta0_var, 10)):
            ttk.Label(manual, text=label).pack(side="left", padx=(0 if var is self.p_var else 10, 3))
            ttk.Entry(manual, textvariable=var, width=width).pack(side="left")
        ttk.Button(manual, text="Пересчитать", command=self._on_recalc).pack(
            side="left", padx=(14, 0))

        self.absolute_var = tk.BooleanVar(value=False)
        ttk.Checkbutton(
            manual, text="абсолютное пропускание (Бланко-кривая; данные всегда отн.)",
            variable=self.absolute_var, command=self._on_recalc).pack(side="left", padx=(14, 0))

        body = ttk.Frame(self, padding=(6, 0, 6, 6))
        body.pack(fill="both", expand=True)

        self.plot_frame = ttk.Frame(body)
        self.plot_frame.pack(side="left", fill="both", expand=True)
        self.canvas = None
        self.r = None    # последний результат analyze()/analyze_road()

        self.text = tk.Text(body, width=48, wrap="word", font=("Consolas", 9))
        self.text.pack(side="right", fill="y")

    # -- ручной пересчёт ---------------------------------------------------
    def _on_recalc(self):
        if self.r is None:
            messagebox.showinfo("Пересчитать", "Сначала нажмите «Посчитать» -- нужны данные серии.")
            return
        try:
            p_um = float(self.p_var.get())
            d_um = float(self.d_var.get())
            loss = float(self.loss_var.get())
            gamma = float(self.gamma_var.get())
            theta0 = float(self.theta0_var.get())
        except ValueError:
            messagebox.showerror("Неверное значение",
                                 "P, D, потери, gamma и theta0 должны быть числами.")
            return
        if p_um <= 0 or d_um <= 0:
            messagebox.showerror("Неверное значение", "P и D должны быть положительными.")
            return

        # theta0 меняет РАЗЛОЖЕНИЕ на каналы (данные), P/D/потери/gamma --
        # модельную кривую поверх них; порядок обязателен именно такой.
        self.r["I_perp"], self.r["I_par"] = channels_at_theta0(self.r["spec_fit"], theta0)
        self.r["P_fit"], self.r["D_fit"] = p_um, d_um
        self.r["loss_fit"], self.r["gamma_fit"] = loss, gamma
        self.r["theta0_manual"] = theta0
        self.r["angular_absolute"] = self.absolute_var.get()

        model_perp, model_par = blanco_channel_model(self.r["freqs_band"], p_um, d_um, loss,
                                                      gamma=gamma)
        rms_perp = float(np.sqrt(np.mean((model_perp - self.r["I_perp"]) ** 2)))
        rms_par = float(np.sqrt(np.mean((model_par - self.r["I_par"]) ** 2)))
        self.status_var.set(
            f"ручной P={p_um:.3f} D={d_um:.3f} потери={loss:.4f} gamma={gamma:.3f} "
            f"theta0={theta0:+.3f} -- RMS(t_perp)={rms_perp:.4f}, RMS(t_par)={rms_par:.2e}")
        self._redraw()

    # -- загрузка серии ------------------------------------------------
    def _run(self):
        dataset = self.dataset_var.get()
        try:
            band = (float(self.band_lo.get()), float(self.band_hi.get()))
        except ValueError:
            messagebox.showerror("Неверная полоса", "полоса должна быть числами, ТГц")
            return
        self.status_var.set(f"считаю {dataset}...")
        self.update_idletasks()
        try:
            if dataset.startswith(ROAD_PREFIX):
                r = analyze_road(DATA_ROOT / dataset[len(ROAD_PREFIX):], band)
            else:
                r = analyze(dataset, DATA_ROOT, band)
        except ValueError as e:
            messagebox.showerror("Нет данных", str(e))
            self.status_var.set("ошибка -- см. диалог")
            return
        except Exception as e:                                     # noqa: BLE001
            messagebox.showerror("Расчёт не удался", f"{type(e).__name__}: {e}")
            self.status_var.set("ошибка -- см. диалог")
            return

        self.r = r
        self.text.delete("1.0", "end")
        self.text.insert("1.0", summary_text(r))

        # поля стартуют с результата авто-фита (P/D/потери) и с офсета,
        # который сам нашёл угловой фит Малюса (theta0_fit) -- дальше можно
        # трогать руками и жать «Пересчитать»
        self.p_var.set(f"{r['P_fit']:.3f}")
        self.d_var.set(f"{r['D_fit']:.3f}")
        self.loss_var.set(f"{r['loss_fit']:.4f}")
        self.gamma_var.set(f"{r.get('gamma_fit', GAMMA_FIXED):.3f}")
        self.theta0_var.set(f"{r['theta0_fit']:.3f}")
        self._on_recalc()
        self.status_var.set(f"готово: {dataset}")

    def _redraw(self):
        fig = build_figure(self.r)
        if self.canvas is not None:
            self.canvas.get_tk_widget().destroy()
        self.canvas = FigureCanvasTkAgg(fig, master=self.plot_frame)
        self.canvas.draw()
        self.canvas.get_tk_widget().pack(fill="both", expand=True)


def main():
    app = MeasuredCurveGUI()
    app.mainloop()


if __name__ == "__main__":
    sys.exit(main() or 0)
