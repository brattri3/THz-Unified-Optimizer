# -*- coding: utf-8 -*-
"""Численная проверка композита владельца: гребёнка сверху, проволоки снизу."""
import os
import sys

import numpy as np
from PIL import Image

os.chdir(r"C:\Users\pop\Antigravity IDE\THz-Unified-Optimizer")
sys.path.insert(0, "research/experiments")
import micrograph_period as mp

SRC = r"C:\Users\pop\.claude\uploads\06cadc2f-5f09-4afb-9788-b28982710213\13babce2-1000000516.jpg"
im = Image.open(SRC).convert("L")
A = np.asarray(im, float)
print("скриншот:", im.size)

# Области в координатах ОРИГИНАЛА (displayed x1.40).
# верх — гребёнка, низ — проволоки; берём полосы подальше от подписей.
BANDS = {
    "гребёнка (верх)": (585, 1010, 210, 1630),      # y0,y1, x0,x1
    "проволоки (низ)": (1080, 1560, 210, 1630),
}


def period_of(y0, y1, x0, x1, pmin=6.0):
    sub = A[y0:y1, x0:x1]
    prof = sub.mean(axis=0)
    per, _s, _p = mp.fft_fundamental(prof, pmin_px=pmin)
    peaks, spac = mp.peak_spacings(prof, per)
    return prof, per, peaks, spac


res = {}
for tag, (y0, y1, x0, x1) in BANDS.items():
    prof, per, peaks, spac = period_of(y0, y1, x0, x1)
    med = float(np.median(spac)) if len(spac) else float("nan")
    res[tag] = (prof, per, peaks)
    print(f"  {tag:<20} период FFT {per:7.2f} px, по пикам {med:7.2f} px, "
          f"найдено {len(peaks)} штрихов")

per_comb = res["гребёнка (верх)"][1]
per_wire = res["проволоки (низ)"][1]
scale = 10.0 / per_comb                      # мкм/px в системе КОМПОЗИТА
print(f"\nмасштаб композита по гребёнке: {scale:.4f} мкм/px "
      f"(1 деление = 10 мкм = {per_comb:.2f} px)")
print(f"период решётки в композите:   {per_wire:.2f} px -> **{per_wire * scale:.2f} мкм**")

# ширина проволок на 10 % в композите
prof_w, per_w, peaks_w = res["проволоки (низ)"]
w10 = mp.peak_widths_at(prof_w, peaks_w, per_w, frac=0.10)
w50 = mp.peak_widths_at(prof_w, peaks_w, per_w, frac=0.50)
print(f"\nширина проволок в композите ({len(w10)} шт.):")
print(f"  на 10 % от пика: " + ", ".join(f"{x * scale:.2f}" for x in np.sort(w10)) + " мкм")
print(f"  на 50 % (FWHM):  " + ", ".join(f"{x * scale:.2f}" for x in np.sort(w50)) + " мкм")
print(f"  медиана 10 %: {np.median(w10) * scale:.2f} мкм | "
      f"медиана FWHM: {np.median(w50) * scale:.2f} мкм")
print(f"  D/P по 10 %: {np.median(w10) / per_w:.3f}")

print("\n--- сверка с моими числами по purewave ---")
print(f"  период:  композит {per_wire * scale:6.2f} | моё (среднее увел.) 24.72 | "
      f"моё (общая калибровка большого) 26.08")
print(f"  диаметр: композит {np.median(w10) * scale:6.2f} | моё по всем кадрам 9.83 | "
      f"владелец вручную ~10.0")
