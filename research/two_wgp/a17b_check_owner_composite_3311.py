# -*- coding: utf-8 -*-
"""Численный разбор композита владельца по test_grid_33_11."""
import os
import sys

import numpy as np
from PIL import Image

os.chdir(r"C:\Users\pop\Antigravity IDE\THz-Unified-Optimizer")
sys.path.insert(0, "research/experiments")
import micrograph_period as mp

SRC = (r"C:\Users\pop\.claude\uploads\06cadc2f-5f09-4afb-9788-b28982710213"
       r"\966f8d74-17860503912167399274243632846908.png")
im = Image.open(SRC).convert("L")
A = np.asarray(im, float)
H, W = A.shape
print(f"композит: {W} x {H}")

# граница между гребёнкой (сверху) и образцом (снизу) — по резкому изменению
# среднего периода строк; ищем по яркости: гребёнка светлее в среднем
rowmean = A.mean(axis=1)
split = None
for y in range(int(0.15 * H), int(0.5 * H)):
    if rowmean[y - 3:y].mean() - rowmean[y:y + 3].mean() > 8:
        split = y
        break
if split is None:
    split = int(0.27 * H)
print(f"граница гребёнка/образец: y = {split}")


def per_of(y0, y1, pmin):
    prof = A[y0:y1, :].mean(axis=0)
    per, _s, _p = mp.fft_fundamental(prof, pmin_px=pmin)
    peaks, spac = mp.peak_spacings(prof, per)
    return prof, per, peaks, (float(np.median(spac)) if len(spac) else float("nan"))


prof_c, per_c, pk_c, med_c = per_of(6, max(split - 6, 20), 4.0)
prof_w, per_w, pk_w, med_w = per_of(split + 20, H - 6, 8.0)
print(f"  гребёнка: период FFT {per_c:.2f} px (по пикам {med_c:.2f}), штрихов {len(pk_c)}")
print(f"  образец : период FFT {per_w:.2f} px (по пикам {med_w:.2f}), проволок {len(pk_w)}")

scale = 10.0 / per_c
print(f"\nмасштаб композита: {scale:.4f} мкм/px (1 деление = 10 мкм = {per_c:.2f} px)")
print(f"ПЕРИОД РЕШЁТКИ:    {per_w:.2f} px -> **{per_w * scale:.2f} мкм**")

w10 = mp.peak_widths_at(prof_w, pk_w, per_w, frac=0.10)
w50 = mp.peak_widths_at(prof_w, pk_w, per_w, frac=0.50)
good = w10[(w10 > 0.2 * np.median(w10)) & (w10 < 2.5 * np.median(w10))]
print(f"\nширины проволок в композите ({len(good)} шт.):")
print("  на 10 %: " + ", ".join(f"{x * scale:.2f}" for x in np.sort(good)) + " мкм")
print(f"  медиана 10 %: {np.median(good) * scale:.2f} мкм | "
      f"медиана FWHM: {np.median(w50) * scale:.2f} мкм")
print(f"  D/P по 10 %: {np.median(good) / per_w:.3f}")

owner = np.array([11.86, 11.37, 10.70, 11.42, 10.93])
print(f"\nваши ручные: {', '.join(f'{x:.2f}' for x in np.sort(owner))} мкм; "
      f"медиана {np.median(owner):.2f}, среднее {owner.mean():.2f}")
print(f"вы приняли: 11.2 мкм")
print("\n--- сверка по test_grid_33_11 ---")
print(f"  период:  композит {per_w * scale:6.2f} | моё (среднее увел.) 31.82 | "
      f"моё (большое, общая калибровка) 31.44")
print(f"  диаметр: композит {np.median(good) * scale:6.2f} | ваш ручной {np.median(owner):.2f} | "
      f"моё A16 12.60 | принято 11.2")
