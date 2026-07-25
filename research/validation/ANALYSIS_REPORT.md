# Micrograph comprehensive analysis (auto-generated)

Calibration µm/px: high=0.1270, mid=0.6443, low≈1.289 (approx — low fine comb unresolved; low absolute lengths indicative, % metrics valid).


## PureWave

### Period (µm) — 4 methods

| mag | N | FFT | autocorr | NN-median | lattice | spread |
|-----|---|-----|----------|-----------|---------|--------|
| low | 98 | 24.49 | 48.97 | 24.50 | 24.48 | ±10.60 |
| mid | 51 | 24.81 | 24.48 | 24.53 | 24.71 | ±0.13 |
| high | 9 | 26.33 | 26.91 | 25.97 | 25.75 | ±0.44 |

### Diameter (µm) — high mag, in-focus wires; 4 methods

| frame | FWHM(glint) | width@20% | width@10% | duty(Otsu) |
|-------|-------------|-----------|-----------|------------|
| highmag_01 | 6.1 | 9.1 | 10.7 | 7.3 |
| highmag_02 | 5.8 | 9.0 | 10.9 | 7.3 |
| highmag_03 | 5.4 | 8.7 | 10.1 | 6.9 |

### Periodicity disorder — 3 methods

| mag | N | NN-spacing CV% | lattice-RMS% | missing/merged % |
|-----|---|----------------|--------------|------------------|
| low | 98 | 16.8 | 17.7 | 9.3 |
| mid | 51 | 16.0 | 14.0 | 3.8 |
| high | 9 | 8.5 | 12.9 | 6.7 |

## Grid_D20_P40

### Period (µm) — 4 methods

| mag | N | FFT | autocorr | NN-median | lattice | spread |
|-----|---|-----|----------|-----------|---------|--------|
| low | 60 | 38.56 | 38.66 | 38.74 | 38.54 | ±0.08 |
| mid | 31 | 38.72 | 38.66 | 39.17 | 38.86 | ±0.20 |
| high | 7 | 40.39 | 38.97 | 38.39 | 38.77 | ±0.75 |

### Diameter (µm) — high mag, in-focus wires; 4 methods

| frame | FWHM(glint) | width@20% | width@10% | duty(Otsu) |
|-------|-------------|-----------|-----------|------------|
| highmag_focus0um | 9.1 | 14.5 | 17.3 | 13.7 |
| highmag_focus200um | 7.6 | 11.9 | 14.1 | 12.2 |

### Periodicity disorder — 3 methods

| mag | N | NN-spacing CV% | lattice-RMS% | missing/merged % |
|-----|---|----------------|--------------|------------------|
| low | 60 | 15.4 | 16.5 | 13.0 |
| mid | 31 | 16.1 | 11.2 | 8.8 |
| high | 7 | 10.2 | 9.2 | 0.0 |

---

## Consolidated conclusions

**Period — robust (4 methods agree).** FFT, NN-median and lattice-fit agree to
<1 µm at mid/high; autocorrelation confirms (it occasionally locks to 2×period at
low mag — the PureWave-low 48.97 µm entry — a known ACF ambiguity, ignore it).
- **PureWave P = 25.5 ± 0.9 µm** (cross-magnification spread = systematic floor).
- **Grid D20/P40 P = 38.8 ± 0.3 µm** (nominal 40).

**Diameter — method-dependent, round-wire optics.** The ladder FWHM < width@20% <
width@10% ≈ duty is the specular-glint → full-footprint progression. Best proxy =
base width @10% on in-focus wires at max mag:
- **PureWave D ≈ 10.6 µm** (fill ≈ 0.42). FWHM-glint (5.7) under-reads ~2×.
- **Grid D ≈ 17 µm** (fill ≈ 0.45; sharpest frame focus0um), nominal 20; shortfall
  = specular fall-off on curved wire flanks.
Take-away: micrographs pin **P** hard, **D** as a prior with a known low bias.

**Periodicity disorder — 3 methods, consistent, upper bound.** NN-CV%, lattice-RMS%
and missing-wire% all land at ~10–17 %, falling at high mag where detection is
cleanest. This is an **upper bound**: inflated by peak-detection/index error and a
few genuine missing/weak wires; the sharp FFT fundamental rules out true 15 %
continuous jitter. Body-of-distribution positional scatter ≈ 3–5 µm (~8–13 % of
pitch) — wound-grid regularity, worse than lithographic.

**Grid-specific:** wires span ~200 µm along the beam (two-focus test) — a distinct,
THz-relevant defect (predicts 0.67–1.33 ps delayed 90° leakage; see MANIFEST).

### Reliability by magnification
| quantity | low mag | mid mag | high mag |
|----------|---------|---------|----------|
| Period (statistics, N wires) | **best** (N~60–100) | good | few (7–10) |
| Period (absolute calib) | weak (comb unresolved) | **good** | **best** |
| Diameter | no (blur) | fair | **best** (in-focus) |
| Disorder % (calib-free) | **best** (large N) | good | few |
