# External validation — micrograph ground truth & measurement index

Geometric ground truth (true period **P** and wire diameter **D**) measured from
microscope images, to validate the best model ideas on independent samples —
above all a direct check of the effective-geometry law `D_eff` (T5b) and
cross-sample transfer of the leakage cluster (HN2/HN7/HN8) and non-Drude γ.

## Naming convention (micrographs)
`research/validation/micrographs/<sampleID>/`
```
<sampleID>__cal__<low|mid|high>mag__div0p01mm.bmp   # graticule, 1 div = 0.01 mm = 10 um
<sampleID>__wgp__<low|mid|high>mag[_NN].bmp          # the polarizer, same mag as its cal
<sampleID>__ANALYSIS_overlay.png                     # produced by micrograph_period.py
```
A cal and a wgp image at the **same** `mag` token share the µm/px calibration.

## Method
`research/experiments/micrograph_period.py`. Calibration: FFT fundamental of the
graticule fine comb = 10 µm → µm/px. WGP: FFT fundamental → period; FWHM of the
bright wire stripes → D; ideal-lattice fit of wire centres → positional disorder.
Period is confirmed by cross-magnification agreement (built-in validation).

**Diameter method.** Reflection microscopy on round wires shows a narrow specular
glint over a wider bell; the wire diameter is the **full reflecting footprint**,
so D is the base width at ~10 % of peak height (FWHM/50 % captures only the glint
and under-reads by ~2×). Always measure D at **maximum magnification on in-focus
wires only** — for non-coplanar samples pass `in_focus_only=True` so defocused
(broadened) wires are excluded.

## Calibration achieved (PureWave optics)
| mag  | fine-comb spacing | µm/px  | note |
|------|-------------------|--------|------|
| high | 78.77 px = 10 µm  | 0.1270 | clean, harmonics at 39/26/20 px |
| mid  | 15.52 px = 10 µm  | 0.6443 | 85 teeth, clean |
| low  | ~7.6 px = 10 µm   | —      | fine comb below optical resolution → not self-calibrating; use for ratio only |

## Ground-truth table
| sampleID | nominal D/P | **P meas ±σ (µm)** | **D meas ±σ (µm)** | fill D/P | disorder (RMS) | measurements | angles |
|----------|-------------|--------------------|--------------------|----------|----------------|--------------|--------|
| purewave | (per model) | **25.5 ± 0.9**     | **10.6 ± 0.5**     | 0.42     | ≤~13% + sparse missing wires | `data_pool/purewave/` (TODO) | TODO |
| specac   | (per model) | **24.9 ± 0.9** (midmag FFT+lattice) | **~14 ± 3** (glint-underread ~8–12 × ~1.9 corr.) | ~0.55 | ~18% upper / 5 missing | **`data_pool/specac/` READY** (13 ang 0–100°) | 0,10,…,80,84,90,96,100 |

**Specac (2026-07-25):** THz data organized to canonical names `specac_<angle>deg_rep1_{sig,bg}.txt` (13 angles; drift-aware bg groups bg1–bg7, one bg per ~2 angles; fine sampling 84/90/96 in the shadow). Model-independent leakage floor **η=0.0357, deep-shadow SNR 68×** (best of all 6 samples). Micrograph P reliable (24.9); D is a wide prior (glint under-reads; run `analyze_wgp(..., diam_frac=0.10, in_focus_only=True)` on highmag for a firmer D). No Specac graticule shot — PureWave calibration transferred (same optics), validated by P≈25.
| specac   | (per model) | TODO               | TODO               | TODO     | TODO           | `data_pool/specac/` (TODO)   | TODO |
| 356att   | D11 / P16   | TODO (add micros)  | TODO               | TODO     | TODO           | `data_pool/356att…`          | ✓ |
| test_grid_40_20 | D20 / P40 | **38.8 ± 0.3** | **~17 (in-focus, nom 20)** | 0.45 | in-plane RMS ~3–5 µm (~8–13%); **+ wires span ~200 µm in z** | `data_pool/test_grid…` | ✓ |

### PureWave notes (measured 2026-07-24)
- **P = 25.5 ± 0.9 µm** — dominant lattice period; σ is the cross-magnification
  calibration spread (mid 24.7 µm vs high 25.8 µm), the systematic floor. FFT
  peak is sharp with clean harmonics ⇒ well-defined mean lattice.
- **D = 10.6 ± 0.5 µm** (fill factor D/P ≈ 0.42) — full reflecting-footprint
  width of the wires at high mag (base width at 10 % of peak; see "Diameter
  method" below). The earlier 5.7 µm figure was the specular-glint FWHM (half
  height) and under-read D by ~2×; superseded.
- **Disorder** — nearest-neighbour scatter ~10–14% at the detected-peak level is
  an **upper bound**: dominated by peak-detection/index error plus a few genuine
  missing/weak wires (visible ~165–225 µm in `highmag_01`; ~2 missing per 50 in
  midmag). Not continuous jitter — the sharp FFT rules that out.
- See `micrographs/purewave/purewave__ANALYSIS_overlay.png`.

### test_grid_40_20 (D20/P40) notes (2026-07-24)
- **No calibration slide captured** for this sample. PureWave's µm/px was
  transferred (same mag tokens) and yields P = 39.2 / 38.7 / 40.4 µm at
  low/mid/high — all ≈ nominal 40 µm, which both confirms the pitch and
  validates re-using the PureWave calibration. Still, capture a graticule at
  each mag for this sample to remove the assumption.
- **P = 38.8 ± 0.3 µm** (nominal 40) — cross-mag consistent (low 38.5 / mid 38.9 /
  high 38.7); PureWave calibration transferred (optics unchanged), and the fact
  that P lands on 40 µm validates the transfer.
- **D ≈ 17 µm (fill factor ≈ 0.45; nominal 20).** Measured at max magnification
  on the **in-focus wires only** (`in_focus_only=True`), base width at 10 % of
  peak: sharpest frame `focus0um` gives 17.3 ± 0.9 µm, `focus200um` 14.1 µm (its
  sharp wires are slightly softer). The small shortfall vs 20 µm is the round-
  wire specular fall-off (curved flanks reflect below threshold). Still: use
  micro-P as hard ground truth and D as a prior.
- **Winding quality (in-plane):** wire centres scatter RMS ≈ 3–5 µm (~8–13 % of
  pitch) about the ideal lattice — the body is within ±3–4 µm with a few larger
  local defects/missed wires (see `…__ANALYSIS_winding.png`). This is
  wound-wire-grid regularity, markedly less perfect than a lithographic grid.
- **Wires are NOT coplanar.** Two high-mag frames focused at f = 0 µm and
  f = 200 µm each bring a *different* subset of wires into focus ⇒ wires span
  ≈ 200 µm along the beam axis z (two layers, or a bowed/wavy single grid over
  the ~260 µm field). Reflection microscopy, bright = wire.
- **THz consequence (Δz = 200 µm):** extra optical path per plane → delay
  τ = Δz/c = **0.67 ps** (single pass) or 2Δz/c = **1.33 ps** (inter-plane
  etalon); spectral ripple 1/τ ≈ **1.5 THz** (or 0.75 THz etalon). Predicts a
  **time-delayed residual at 90°** — a signature that separates *sample* leakage
  (delayed) from *detector* cross-polarization (zero delay) → direct test for
  the HN2/HN7/HN8 leakage cluster (T12). Also a candidate contributor to the
  angle-dependent D_eff (T5b): under any residual beam/grid tilt θ the two
  depth planes shift laterally by Δz·tanθ (≈17 µm at 5° ≈ half the 40 µm pitch),
  modulating the effective fill factor with angle. See files
  `micrographs/wgp_grid_D20_P40/…focus0um.bmp` / `…focus200um.bmp`.
