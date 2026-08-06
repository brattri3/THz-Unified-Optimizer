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
| test_grid_33_11 | D11 / P33 (nominal, owner) | **31.75 ± 0.65** (cross-mag: mid 31.82 / high 31.69, 0.4% apart) | **10.91 ± 0.41** (defocus-corrected, A14b; the raw in-focus-filtered value was 12.55 — withdrawn) | **0.344** (nominal 0.333) | ~22% (median over frames; > purewave ≤13%) | `data_pool/test_grid_33_11/` ✓ 37 angles −100…+100 | ✓ |

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

### test_grid_33_11 notes (measured 2026-08-06 from shots taken 2026-08-04)
- Source: `I:\att\WGP-11-33\` (owner's USB), 9 frames at three magnifications
  plus a focus pair 0 / 0.12 mm. Filed under the naming convention above.
- ⚠ **No graticule was shot in that session** — the scale is TRANSFERRED from 2026-07-24
  (the same purewave graticule frames), as was done for specac. Transfer validated twice,
  independently: (a) pixel-period ratio between zoom steps high/mid = 5.019 vs July's scale
  ratio 5.098 (1.5% apart), mid/low = 2.027 vs the expected 2; (b) after scaling, the period
  in microns agrees across magnifications to 0.4%. Both checks must be redone if the
  microscope zoom is ever reset.
- **The measured period beat the nominal one on the THz data too**: at equal parameter count
  the M5 multistart fit scores 7.1 AIC better with P = 31.75 than with P = 33
  (`research/results/two_wgp/a14_refit_measured_geometry.json`). Two independent channels agree.
- Consequence for the D_eff law: D_eff/D_phys = **0.576** (was 0.694 on nominal geometry),
  H5 predicts 0.664 at the measured D/P ⇒ deviation **−0.088**, the largest negative of the
  eight samples. The earlier claim that H5 passed an extrapolative test is withdrawn.
- Wires are **non-coplanar**: the 0 / 0.12 mm focus pair brings different wires into focus
  (4 each). Same caveat as test_grid_40_20 (~200 µm span) ⇒ relevant to HN12.
- Analysis: `research/two_wgp/a14_micrograph_33_11.py`, artefact
  `research/results/two_wgp/a14_micrograph_33_11.json`.

### ⚠ Diameter is defocus-sensitive — correction added 2026-08-06 (A14b)

The owner pointed out that at maximum magnification not all wires are in focus and the image
smears; measuring `D` on such a frame over-reads. Confirmed, and the cause is in the tooling:
`micrograph_period.in_focus_mask` selects by **quantiles within the frame** (top 45 % by edge
steepness, top 60 % by brightness), so it always keeps about half the wires **even when all of
them are blurred**. Fine for a coplanar sample, wrong for one whose wires sit at different depths.

**Correction procedure** (`research/two_wgp/a14b_diameter_defocus.py`). Per wire measure
`w10` (width at 10 % of peak) and `b = peak height / max |gradient|` — the inverse normalised
edge steepness, which grows linearly with the blur width. Regress `w10` on `b` and evaluate at
the **instrument's own blur floor**, measured from the graticule (its edges are physically sharp,
so all blur there is instrumental): **8.69 px = 1.10 µm** at high mag. Do **not** extrapolate to
`b → 0` — that removes the instrument PSF as well and under-reads.

**The slope is instrumental, not per-sample:** +0.211 µm/px on test_grid_33_11 versus +0.195 on
purewave, 8 % apart. What differs between samples is the blur SPREAD (purewave 13–27 px,
test_grid_33_11 10–49 px), i.e. how coplanar the winding is.

| sample | slope, µm/px | D published | **D defocus-corrected** |
|---|---|---|---|
| test_grid_33_11 | +0.211 | (12.55, withdrawn) | **10.91 ± 0.41** |
| purewave | +0.195 | 10.6 ± 0.5 | 9.49 ± 2.08 |

⚠ **Consequence for the whole table.** Every published `D` here carries a defocus contribution
that differs per sample, so the `D/P` column is not strictly comparable across rows. Recomputing
all of them would also retune the empirical `D_eff` law (H5), which was fitted on these very
`D/P` values. Only test_grid_33_11 has been corrected so far — see `research/two_wgp/state.json`,
task `A16_defocus_recompute_all_samples` (proposed, needs the owner's decision).
