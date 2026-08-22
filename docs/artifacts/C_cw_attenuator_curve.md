# C — CW angular attenuation curve for `attenuator_app`

Task report (canonical copy). Owner-approved 2026-08-22. Role C, Zone C.

## Goal

A CW (continuous-wave) variant of `attenuator_app`: a **non-interactive** CLI
command that, for a monochromatic CW source, renders a PNG of the angular
dependence of attenuation A(θ). Base = the **client core** (`attenuator_app/core/*`),
not the service track. Reuse existing physics; write no new physics. Output PNG
only (matplotlib Agg); no CSV; slope panel off by default.

## What was built

* `attenuator_app/cw_curve.py` — new non-interactive `argparse` CLI + a numeric
  self-check (`--selftest`). Pure assembly over `core/forward.py`,
  `core/passport.py`, `core/plots.py`. The interactive `cli.py` is untouched.
* `docs/attenuator_app/09_CW_CURVE.md` — usage note (toggles, angle convention).
* `docs/attenuator_app/10_CALIBRATION_DESIGN.md` — DESIGN-ONLY calibration
  algorithm (angle set, DR determination, deep-extinction handling). No code.
* `research/results/cw_attenuator/*.png` — six demonstration PNGs.
* `coordination/QUESTIONS.md` — rows C-7..C-10 (all non-blocking, defaults taken).
* `attenuator_app/STATE.md` — updated.

## Interface

```
python -m attenuator_app.cw_curve --freq 0.2 --detector coherent --analyzer H
python -m attenuator_app.cw_curve --detector power --source depolarized
python -m attenuator_app.cw_curve --selftest
```

Two independent toggles: **detector** (`coherent` = polarization-sensitive, has
analyzer axis / `power` = insensitive, total power) and **source** (`linear`
DOP=1 / `depolarized` DOP=0). Analyzer presets H=0°, V=90°, or arbitrary degrees.
Sweep −90..+90°, step 1° (all parametrized). Default passport
`ATT-11-16-CA85_02721.json` (same as `cli.py`).

## Angle convention (two rotators)

Chosen default, straight from the passport `calibration_note` + code: first WGP
fixed at the source axis (θ₁=0), second (analyzer) rotates, attenuation set by the
**relative angle θ=θ₂−θ₁**; θ is swept. Calibrated scheme S0 ≡ two-WGP S1 by the
alignment theorem (`selftest` check 6). Recorded for confirmation as C-7.

## Source toggle implementation

Reused the existing DOP mechanism (`Setup.dop` + unpolarized branch of
`forward.intensity`). `linear→dop=1.0`, `depolarized→dop=0.0`. The mechanism was
**sufficient**; no new physics. Only the two extreme DOP values are exposed (C-8).

## Verification

* `python -m attenuator_app.selftest` → **13/13** (no regression).
* `python -m attenuator_app.cw_curve --selftest` → **5/5**:
  1. power+linear = Malus cos² (=0.5×ideal cos⁴), max|dA| 0.12 dB over |θ|≤60°.
  2. coherent+H+linear = ideal cos⁴, max|dA| 0.15 dB over |θ|≤45°, A(0)=0.
  3. coherent+V+linear = deep extinction at θ=0: A(0)=42.1 dB vs A(45)=6.0 dB.
  4. depolarized+power = flat (rotational invariance), max|A|=3.9e-15 dB.
  5. CW point == broadband spectrum at same freq, max|dA|=0 dB.

## Demonstration PNGs (`research/results/cw_attenuator/`)

| file | detector | source | analyzer | shows |
|---|---|---|---|---|
| `cw_sensitive_linear_H_200GHz.png` | coherent | linear | H | classic cos⁴ extinction, floor ~53 dB at ±90° |
| `cw_sensitive_linear_V_200GHz.png` | coherent | linear | V (abs) | crossed analyzer: extinction at θ=0, peak near ±45° |
| `cw_insensitive_linear_200GHz.png` | power | linear | — | Malus cos² law |
| `cw_insensitive_depolarized_200GHz.png` | power | depolarized | — | flat 0 dB (rotational invariance) |
| `cw_sensitive_depolarized_H_200GHz.png` | coherent | depolarized | H | cos² (depolarization removes one Malus factor) |
| `cw_sensitive_linear_H_200GHz_slope.png` | coherent | linear | H | same as #1 with the optional `|dA/dθ|` slope panel |

All at CW 0.2 THz on the default passport.

## Calibration design note (summary)

`10_CALIBRATION_DESIGN.md`, DESIGN-ONLY. Standard non-uniform angle set (coarse
0..±60°, fine ±70..±90° where slope is largest, ≥3 repeats for between-session σ);
θ₀ from the transmission maximum; A_meas = −10log10(P/P₀). **Dynamic range**
measured in the same session as `DR=10·log10(P₀/P_dark)` (P_dark = blocked-beam
reading, the same quantity `limits.dynamic_range_db` consumes). Deep-extinction
points above `DR−MARGIN` (or within K·σ of the dark floor) are reported as **lower
bounds**, not measurements, and compared against the model's own finite floor
(~40 dB, `auto_cross`). DR-extension options (averaging, two-range source-power
stitching, coherent leakage subtraction) noted, not auto-applied.

## Open questions

C-7 (angle convention), C-8 (depolarized source model / partial DOP), C-9 (H/V
presets and −90..+90 range), C-10 (calibration DR thresholds). All **non-blocking**,
defaults documented in-line in QUESTIONS.md and in the code/docs.
