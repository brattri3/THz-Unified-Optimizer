# PureWave angular run — roadmap (TIME-OPTIMIZED for one working day)

One trace = **6 min 15 s** (130 ps window). This plan: **54 signal** + ~9–19 background traces.

- Est. time @6:15/trace: **~6.6 h** (fits a day).
- If you SHORTEN the scan window to ~65 ps (drops only the far water-line tail we mask anyway; keeps 0.2–1.5 THz + echoes): ~3:10/trace -> **~3.3 h** (then raise crossed reps to 6–8).

## Zones & reps

- **bright (coarse, 1 rep)**: [0, 20, 40, 55, 65, 72]
- **crossed (2 deg, 4 reps)**: [80, 82, 84, 86, 88, 90, 92, 94, 96, 98, 100]
- **beyond crossed (2 reps)**: [105, 110]

## Background policy

Measure background (empty beam) about **every 2 angles** or **every ~30 min**, and copy that file into each `_bg` slot it covers. Backgrounds drift slowly — no need to pair one per signal. (Files provided: one `*_rep1_bg.txt` per angle.)

## Suggested order (drift-robust)

Do the whole set of angles for rep1 first (ascending), then rep2 (descending), etc. — so each angle's reps are spread across the session and laser drift decorrelates from angle. Crossed angles (80–100) get 4 reps; bright/beyond fewer.

## Per-angle steps

1) Rotate WGP. 2) Signal -> `purewave_<angle>deg_rep<N>_sig.txt`. 3) Background every ~2 angles -> `purewave_<angle>deg_rep1_bg.txt` (copy to neighbours). 4) Note the deep-shadow (near-crossed) SNR stays >10x; add reps there if weak. 5) True crossed minimum is found later from the data (Specac was ~84°, not 90°).
