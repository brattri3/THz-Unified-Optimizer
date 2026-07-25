"""Prepare Monday's PureWave angular run: empty measurement files + order roadmap.

TIME-OPTIMIZED design (one trace = 6 min 15 s at the 130 ps window).
A full uniform 5-rep sweep (210 traces) = ~22 h (3 days) -> infeasible.
This concentrates the measurement budget where the physics is (near-crossed) and
uses periodic backgrounds, fitting a single working day.

Zones (per-zone reps):
  - bright 0..72 deg : 1 rep  (high SNR, few points constrain the Malus curve)
  - crossed 80..100  : 4 reps (2-deg; low SNR + all the D_eff<->eta<->gamma
                               degeneracy lives here -> need per-point sigma)
  - beyond 105,110   : 2 reps (extra t_par leverage)
Backgrounds: ONE per angle (measured periodically / every ~2 angles and reused
across that angle's reps by copying) -> ~19 bg instead of 1-per-trace.
Total ~54 sig + ~19 bg placeholders; ~64 actual traces if bg every ~2 angles
=> ~6.7 h at 6:15/trace, or ~3.4 h if the scan window is shortened to ~65 ps.
"""
from pathlib import Path

DS = "purewave"
POOL = Path(__file__).resolve().parents[2] / "data_pool" / DS
POOL.mkdir(parents=True, exist_ok=True)

ZONES = [  # (angles, reps, label)
    ([0, 20, 40, 55, 65, 72], 1, "bright (coarse, 1 rep)"),
    ([80, 82, 84, 86, 88, 90, 92, 94, 96, 98, 100], 4, "crossed (2 deg, 4 reps)"),
    ([105, 110], 2, "beyond crossed (2 reps)"),
]
TRACE_MIN = 6.25   # min per trace @130 ps window


def main():
    # clear stale EMPTY placeholders we created before (never touch non-empty = real data)
    removed = 0
    for f in POOL.glob(f"{DS}_*deg_rep*_*.txt"):
        if f.stat().st_size == 0:
            f.unlink(); removed += 1

    made = 0; n_sig = 0; n_bg = 0
    order_rows = []
    for angles, reps, label in ZONES:
        for ang in angles:
            for rep in range(1, reps + 1):
                s = POOL / f"{DS}_{ang}deg_rep{rep}_sig.txt"
                if not s.exists():
                    s.write_text("", encoding="utf-8"); made += 1; n_sig += 1
            # one background per angle (reuse across reps; measure every ~2 angles)
            b = POOL / f"{DS}_{ang}deg_rep1_bg.txt"
            if not b.exists():
                b.write_text("", encoding="utf-8"); made += 1; n_bg += 1
            order_rows.append((ang, reps, label))

    n_traces = n_sig + max(1, n_bg // 2)   # bg every ~2 angles in practice
    hrs_full = n_sig * TRACE_MIN / 60 + (n_bg // 2) * TRACE_MIN / 60
    hrs_short = hrs_full * (65.0 / 130.0)

    L = []
    L.append("# PureWave angular run — roadmap (TIME-OPTIMIZED for one working day)\n")
    L.append(f"One trace = **6 min 15 s** (130 ps window). This plan: **{n_sig} signal** + "
             f"~{n_bg//2}–{n_bg} background traces.\n")
    L.append(f"- Est. time @6:15/trace: **~{hrs_full:.1f} h** (fits a day).")
    L.append(f"- If you SHORTEN the scan window to ~65 ps (drops only the far water-line tail we mask "
             f"anyway; keeps 0.2–1.5 THz + echoes): ~3:10/trace -> **~{hrs_short:.1f} h** "
             f"(then raise crossed reps to 6–8).\n")
    L.append("## Zones & reps\n")
    for angles, reps, label in ZONES:
        L.append(f"- **{label}**: {angles}")
    L.append("")
    L.append("## Background policy\n")
    L.append("Measure background (empty beam) about **every 2 angles** or **every ~30 min**, and copy "
             "that file into each `_bg` slot it covers. Backgrounds drift slowly — no need to pair one "
             "per signal. (Files provided: one `*_rep1_bg.txt` per angle.)\n")
    L.append("## Suggested order (drift-robust)\n")
    L.append("Do the whole set of angles for rep1 first (ascending), then rep2 (descending), etc. — so "
             "each angle's reps are spread across the session and laser drift decorrelates from angle. "
             "Crossed angles (80–100) get 4 reps; bright/beyond fewer.\n")
    L.append("## Per-angle steps\n")
    L.append("1) Rotate WGP. 2) Signal -> `purewave_<angle>deg_rep<N>_sig.txt`. 3) Background every "
             "~2 angles -> `purewave_<angle>deg_rep1_bg.txt` (copy to neighbours). 4) Note the deep-"
             "shadow (near-crossed) SNR stays >10x; add reps there if weak. 5) True crossed minimum is "
             "found later from the data (Specac was ~84°, not 90°).\n")
    (POOL / "MEASUREMENT_ORDER.md").write_text("\n".join(L), encoding="utf-8")

    print(f"removed {removed} old empty placeholders")
    print(f"created {made} files: {n_sig} sig + {n_bg} bg (1 bg/angle)")
    print(f"~{n_traces} actual traces -> ~{hrs_full:.1f} h @6:15  |  ~{hrs_short:.1f} h @~3:10 (65 ps window)")
    print(f"roadmap -> {POOL/'MEASUREMENT_ORDER.md'}")


if __name__ == "__main__":
    main()
