"""T5 (step 5b) — extract effective diameter D_eff by matching the Blanco model
to the validated lattice-sum reference, and test the empirical law
    D_eff / D_phys = 1 - 0.85 (D/P).

Method: the perpendicular channel t_perp is ~1 and D-insensitive in the
sub-wavelength band (useless for extraction); the parallel channel t_par is the
D-sensitive (inductive/reflective) one. For each (D/P, freq) we compute the
exact |t_par| with the PEC-cylinder lattice-sum solver (physical D), then find
the effective ratio x=D_eff/P for which Blanco's compute_t_par matches, and
report D_eff/D_phys vs the empirical law.

Extraction is well-conditioned for D/P >= 0.5 (|t_par| monotone increasing);
below that |t_par| is non-monotone (inductive/capacitive crossover) so those
points are flagged low-confidence.

Read-only on originals; output under research/results/latticesum/.
"""
import sys, json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.optimize import minimize_scalar

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import t5_latticesum as t5
sys.path.insert(0, str(HERE.parents[1]))   # repo root
from unified_optimizer import model_blanco as mb

OUT = HERE.parents[0] / "results" / "latticesum"
OUT.mkdir(parents=True, exist_ok=True)

C = 3e8
P = 15.5e-6
FREQS_THZ = [0.25, 0.5, 1.0, 1.5]
DP_GRID = [0.3, 0.4, 0.5, 0.6, 0.7, 0.8]
LAW = lambda dp: 1.0 - 0.85*dp


def exact_tpar(dp, nu_thz):
    k = 2*np.pi*nu_thz*1e12/C
    a = (dp*P)/2.0
    T, _ = t5.grating_TR(k, P, a, "par", L=8)
    return abs(T)


def blanco_tpar_abs(x_dp, pol_lam, nu_thz):
    return abs(mb.compute_t_par(pol_lam, x_dp, freq_thz=nu_thz, use_drude=False))


def extract_deff(dp, nu_thz):
    lam = C/(nu_thz*1e12)
    pol_lam = P/lam
    target = exact_tpar(dp, nu_thz)
    # search x=D_eff/P minimizing |blanco - target|
    res = minimize_scalar(lambda x: (blanco_tpar_abs(x, pol_lam, nu_thz)-target)**2,
                          bounds=(0.03, 0.90), method="bounded",
                          options={"xatol": 1e-4})
    x = float(res.x); miss = float(np.sqrt(res.fun))
    return {"D_over_P": dp, "nu_THz": nu_thz, "P_over_lambda": pol_lam,
            "exact_tpar": float(target), "blanco_tpar_atfit": float(blanco_tpar_abs(x, pol_lam, nu_thz)),
            "Deff_over_P": x, "Deff_over_Dphys": float(x/dp), "residual": miss,
            "law_value": float(LAW(dp)), "well_conditioned": bool(dp >= 0.5)}


def main():
    rows = []
    for nu in FREQS_THZ:
        for dp in DP_GRID:
            rows.append(extract_deff(dp, nu))

    # figure: D_eff/D_phys vs D/P per frequency + empirical law
    fig, ax = plt.subplots(figsize=(7.5, 5))
    dp_fine = np.linspace(0.25, 0.85, 100)
    ax.plot(dp_fine, LAW(dp_fine), "k--", lw=2, label="empirical law 1−0.85(D/P)")
    for nu in FREQS_THZ:
        sub = [r for r in rows if r["nu_THz"] == nu]
        dps = [r["D_over_P"] for r in sub]
        rat = [r["Deff_over_Dphys"] for r in sub]
        ax.plot(dps, rat, "o-", ms=5, label=f"exact-matched, {nu} THz")
    ax.axhline(1.0, color="gray", lw=0.6)
    ax.set_xlabel("D/P (physical density)"); ax.set_ylabel("D_eff / D_phys")
    ax.set_title("D_eff from lattice-sum matching vs empirical law")
    ax.legend(fontsize=8); ax.grid(alpha=0.3); fig.tight_layout()
    fig.savefig(OUT / "deff_vs_law.png", dpi=120); plt.close(fig)

    # summary numbers (well-conditioned D/P>=0.5, averaged over freq)
    wc = [r for r in rows if r["well_conditioned"]]
    by_dp = {}
    for dp in sorted({r["D_over_P"] for r in wc}):
        vals = [r["Deff_over_Dphys"] for r in wc if r["D_over_P"] == dp]
        by_dp[dp] = {"mean_Deff_over_Dphys": float(np.mean(vals)),
                     "std": float(np.std(vals)), "law": float(LAW(dp))}

    out = {"geometry": {"P_um": P*1e6}, "rows": rows, "by_dp_wellcond": by_dp,
           "conclusion_metric": {
               "note": "compare mean exact-matched D_eff/D_phys to empirical law",
           }}
    (OUT / "deff_extraction.json").write_text(json.dumps(out, ensure_ascii=False, indent=2),
                                              encoding="utf-8")

    print("=== T5b: D_eff extraction (parallel channel matching) ===")
    print(f"{'D/P':>5} {'nu':>5} {'exact|tpar|':>11} {'Deff/Dphys':>11} {'law':>7} {'cond':>5}")
    for r in rows:
        print(f"{r['D_over_P']:>5.2f} {r['nu_THz']:>5.2f} {r['exact_tpar']:>11.4f} "
              f"{r['Deff_over_Dphys']:>11.3f} {r['law_value']:>7.3f} "
              f"{'ok' if r['well_conditioned'] else 'low':>5}")
    print("\nWell-conditioned (D/P>=0.5), freq-averaged:")
    print(f"{'D/P':>5} {'Deff/Dphys(exact)':>18} {'law':>7} {'gap':>8}")
    for dp, v in by_dp.items():
        print(f"{dp:>5.2f} {v['mean_Deff_over_Dphys']:>18.3f} {v['law']:>7.3f} "
              f"{v['mean_Deff_over_Dphys']-v['law']:>+8.3f}")
    print(f"\nartifacts -> {OUT}")


if __name__ == "__main__":
    main()
