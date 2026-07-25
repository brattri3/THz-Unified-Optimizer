"""T15c/HN11 — is the non-Drude gamma a SHARED material constant across samples?

T11 found gamma=0.91 (356att dense) vs 2.08 (test_grid sparse). With the 3rd sample
(specac) we test identifiability of a shared physical constant:
 (1) per-sample free-gamma fit (M3-style: D free, NO leakage -> gamma is the clean
     material/scattering exponent);
 (2) scan a SHARED gamma across a grid, refit each sample's geometry at fixed gamma,
     sum chisqr -> joint-optimal shared gamma and the AIC cost of forcing sharing.
If sharing costs little, gamma is a genuine shared constant; if a lot, samples
differ (microstructure). Read-only; output results/HN11/.
"""
import sys, json
from pathlib import Path
import numpy as np
import lmfit

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import model_m5 as m5
from t6_hn8 import compute_grid_hn8

OUT = HERE.parents[0] / "results" / "HN11"
OUT.mkdir(parents=True, exist_ok=True)
DATASETS = ["356att", "test_grid_40_20", "specac"]


def fit_one(dataset, gamma_fixed=None):
    angles, freqs, exp, valid = m5.build_experiment_masked(dataset)
    p_um = m5.fit_lib.GEOMETRY[dataset]["P_um"]
    d_init = m5.fit_lib.GEOMETRY[dataset]["D_phys_um"]
    P = lmfit.Parameters()
    P.add("P_um", value=p_um, vary=False)
    P.add("D_um", value=min(d_init, p_um-1.0), min=0.1, max=p_um-0.5)
    P.add("loss_factor", value=0.3, min=0.0, max=5.0)
    P.add("gamma", value=(gamma_fixed if gamma_fixed else 1.5),
          min=0.2, max=5.0, vary=(gamma_fixed is None))
    P.add("angle_offset", value=0.0, min=-10, max=10)
    P.add("tau_ps", value=0.0, min=-10, max=10)
    P.add("tau_par_ps", value=0.0, min=-5, max=5)
    P.add("eta0", value=0.0, vary=False)   # no leakage: clean gamma
    P.add("eta_exp", value=1.0, vary=False)
    P.add("delta0", value=0.0, vary=False)
    P.add("tau_leak", value=0.0, vary=False)
    res = lmfit.Minimizer(m5.residual, P, fcn_args=(angles, freqs, exp, valid)).minimize(method="leastsq")
    return {"gamma": float(res.params["gamma"].value), "D_um": float(res.params["D_um"].value),
            "chisqr": float(res.chisqr), "ndata": int(res.ndata), "aic": float(res.aic)}


def main():
    # (1) per-sample free gamma
    per = {ds: fit_one(ds) for ds in DATASETS}
    sum_aic_free = sum(per[ds]["aic"] for ds in DATASETS)

    # (2) shared-gamma scan
    grid = [0.6, 0.8, 1.0, 1.2, 1.5, 1.8, 2.1]
    scan = []
    for g in grid:
        fits = {ds: fit_one(ds, gamma_fixed=g) for ds in DATASETS}
        tot_chi = sum(fits[ds]["chisqr"] for ds in DATASETS)
        tot_aic = sum(fits[ds]["aic"] for ds in DATASETS)
        scan.append({"gamma_shared": g, "total_chisqr": tot_chi, "total_aic": tot_aic,
                     "D_um": {ds: fits[ds]["D_um"] for ds in DATASETS}})
    best = min(scan, key=lambda s: s["total_aic"])
    # +2 free params saved by sharing (2 gammas fixed to 1) -> compare AIC fairly
    cost = best["total_aic"] - sum_aic_free

    out = {"per_sample_free_gamma": per, "sum_aic_free": sum_aic_free,
           "shared_gamma_scan": scan, "best_shared_gamma": best["gamma_shared"],
           "shared_total_aic": best["total_aic"], "aic_cost_of_sharing": cost}
    (OUT / "hn11_shared_gamma.json").write_text(json.dumps(out, ensure_ascii=False, indent=2), encoding="utf-8")

    print("=== T15c/HN11: is gamma a shared material constant? ===")
    print("per-sample free gamma:")
    for ds in DATASETS:
        print(f"  {ds:<16} gamma={per[ds]['gamma']:.3f}  D_eff={per[ds]['D_um']:.3f}  AIC={per[ds]['aic']:.1f}")
    print(f"\nshared-gamma scan (sum AIC over 3 samples):")
    for s in scan:
        print(f"  gamma={s['gamma_shared']:.1f}  total_AIC={s['total_aic']:.1f}")
    print(f"\nbest shared gamma = {best['gamma_shared']:.1f}; sum-AIC(free)={sum_aic_free:.1f} "
          f"sum-AIC(shared)={best['total_aic']:.1f}; cost of sharing = {cost:+.1f}")
    print(f"artifacts -> {OUT}")


if __name__ == "__main__":
    main()
