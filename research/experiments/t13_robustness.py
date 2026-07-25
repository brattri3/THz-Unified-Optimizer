"""T13 — robustness / overfitting guard (only 2 datasets -> be careful).

(1) HELD-OUT-ANGLE cross-validation: fit on non-shadow angles (theta<75), PREDICT
    the deep-shadow angles (theta>=75) that drive the leakage. If M5 (with leakage)
    predicts held-out deep shadow BETTER than M3 (no leakage), the leakage is a
    genuine predictive feature, not overfitting.
(2) BOOTSTRAP CI on eta0 / tau_leak (D pinned): resample frequency columns with
    replacement, refit, report 16-84 percentile intervals.
Dense 356att (leakage regime) + sparse control. Read-only; results/robustness/.
"""
import sys, json
from pathlib import Path
import numpy as np
import lmfit

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import model_m5 as m5
from t6_hn8 import compute_grid_hn8

OUT = HERE.parents[0] / "results" / "robustness"
OUT.mkdir(parents=True, exist_ok=True)


def params(dataset, leakage, Dfix):
    p_um = m5.fit_lib.GEOMETRY[dataset]["P_um"]
    P = lmfit.Parameters()
    P.add("P_um", value=p_um, vary=False)
    P.add("D_um", value=Dfix, vary=False)
    P.add("loss_factor", value=0.3, min=0.0, max=5.0)
    P.add("gamma", value=2.0, min=0.2, max=5.0)   # free gamma (HN5)
    P.add("angle_offset", value=0.0, min=-10, max=10)
    P.add("tau_ps", value=0.0, min=-10, max=10)
    P.add("tau_par_ps", value=0.0, min=-5, max=5)
    P.add("eta0", value=0.02 if leakage else 0.0, min=0.0, max=1.0, vary=leakage)
    P.add("eta_exp", value=1.0, min=-2, max=6, vary=leakage)
    P.add("delta0", value=0.0, min=-np.pi, max=np.pi, vary=leakage)
    P.add("tau_leak", value=0.1, min=-5, max=5, vary=leakage)
    return P


def predict_db(pv, angles, freqs, exp, mask):
    theo = compute_grid_hn8(angles, freqs, pv["P_um"].value*1e-6, pv["D_um"].value*1e-6,
        pv["loss_factor"].value, pv["angle_offset"].value, pv["tau_ps"].value, pv["gamma"].value,
        pv["tau_par_ps"].value, pv["eta0"].value, pv["eta_exp"].value, pv["delta0"].value, pv["tau_leak"].value)
    d = (20*np.log10(np.maximum(np.abs(exp),1e-12))-20*np.log10(np.maximum(np.abs(theo),1e-12)))
    d = d[mask]; d = d[np.isfinite(d)]
    return float(np.sqrt(np.mean(d**2))), float(np.mean(d))


def heldout(dataset):
    angles, freqs, exp, valid = m5.build_experiment_masked(dataset)
    Dfix = m5.M3_D[dataset]
    train = valid & (angles[:, None] < 75)
    test = valid & (angles[:, None] >= 75)
    res = {}
    for name, leak in (("M3_noleak", False), ("M5_leak", True)):
        r = lmfit.Minimizer(m5.residual, params(dataset, leak, Dfix),
                            fcn_args=(angles, freqs, exp, train)).minimize(method="leastsq")
        rmse_tr, _ = predict_db(r.params, angles, freqs, exp, train)
        rmse_te, bias_te = predict_db(r.params, angles, freqs, exp, test)
        res[name] = {"train_rmse_db": rmse_tr, "heldout_shadow_rmse_db": rmse_te,
                     "heldout_shadow_bias_db": bias_te, "eta0": float(r.params["eta0"].value)}
    return res


def bootstrap(dataset, nboot=40, seed=1):
    angles, freqs, exp, valid = m5.build_experiment_masked(dataset)
    Dfix = m5.M3_D[dataset]
    rng = np.random.default_rng(seed)
    nf = len(freqs); ets, tls, gms = [], [], []
    for _ in range(nboot):
        cols = rng.integers(0, nf, nf)         # resample frequency columns w/ replacement
        fb = freqs[cols]; eb = exp[:, cols]; vb = valid[:, cols]
        r = lmfit.Minimizer(m5.residual, params(dataset, True, Dfix),
                            fcn_args=(angles, fb, eb, vb)).minimize(method="leastsq")
        ets.append(r.params["eta0"].value); tls.append(r.params["tau_leak"].value); gms.append(r.params["gamma"].value)
    def ci(a):
        a = np.array(a); return [float(np.percentile(a, 16)), float(np.median(a)), float(np.percentile(a, 84))]
    return {"eta0_ci_16_50_84": ci(ets), "tau_leak_ci_16_50_84": ci(tls), "gamma_ci_16_50_84": ci(gms), "nboot": nboot}


def main():
    out = {}
    for ds in ("356att", "test_grid_40_20"):
        out[ds] = {"heldout_angle": heldout(ds), "bootstrap": bootstrap(ds)}
    (OUT / "robustness.json").write_text(json.dumps(out, ensure_ascii=False, indent=2), encoding="utf-8")

    print("=== T13: robustness (held-out-angle + bootstrap) ===")
    for ds in out:
        ho = out[ds]["heldout_angle"]; bs = out[ds]["bootstrap"]
        print(f"\n{ds}:")
        print(f"  HELD-OUT deep shadow (train theta<75, predict theta>=75):")
        for name in ("M3_noleak", "M5_leak"):
            h = ho[name]
            print(f"    {name:<10} train_rmse={h['train_rmse_db']:.3f}dB  heldout_shadow_rmse={h['heldout_shadow_rmse_db']:.3f}dB "
                  f"bias={h['heldout_shadow_bias_db']:+.3f} eta0={h['eta0']:.4f}")
        imp = ho["M3_noleak"]["heldout_shadow_rmse_db"] - ho["M5_leak"]["heldout_shadow_rmse_db"]
        print(f"    -> leakage improves held-out shadow rmse by {imp:+.3f} dB (positive = predictive, not overfit)")
        print(f"  BOOTSTRAP CI (16/50/84): eta0={bs['eta0_ci_16_50_84']}  tau_leak={bs['tau_leak_ci_16_50_84']}  gamma={bs['gamma_ci_16_50_84']}")
    print(f"\nartifacts -> {OUT}")


if __name__ == "__main__":
    main()
