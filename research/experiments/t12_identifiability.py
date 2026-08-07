"""T12 — identifiability: map the D_eff<->eta degeneracy and separate SAMPLE
leakage from DETECTOR cross-pol.

(1) Parameter correlation matrix from the free M5 fit (which params trade off).
(2) Profile likelihood: scan D_um on a grid, at each pin D and refit eta0(+others);
    the chisqr(D) curve maps the degeneracy valley (flat => unidentified).
(3) Sample vs detector: sample leakage enters as sin^2(theta)*t_par (grows with
    grid rotation); a detector cross-pol eps_cross is angle-INDEPENDENT. Add an
    angle-independent complex offset eps_cross to E_out and check whether it is
    separately identified from eta (correlation + value), using the angular signature.
Dense 356att (leakage regime). Read-only; output results/identifiability/.
"""
import sys, json
from pathlib import Path
import numpy as np
import lmfit

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import model_m5 as m5
from t6_hn8 import compute_grid_hn8

OUT = HERE.parents[0] / "results" / "identifiability"
OUT.mkdir(parents=True, exist_ok=True)
DS = "att-11-16-356"


def base_params(Dfix=None, gamma_free=True):
    p_um = m5.fit_lib.GEOMETRY[DS]["P_um"]
    P = lmfit.Parameters()
    P.add("P_um", value=p_um, vary=False)
    if Dfix is None:
        P.add("D_um", value=m5.M3_D[DS], min=0.5, max=p_um-0.5)
    else:
        P.add("D_um", value=Dfix, vary=False)
    P.add("loss_factor", value=0.3, min=0.0, max=5.0)
    P.add("gamma", value=2.0, min=0.2, max=5.0, vary=gamma_free)
    P.add("angle_offset", value=0.0, min=-10, max=10)
    P.add("tau_ps", value=0.0, min=-10, max=10)
    P.add("tau_par_ps", value=0.0, min=-5, max=5)
    P.add("eta0", value=0.02, min=0.0, max=1.0)
    P.add("eta_exp", value=1.0, min=-2, max=6)
    P.add("delta0", value=0.0, min=-np.pi, max=np.pi)
    P.add("tau_leak", value=0.1, min=-5, max=5)
    return P


def main():
    angles, freqs, exp, valid = m5.build_experiment_masked(DS)
    args = (angles, freqs, exp, valid)

    # (1) free fit + correlation matrix
    res = lmfit.Minimizer(m5.residual, base_params(), fcn_args=args).minimize(method="leastsq")
    correl = {}
    if res.params["D_um"].correl:
        for k, v in res.params["D_um"].correl.items():
            correl[f"D_um~{k}"] = float(v)
    free_fit = {"D_um": float(res.params["D_um"].value), "eta0": float(res.params["eta0"].value),
                "gamma": float(res.params["gamma"].value), "redchi": float(res.redchi),
                "D_um_stderr": (float(res.params["D_um"].stderr) if res.params["D_um"].stderr else None),
                "eta0_stderr": (float(res.params["eta0"].stderr) if res.params["eta0"].stderr else None),
                "correl_with_D": correl}

    # (2) profile likelihood chisqr(D)
    Dgrid = np.linspace(2.0, 11.0, 19)
    profile = []
    for Dv in Dgrid:
        r = lmfit.Minimizer(m5.residual, base_params(Dfix=Dv), fcn_args=args).minimize(method="leastsq")
        profile.append({"D_um": float(Dv), "chisqr": float(r.chisqr), "aic": float(r.aic),
                        "eta0": float(r.params["eta0"].value)})
    chis = np.array([p["chisqr"] for p in profile]); dmin = Dgrid[int(np.argmin(chis))]
    # curvature / flatness: range of chisqr over grid relative to min
    flatness = float((chis.max()-chis.min())/chis.min())

    # (3) add angle-independent detector cross-pol eps_cross, test separability
    def resid_eps(params, angles_val, fr, ex, val):
        p = params["P_um"].value*1e-6; d = params["D_um"].value*1e-6
        if d >= p: return np.ones(int(np.sum(val))*2)*1e6
        E = compute_grid_hn8(angles_val, fr, p, d, params["loss_factor"].value,
            params["angle_offset"].value, params["tau_ps"].value, params["gamma"].value,
            params["tau_par_ps"].value, params["eta0"].value, params["eta_exp"].value,
            params["delta0"].value, params["tau_leak"].value)
        epsr = params["eps_re"].value; epsi = params["eps_im"].value
        E = E + (epsr + 1j*epsi)                       # angle-INDEPENDENT detector offset
        em, tm = ex[val], E[val]
        amp = np.abs(em)-np.abs(tm); ph = np.angle(em)-np.angle(tm)
        ph = np.arctan2(np.sin(ph), np.cos(ph))
        return np.concatenate([amp*m5.W_AMP, ph*m5.W_PHASE])
    P2 = base_params()
    P2.add("eps_re", value=0.0, min=-0.2, max=0.2)
    P2.add("eps_im", value=0.0, min=-0.2, max=0.2)
    r2 = lmfit.Minimizer(resid_eps, P2, fcn_args=args).minimize(method="leastsq")
    eps_mag = float(np.hypot(r2.params["eps_re"].value, r2.params["eps_im"].value))
    sep = {"eps_cross_mag": eps_mag, "eta0_with_eps": float(r2.params["eta0"].value),
           "aic_with_eps": float(r2.aic), "aic_no_eps": float(res.aic),
           "dAIC_eps": float(r2.aic-res.aic),
           "eta0_eps_correl": (float(r2.params["eta0"].correl.get("eps_re")) if (r2.params["eta0"].correl and "eps_re" in r2.params["eta0"].correl) else None)}

    out = {"dataset": DS, "free_fit": free_fit, "profile_Dscan": profile,
           "profile_min_D": float(dmin), "profile_flatness": flatness, "sample_vs_detector": sep}
    (OUT / "identifiability.json").write_text(json.dumps(out, ensure_ascii=False, indent=2), encoding="utf-8")

    print("=== T12: identifiability (356att) ===")
    print(f"(1) free fit: D={free_fit['D_um']:.3f}(+-{free_fit['D_um_stderr']}) "
          f"eta0={free_fit['eta0']:.4f}(+-{free_fit['eta0_stderr']}) gamma={free_fit['gamma']:.3f}")
    print(f"    correlations with D_um: {json.dumps({k: round(v,3) for k,v in correl.items()})}")
    print(f"(2) profile chisqr(D): min at D={dmin:.2f}um, flatness (max-min)/min = {flatness:.4f}")
    print(f"    eta0 rises as D falls: " + ", ".join(f"D{p['D_um']:.1f}:eta{p['eta0']:.3f}" for p in profile[::4]))
    print(f"(3) sample vs detector: eps_cross_mag={eps_mag:.4f} dAIC(add eps)={sep['dAIC_eps']:+.1f} "
          f"eta0(with eps)={sep['eta0_with_eps']:.4f} correl(eta0,eps_re)={sep['eta0_eps_correl']}")
    print(f"\nartifacts -> {OUT}")


if __name__ == "__main__":
    main()
