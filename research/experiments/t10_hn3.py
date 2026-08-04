"""T10 — HN3 Fabry-Perot ripple tested against the POST-M5 residual.

A finite substrate/inter-grid gap adds a weak etalon ripple to the transmission:
    F(nu) = 1/(1 - r * exp(-j 2 pi nu tau_rt)),   r small, tau_rt = round-trip delay.
Applied multiplicatively to the M5 field. Ablation M5 vs M5+FP on both datasets;
verdict by dAIC vs M5 + whether tau_rt maps to a plausible physical thickness.
T9 found FP ripple prominence only 0.083 (weak) -> expect small effect.
Read-only; output results/HN3/.
"""
import sys, json
from pathlib import Path
import numpy as np
import lmfit

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import model_m5 as m5
from t6_hn8 import compute_grid_hn8

OUT = HERE.parents[0] / "results" / "HN3"
OUT.mkdir(parents=True, exist_ok=True)


def residual_fp(params, angles_val, freqs, exp, valid):
    p = params["P_um"].value*1e-6; d = params["D_um"].value*1e-6
    if d >= p:
        return np.ones(int(np.sum(valid))*2)*1e6
    E = compute_grid_hn8(angles_val, freqs, p, d, params["loss_factor"].value,
        params["angle_offset"].value, params["tau_ps"].value, params["gamma"].value,
        params["tau_par_ps"].value, params["eta0"].value, params["eta_exp"].value,
        params["delta0"].value, params["tau_leak"].value)
    r = params["r_fp"].value; tau = params["tau_fp"].value
    F = 1.0/(1.0 - r*np.exp(-1j*2*np.pi*freqs*tau))[None, :]
    theo = E * F
    em, tm = exp[valid], theo[valid]
    amp = np.abs(em)-np.abs(tm)
    ph = np.angle(em)-np.angle(tm); ph = np.arctan2(np.sin(ph), np.cos(ph))
    return np.concatenate([amp*m5.W_AMP, ph*m5.W_PHASE])


def fit(dataset, fp):
    angles, freqs, exp, valid = m5.build_experiment_masked(dataset)
    p_um = m5.fit_lib.GEOMETRY[dataset]["P_um"]
    P = lmfit.Parameters()
    P.add("P_um", value=p_um, vary=False)
    P.add("D_um", value=m5.M3_D[dataset], vary=False)
    P.add("loss_factor", value=0.3, min=0.0, max=5.0)
    P.add("gamma", value=2.0, vary=False)
    P.add("angle_offset", value=0.0, min=-10, max=10)
    P.add("tau_ps", value=0.0, min=-10, max=10)
    P.add("tau_par_ps", value=0.0, min=-5, max=5)
    P.add("eta0", value=0.02, min=0.0, max=1.0)
    P.add("eta_exp", value=1.0, min=-2, max=6)
    P.add("delta0", value=0.0, min=-np.pi, max=np.pi)
    P.add("tau_leak", value=0.1, min=-5, max=5)
    P.add("r_fp", value=0.05 if fp else 0.0, min=0.0, max=0.6, vary=fp)
    P.add("tau_fp", value=5.0, min=0.5, max=25.0, vary=fp)   # ps round-trip
    res = lmfit.Minimizer(residual_fp, P, fcn_args=(angles, freqs, exp, valid)).minimize(method="leastsq")
    return {"dataset": dataset, "fp": fp, "aic": float(res.aic), "bic": float(res.bic),
            "r_fp": float(res.params["r_fp"].value), "tau_fp_ps": float(res.params["tau_fp"].value),
            "nvarys": int(res.nvarys)}


def main():
    rows = []
    for ds in ("att-11-16-356", "test_grid_40_20"):
        m = fit(ds, False); f = fit(ds, True)
        f["dAIC_vs_M5"] = f["aic"]-m["aic"]
        # implied substrate thickness for n=3.4 (HR-Si): L = c*tau/(2n)
        f["implied_L_um_Si"] = 3e8*f["tau_fp_ps"]*1e-12/(2*3.4)*1e6
        rows.extend([m, f])
    (OUT / "hn3_fabryperot.json").write_text(json.dumps(rows, ensure_ascii=False, indent=2), encoding="utf-8")
    print("=== T10 HN3: Fabry-Perot vs post-M5 residual ===")
    for ds in ("att-11-16-356", "test_grid_40_20"):
        m = [r for r in rows if r["dataset"] == ds and not r["fp"]][0]
        f = [r for r in rows if r["dataset"] == ds and r["fp"]][0]
        print(f"\n{ds}: M5 AIC={m['aic']:.1f} | M5+FP AIC={f['aic']:.1f} dAIC={f['dAIC_vs_M5']:+.1f}")
        print(f"   r_fp={f['r_fp']:.4f}  tau_fp={f['tau_fp_ps']:.2f}ps  -> implied Si thickness {f['implied_L_um_Si']:.0f} um "
              f"(ripple period {1/f['tau_fp_ps']:.3f} THz)")
    print(f"\nartifacts -> {OUT}")


if __name__ == "__main__":
    main()
