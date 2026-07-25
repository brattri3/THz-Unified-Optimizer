"""HN10 — grazing-angle aperture apodization vs the POST-M5 angular residual.

T9 found a weak grazing angular slope in the M5 residual. HN10: a geometric
aperture/footprint effect could apodize the throughput at grazing rotation as
A(theta)=1 - a*(1-cos theta) (a free; a=0 no effect, a>0 reduces grazing, a<0
enhances). Applied multiplicatively to the M5 field (NOT a bare cos-theta, which
would wrongly null the real deep-shadow leakage floor). Ablation M5 vs M5+apod on
356att (dense, residual remains) and specac. Read-only; output results/HN10/.
"""
import sys, json
from pathlib import Path
import numpy as np
import lmfit

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import model_m5 as m5
from t6_hn8 import compute_grid_hn8

OUT = HERE.parents[0] / "results" / "HN10"
OUT.mkdir(parents=True, exist_ok=True)
DEFF = {"356att": 4.683, "test_grid_40_20": 11.447, "specac": 6.240}


def residual_apod(params, angles_val, freqs, exp, valid):
    p = params["P_um"].value*1e-6; d = params["D_um"].value*1e-6
    if d >= p:
        return np.ones(int(np.sum(valid))*2)*1e6
    E = compute_grid_hn8(angles_val, freqs, p, d, params["loss_factor"].value,
        params["angle_offset"].value, params["tau_ps"].value, params["gamma"].value,
        params["tau_par_ps"].value, params["eta0"].value, params["eta_exp"].value,
        params["delta0"].value, params["tau_leak"].value)
    a = params["apod"].value
    A = (1.0 - a*(1.0 - np.cos(np.deg2rad(angles_val))))[:, None]
    theo = A * E
    em, tm = exp[valid], theo[valid]
    amp = np.abs(em)-np.abs(tm)
    ph = np.angle(em)-np.angle(tm); ph = np.arctan2(np.sin(ph), np.cos(ph))
    return np.concatenate([amp*m5.W_AMP, ph*m5.W_PHASE])


def fit(dataset, apod):
    angles, freqs, exp, valid = m5.build_experiment_masked(dataset)
    p_um = m5.fit_lib.GEOMETRY[dataset]["P_um"]
    P = lmfit.Parameters()
    P.add("P_um", value=p_um, vary=False)
    P.add("D_um", value=DEFF[dataset], vary=False)
    P.add("loss_factor", value=0.3, min=0.0, max=5.0)
    P.add("gamma", value=2.0, min=0.2, max=5.0)
    P.add("angle_offset", value=0.0, min=-10, max=10)
    P.add("tau_ps", value=0.0, min=-10, max=10)
    P.add("tau_par_ps", value=0.0, min=-5, max=5)
    P.add("eta0", value=0.02, min=0.0, max=1.0)
    P.add("eta_exp", value=1.0, min=-2, max=6)
    P.add("delta0", value=0.0, min=-np.pi, max=np.pi)
    P.add("tau_leak", value=0.1, min=-5, max=5)
    P.add("apod", value=0.0, min=-0.5, max=0.5, vary=apod)
    res = lmfit.Minimizer(residual_apod, P, fcn_args=(angles, freqs, exp, valid)).minimize(method="leastsq")
    return {"dataset": dataset, "apod_on": apod, "aic": float(res.aic), "bic": float(res.bic),
            "apod": float(res.params["apod"].value),
            "apod_err": (float(res.params["apod"].stderr) if res.params["apod"].stderr else None)}


def main():
    rows = []
    for ds in ("356att", "specac"):
        m = fit(ds, False); a = fit(ds, True)
        a["dAIC_vs_M5"] = a["aic"]-m["aic"]
        rows.extend([m, a])
    (OUT / "hn10_apodization.json").write_text(json.dumps(rows, ensure_ascii=False, indent=2), encoding="utf-8")
    print("=== HN10: grazing-angle apodization vs post-M5 residual ===")
    for ds in ("356att", "specac"):
        m = [r for r in rows if r["dataset"] == ds and not r["apod_on"]][0]
        a = [r for r in rows if r["dataset"] == ds and r["apod_on"]][0]
        print(f"\n{ds}: M5 AIC={m['aic']:.1f} | M5+apod AIC={a['aic']:.1f} dAIC={a['dAIC_vs_M5']:+.1f}")
        print(f"   apod a={a['apod']:.4f} +- {a['apod_err']}  (0=no effect; A(90)=1-a)")
    print(f"\nartifacts -> {OUT}")


if __name__ == "__main__":
    main()
