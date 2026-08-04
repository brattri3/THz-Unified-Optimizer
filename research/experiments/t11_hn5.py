"""T11 — HN5 non-Drude / physical gamma vs the POST-M5 dispersion residual.

T9 found a smooth -0.62 dB/THz slope in the dense post-M5 residual = a
mis-modelled frequency dependence of loss. M5 fixes the scattering exponent
gamma=2 (Rayleigh-like). HN5 (Laman-Grischkowsky / Naftaly): THz conductivity
departs from Drude (grain boundaries, Fuchs-Sondheimer) -> fractional gamma.
Test: free gamma vs fixed gamma=2. Verdict by dAIC vs M5 + whether the residual
dispersion slope flattens + physical plausibility of gamma (skin ~0.5 .. Rayleigh 4).
Read-only; output results/HN5/.
"""
import sys, json
from pathlib import Path
import numpy as np
import lmfit

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import model_m5 as m5
from t6_hn8 import compute_grid_hn8

OUT = HERE.parents[0] / "results" / "HN5"
OUT.mkdir(parents=True, exist_ok=True)


def fit(dataset, gamma_free):
    angles, freqs, exp, valid = m5.build_experiment_masked(dataset)
    p_um = m5.fit_lib.GEOMETRY[dataset]["P_um"]
    P = lmfit.Parameters()
    P.add("P_um", value=p_um, vary=False)
    P.add("D_um", value=m5.M3_D[dataset], vary=False)
    P.add("loss_factor", value=0.3, min=0.0, max=5.0)
    P.add("gamma", value=2.0, min=0.2, max=5.0, vary=gamma_free)
    P.add("angle_offset", value=0.0, min=-10, max=10)
    P.add("tau_ps", value=0.0, min=-10, max=10)
    P.add("tau_par_ps", value=0.0, min=-5, max=5)
    P.add("eta0", value=0.02, min=0.0, max=1.0)
    P.add("eta_exp", value=1.0, min=-2, max=6)
    P.add("delta0", value=0.0, min=-np.pi, max=np.pi)
    P.add("tau_leak", value=0.1, min=-5, max=5)
    res = lmfit.Minimizer(m5.residual, P, fcn_args=(angles, freqs, exp, valid)).minimize(method="leastsq")
    pv = res.params
    theo = compute_grid_hn8(angles, freqs, pv["P_um"].value*1e-6, pv["D_um"].value*1e-6,
        pv["loss_factor"].value, pv["angle_offset"].value, pv["tau_ps"].value, pv["gamma"].value,
        pv["tau_par_ps"].value, pv["eta0"].value, pv["eta_exp"].value, pv["delta0"].value, pv["tau_leak"].value)
    ampdb = np.where(valid, 20*np.log10(np.maximum(np.abs(exp),1e-12))-20*np.log10(np.maximum(np.abs(theo),1e-12)), np.nan)
    spec_mean = np.nanmean(ampdb, axis=0); good = np.isfinite(spec_mean)
    slope = float(np.polyfit(freqs[good], spec_mean[good], 1)[0])
    flat = ampdb[np.isfinite(ampdb)]
    dw = float(np.sum(np.diff(flat)**2)/np.sum(flat**2))
    return {"dataset": dataset, "gamma_free": gamma_free, "aic": float(res.aic), "bic": float(res.bic),
            "gamma": float(pv["gamma"].value), "loss_factor": float(pv["loss_factor"].value),
            "spectral_slope_dB_per_THz": slope, "durbin_watson": dw,
            "rmse_amp_db": float(np.sqrt(np.nanmean(ampdb[np.isfinite(ampdb)]**2)))}


def main():
    rows = []
    for ds in ("att-11-16-356", "test_grid_40_20"):
        a = fit(ds, False); b = fit(ds, True)
        b["dAIC_vs_M5"] = b["aic"]-a["aic"]
        rows.extend([a, b])
    (OUT / "hn5_nondrude.json").write_text(json.dumps(rows, ensure_ascii=False, indent=2), encoding="utf-8")
    print("=== T11 HN5: free gamma (non-Drude dispersion) vs post-M5 residual ===")
    for ds in ("att-11-16-356", "test_grid_40_20"):
        a = [r for r in rows if r["dataset"] == ds and not r["gamma_free"]][0]
        b = [r for r in rows if r["dataset"] == ds and r["gamma_free"]][0]
        print(f"\n{ds}:")
        print(f"  M5 (gamma=2): AIC={a['aic']:.1f} slope={a['spectral_slope_dB_per_THz']:+.3f}dB/THz rmseDB={a['rmse_amp_db']:.3f} DW={a['durbin_watson']:.3f}")
        print(f"  M5+free gamma: AIC={b['aic']:.1f} dAIC={b['dAIC_vs_M5']:+.1f} gamma={b['gamma']:.3f} "
              f"slope={b['spectral_slope_dB_per_THz']:+.3f} rmseDB={b['rmse_amp_db']:.3f} DW={b['durbin_watson']:.3f}")
    print(f"\nartifacts -> {OUT}")


if __name__ == "__main__":
    main()
