"""T9 — characterize the STRUCTURAL residual remaining after M5.

M5 fixes the deep-shadow deficit but the residual is still autocorrelated (DW<<2).
Decompose what's left (spectral / angular / Q-Q) to localize the leftover mechanism:
 - spectral hotspot ~0.246 THz (water-mask leak vs real line, HN9)?
 - periodic ripple in frequency (Fabry-Perot, HN3)?
 - smooth frequency slope (dispersion / non-Drude gamma, HN5)?
 - angular apodization at grazing angles (HN10)?
Uses the M5 config (D pinned=4.68 dense / 11.45 sparse + fitted leakage) and the
masked builder. Read-only; output results/residual_remain/.
"""
import sys, json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy import stats

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import model_m5 as m5
from t6_hn8 import compute_grid_hn8

OUT = HERE.parents[0] / "results" / "residual_remain"
OUT.mkdir(parents=True, exist_ok=True)


def m5_residual_grid(dataset):
    angles, freqs, exp, valid = m5.build_experiment_masked(dataset)
    # refit M5 (D pinned effective + leakage) to get params
    import lmfit
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
    res = lmfit.Minimizer(m5.residual, P, fcn_args=(angles, freqs, exp, valid)).minimize(method="leastsq")
    pv = res.params
    theo = compute_grid_hn8(angles, freqs, pv["P_um"].value*1e-6, pv["D_um"].value*1e-6,
        pv["loss_factor"].value, pv["angle_offset"].value, pv["tau_ps"].value, pv["gamma"].value,
        pv["tau_par_ps"].value, pv["eta0"].value, pv["eta_exp"].value, pv["delta0"].value, pv["tau_leak"].value)
    ampdb = np.where(valid,
        20*np.log10(np.maximum(np.abs(exp),1e-12))-20*np.log10(np.maximum(np.abs(theo),1e-12)), np.nan)
    return angles, freqs, ampdb, valid


def analyze(dataset):
    angles, freqs, ampdb, valid = m5_residual_grid(dataset)
    spec_mean = np.nanmean(ampdb, axis=0); spec_std = np.nanstd(ampdb, axis=0)
    ang_mean = np.nanmean(ampdb, axis=1); ang_std = np.nanstd(ampdb, axis=1)
    # spectral hotspot
    k = int(np.nanargmax(np.abs(spec_mean)))
    # FP ripple: FFT of spectral residual (detrended), look for periodicity
    good = np.isfinite(spec_mean)
    sm = spec_mean.copy(); sm[~good] = np.nanmean(spec_mean[good])
    sm_dt = sm - np.polyval(np.polyfit(freqs[good], spec_mean[good], 2), freqs)  # remove smooth trend
    ft = np.abs(np.fft.rfft(sm_dt*np.hanning(len(sm_dt))))
    ripple_frac = float(np.max(ft[2:])/ (np.sum(ft[1:])+1e-12))  # dominant ripple prominence
    # smooth slope (dispersion): linear fit of spectral residual
    slope = float(np.polyfit(freqs[good], spec_mean[good], 1)[0])
    # angular: grazing-angle trend (apodization signature)
    hi = angles >= 60
    graz_slope = float(np.polyfit(angles[hi], ang_mean[hi], 1)[0]) if np.sum(hi) >= 2 else None
    flat = ampdb[np.isfinite(ampdb)]
    dw = float(np.sum(np.diff(flat)**2)/np.sum(flat**2))
    sh = float(stats.shapiro(flat[:5000]).pvalue)

    # figures
    fig, axs = plt.subplots(1, 2, figsize=(12, 3.6))
    axs[0].plot(freqs, spec_mean, "b-"); axs[0].fill_between(freqs, spec_mean-spec_std, spec_mean+spec_std, alpha=0.2)
    axs[0].axhline(0, color="k", lw=0.5); axs[0].axvline(0.246, color="r", ls=":", label="0.246 THz")
    axs[0].set_xlabel("Freq (THz)"); axs[0].set_ylabel("M5 residual (dB)"); axs[0].set_title(f"{dataset}: spectral"); axs[0].legend(fontsize=7)
    axs[1].plot(angles, ang_mean, "r-o", ms=3); axs[1].fill_between(angles, ang_mean-ang_std, ang_mean+ang_std, alpha=0.2, color="r")
    axs[1].axhline(0, color="k", lw=0.5); axs[1].set_xlabel("Angle (deg)"); axs[1].set_title(f"{dataset}: angular")
    fig.tight_layout(); fig.savefig(OUT / f"{dataset}_m5_residual.png", dpi=110); plt.close(fig)

    return {"dataset": dataset,
            "spectral_hotspot": {"freq_THz": float(freqs[k]), "bias_dB": float(spec_mean[k])},
            "spectral_slope_dB_per_THz": slope, "ripple_prominence": ripple_frac,
            "grazing_angular_slope_dB_per_deg": graz_slope,
            "durbin_watson": dw, "shapiro_p": sh,
            "spec_mean": spec_mean.tolist(), "freqs": freqs.tolist(),
            "ang_mean": ang_mean.tolist(), "angles": angles.tolist()}


def main():
    res = [analyze(ds) for ds in ("att-11-16-356", "test_grid_40_20")]
    (OUT / "residual_remain.json").write_text(json.dumps(res, ensure_ascii=False, indent=2), encoding="utf-8")
    print("=== T9: structural residual remaining after M5 ===")
    for r in res:
        print(f"\n{r['dataset']}:")
        print(f"  spectral hotspot : {r['spectral_hotspot']['freq_THz']:.3f} THz ({r['spectral_hotspot']['bias_dB']:+.3f} dB)")
        print(f"  spectral slope   : {r['spectral_slope_dB_per_THz']:+.3f} dB/THz (dispersion signature)")
        print(f"  ripple prominence: {r['ripple_prominence']:.3f} (Fabry-Perot signature)")
        print(f"  grazing ang slope: {r['grazing_angular_slope_dB_per_deg']:+.4f} dB/deg (apodization signature)")
        print(f"  DW / shapiro_p   : {r['durbin_watson']:.3f} / {r['shapiro_p']:.2e}")
    print(f"\nartifacts -> {OUT}")


if __name__ == "__main__":
    main()
