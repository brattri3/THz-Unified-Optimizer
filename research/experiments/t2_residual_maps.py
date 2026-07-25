"""T2 — residual structure maps for the M3 baseline on both datasets.

Produces, per dataset:
  - 2D residual map amp_res_dB over (angle x frequency)
  - phase residual map (rad)
  - Q-Q plot of amp residual vs normal
  - spectral residual profile: mean +/- std of amp_res_dB across angles vs frequency
  - angular-tail profile: mean +/- std of amp_res_dB across frequency vs angle (focus deep shadow)
And a JSON summary quantifying where the structured residual lives.

Read-only on originals; all output under research/results/residuals/.
"""
import sys, json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy import stats

sys.path.insert(0, str(Path(__file__).resolve().parent))
from fit_lib import fit_model, GEOMETRY

OUT = Path(__file__).resolve().parents[1] / "results" / "residuals"
OUT.mkdir(parents=True, exist_ok=True)


def residual_grid_db(exp, theo, valid):
    """Full 2D amp residual in dB, NaN where invalid."""
    amp = 20*np.log10(np.maximum(np.abs(exp), 1e-12)) - 20*np.log10(np.maximum(np.abs(theo), 1e-12))
    amp = np.where(valid, amp, np.nan)
    return amp


def phase_grid(exp, theo, valid):
    ph = np.angle(exp) - np.angle(theo)
    ph = np.arctan2(np.sin(ph), np.cos(ph))
    ph = np.where(valid, ph, np.nan)
    return ph


def summarize(ds):
    m, arr = fit_model(ds, label="M3", free_gamma=False)
    angles = arr["angles"]; freqs = arr["freqs"]  # already in THz
    exp = arr["exp"]; theo = arr["theo"]; valid = arr["valid"]
    amp = residual_grid_db(exp, theo, valid)
    ph = phase_grid(exp, theo, valid)

    # profiles (nan-aware)
    spec_mean = np.nanmean(amp, axis=0); spec_std = np.nanstd(amp, axis=0)
    ang_mean = np.nanmean(amp, axis=1); ang_std = np.nanstd(amp, axis=1)

    # deep-shadow tail: highest 25% of angles
    hi = angles >= np.percentile(angles, 75)
    tail_bias = float(np.nanmean(amp[hi, :]))
    tail_rmse = float(np.sqrt(np.nanmean(amp[hi, :]**2)))
    all_bias = float(np.nanmean(amp))
    all_rmse = float(np.sqrt(np.nanmean(amp**2)))

    # spectral hot spot
    k = int(np.nanargmax(np.abs(spec_mean)))
    spec_hot = {"freq_THz": float(freqs[k]), "mean_bias_dB": float(spec_mean[k])}
    # angular hot spot
    j = int(np.nanargmax(np.abs(ang_mean)))
    ang_hot = {"angle_deg": float(angles[j]), "mean_bias_dB": float(ang_mean[j])}

    # ---- figures ----
    extent = [freqs.min(), freqs.max(), angles.max(), angles.min()]
    # 2D amp map
    fig, ax = plt.subplots(figsize=(7, 4))
    vmax = np.nanpercentile(np.abs(amp), 98)
    im = ax.imshow(amp, aspect="auto", extent=extent, cmap="RdBu_r", vmin=-vmax, vmax=vmax)
    ax.set_xlabel("Frequency (THz)"); ax.set_ylabel("Angle (deg)")
    ax.set_title(f"{ds}: M3 amp residual (dB) [exp-model]")
    fig.colorbar(im, ax=ax, label="dB"); fig.tight_layout()
    fig.savefig(OUT / f"{ds}_amp_map.png", dpi=110); plt.close(fig)

    # phase map
    fig, ax = plt.subplots(figsize=(7, 4))
    pv = np.nanpercentile(np.abs(ph), 98)
    im = ax.imshow(ph, aspect="auto", extent=extent, cmap="PuOr_r", vmin=-pv, vmax=pv)
    ax.set_xlabel("Frequency (THz)"); ax.set_ylabel("Angle (deg)")
    ax.set_title(f"{ds}: M3 phase residual (rad)")
    fig.colorbar(im, ax=ax, label="rad"); fig.tight_layout()
    fig.savefig(OUT / f"{ds}_phase_map.png", dpi=110); plt.close(fig)

    # Q-Q of amp residual
    flat = amp[np.isfinite(amp)]
    fig, ax = plt.subplots(figsize=(4.5, 4.5))
    stats.probplot(flat, dist="norm", plot=ax)
    ax.set_title(f"{ds}: Q-Q amp residual"); fig.tight_layout()
    fig.savefig(OUT / f"{ds}_qq.png", dpi=110); plt.close(fig)

    # spectral profile
    fig, ax = plt.subplots(figsize=(7, 3.5))
    ax.plot(freqs, spec_mean, "b-", label="mean")
    ax.fill_between(freqs, spec_mean-spec_std, spec_mean+spec_std, alpha=0.2, color="b")
    ax.axhline(0, color="k", lw=0.5); ax.set_xlabel("Frequency (THz)")
    ax.set_ylabel("amp residual (dB)"); ax.set_title(f"{ds}: spectral residual profile")
    ax.legend(); fig.tight_layout(); fig.savefig(OUT / f"{ds}_spectral.png", dpi=110); plt.close(fig)

    # angular profile
    fig, ax = plt.subplots(figsize=(7, 3.5))
    ax.plot(angles, ang_mean, "r-o", ms=3, label="mean")
    ax.fill_between(angles, ang_mean-ang_std, ang_mean+ang_std, alpha=0.2, color="r")
    ax.axhline(0, color="k", lw=0.5); ax.set_xlabel("Angle (deg)")
    ax.set_ylabel("amp residual (dB)"); ax.set_title(f"{ds}: angular residual profile")
    ax.legend(); fig.tight_layout(); fig.savefig(OUT / f"{ds}_angular.png", dpi=110); plt.close(fig)

    out = {
        "dataset": ds,
        "n_angles": int(len(angles)), "n_freqs": int(len(freqs)),
        "angle_range_deg": [float(angles.min()), float(angles.max())],
        "freq_range_THz": [float(freqs.min()), float(freqs.max())],
        "amp_all_bias_dB": all_bias, "amp_all_rmse_dB": all_rmse,
        "amp_deepshadow_bias_dB": tail_bias, "amp_deepshadow_rmse_dB": tail_rmse,
        "spectral_hotspot": spec_hot, "angular_hotspot": ang_hot,
        "durbin_watson_amp": m["residual_normality"]["amp"].get("durbin_watson"),
        "shapiro_p_amp": m["residual_normality"]["amp"].get("shapiro_p"),
        "jb_p_amp": m["residual_normality"]["amp"].get("jb_p"),
        "phase_rmse_rad": m["rmse_phase_rad"],
    }
    np.savez_compressed(OUT / f"{ds}_grids.npz",
                        angles=angles, freqs_THz=freqs, amp_res_db=amp, phase_res=ph,
                        spec_mean=spec_mean, spec_std=spec_std,
                        ang_mean=ang_mean, ang_std=ang_std)
    return out


if __name__ == "__main__":
    results = [summarize(ds) for ds in GEOMETRY]
    (OUT / "residuals_summary.json").write_text(
        json.dumps(results, ensure_ascii=False, indent=2), encoding="utf-8")
    for r in results:
        print(f"\n=== {r['dataset']} ===")
        print(f"  amp bias/rmse (all)   : {r['amp_all_bias_dB']:+.3f} / {r['amp_all_rmse_dB']:.3f} dB")
        print(f"  amp bias/rmse (shadow): {r['amp_deepshadow_bias_dB']:+.3f} / {r['amp_deepshadow_rmse_dB']:.3f} dB")
        print(f"  spectral hotspot      : {r['spectral_hotspot']['freq_THz']:.3f} THz "
              f"({r['spectral_hotspot']['mean_bias_dB']:+.3f} dB)")
        print(f"  angular hotspot       : {r['angular_hotspot']['angle_deg']:.1f} deg "
              f"({r['angular_hotspot']['mean_bias_dB']:+.3f} dB)")
        print(f"  DW/shapiro_p/jb_p     : {r['durbin_watson_amp']:.3f} / "
              f"{r['shapiro_p_amp']:.2e} / {r['jb_p_amp']:.2e}")
    print(f"\nartifacts -> {OUT}")
