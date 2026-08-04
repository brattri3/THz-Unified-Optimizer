"""HN9 — is the post-M5 spectral hotspot ~0.215 THz a real feature or a
band-edge / water-mask leak?

The hotspot (T9: 0.215 THz, +1.74 dB, dense) sits right at the low-frequency band
edge F_MIN=0.2 THz, where THz power is low and relative noise high. Test: refit M5
with progressively raised low-frequency cut; if the hotspot/residual improves and
vanishes, it is a band-edge/low-SNR artifact (not a physical line). Also list
known THz water lines near 0.215 (there are none strong) and the auto-water-mask
coverage. Read-only; output results/HN9/.
"""
import sys, json
from pathlib import Path
import numpy as np
import lmfit

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import model_m5 as m5
from t6_hn8 import compute_grid_hn8

OUT = HERE.parents[0] / "results" / "HN9"
OUT.mkdir(parents=True, exist_ok=True)
# strong rotational water-vapor lines in the THz band (THz), for reference
WATER_LINES = [0.557, 0.752, 0.988, 1.097, 1.113, 1.163, 1.207, 1.229]


def fit_with_fmin(dataset, fmin):
    angles, freqs, exp, valid = m5.build_experiment_masked(dataset)
    valid = valid & (freqs[None, :] >= fmin)          # extra low-frequency cut
    p_um = m5.fit_lib.GEOMETRY[dataset]["P_um"]
    P = lmfit.Parameters()
    P.add("P_um", value=p_um, vary=False)
    P.add("D_um", value=m5.M3_D[dataset], vary=False)
    P.add("loss_factor", value=0.3, min=0.0, max=5.0)
    P.add("gamma", value=2.0, min=0.2, max=5.0)
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
    spec_mean = np.nanmean(ampdb, axis=0)
    good = np.isfinite(spec_mean)
    k = int(np.nanargmax(np.abs(spec_mean)))
    flat = ampdb[np.isfinite(ampdb)]
    dw = float(np.sum(np.diff(flat)**2)/np.sum(flat**2))
    return {"fmin": fmin, "npts": int(np.sum(valid)),
            "hotspot_THz": float(freqs[k]), "hotspot_dB": float(spec_mean[k]),
            "rmse_db": float(np.sqrt(np.nanmean(flat**2))), "durbin_watson": dw}


def main():
    ds = "att-11-16-356"
    rows = [fit_with_fmin(ds, fm) for fm in (0.20, 0.25, 0.30, 0.35, 0.40)]
    near = [w for w in WATER_LINES if abs(w-0.215) < 0.05]
    out = {"dataset": ds, "fmin_scan": rows, "water_lines_near_0.215THz": near,
           "note": "no strong rotational water line within 50 GHz of 0.215 THz" if not near else near}
    (OUT / "hn9_spectral_line.json").write_text(json.dumps(out, ensure_ascii=False, indent=2), encoding="utf-8")
    print("=== HN9: is the 0.215 THz hotspot real or a band-edge artifact? (356att) ===")
    print(f"{'F_MIN':>7}{'npts':>7}{'hotspot_THz':>13}{'hotspot_dB':>12}{'rmse_dB':>9}{'DW':>7}")
    for r in rows:
        print(f"{r['fmin']:>7.2f}{r['npts']:>7}{r['hotspot_THz']:>13.3f}{r['hotspot_dB']:>12.3f}"
              f"{r['rmse_db']:>9.3f}{r['durbin_watson']:>7.3f}")
    print(f"\nStrong water lines within 50 GHz of 0.215 THz: {near if near else 'NONE'}")
    print(f"artifacts -> {OUT}")


if __name__ == "__main__":
    main()
