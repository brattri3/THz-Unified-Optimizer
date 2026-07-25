"""T1 — baseline: fit M0 (pure Blanco) and M3 (Blanco+Drude+Scatter+tau_par) on both datasets."""
import sys, json
from pathlib import Path
import numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parent))
from fit_lib import fit_model, GEOMETRY

OUT = Path(__file__).resolve().parents[1] / "results" / "baseline"
OUT.mkdir(parents=True, exist_ok=True)

MODELS = {
    "M0": dict(use_drude=False, use_scattering=False, free_gamma=False, use_tau_par=False),
    "M3": dict(use_drude=True,  use_scattering=True,  free_gamma=False, use_tau_par=True),
}

summary = []
for ds in GEOMETRY:
    aics = {}
    for label, cfg in MODELS.items():
        m, arr = fit_model(ds, label=label, **cfg)
        (OUT / f"{ds}_{label}.json").write_text(json.dumps(m, ensure_ascii=False, indent=2), encoding="utf-8")
        np.savez_compressed(OUT / f"{ds}_{label}_arrays.npz",
                            angles=arr["angles"], freqs=arr["freqs"],
                            amp_res_db=arr["amp_res_db"])
        aics[label] = m["aic"]
        summary.append({"dataset": ds, "model": label, "D_um": m["params"]["D_um"]["value"],
                        "redchi": m["redchi"], "aic": m["aic"],
                        "rmse_amp_db": m["rmse_amp_db"], "rmse_phase": m["rmse_phase_rad"],
                        "tau_par": m["params"]["tau_par_ps"]["value"],
                        "dw_amp": m["residual_normality"]["amp"].get("durbin_watson")})
    print(f"\n=== {ds} ===  dAIC(M3-M0) = {aics['M3']-aics['M0']:+.1f}")

(OUT / "summary.json").write_text(json.dumps(summary, ensure_ascii=False, indent=2), encoding="utf-8")
print("\n{:<16} {:<4} {:>8} {:>10} {:>10} {:>10} {:>9} {:>8}".format(
      "dataset","mdl","D_um","redchi","AIC","rmse_dB","tau_par","DW_amp"))
for s in summary:
    print("{:<16} {:<4} {:>8.3f} {:>10.5f} {:>10.1f} {:>10.3f} {:>9.4f} {:>8.3f}".format(
          s["dataset"], s["model"], s["D_um"], s["redchi"], s["aic"],
          s["rmse_amp_db"], s["tau_par"], s["dw_amp"] or float("nan")))
print(f"\nartifacts -> {OUT}")
