"""T6 test HN8 (user-added) — group delay of the leakage channel.

Refines HN2: the parallel leakage may carry its OWN temporal delay, distinct
from the Blanco-channel tau_par_ps:
    t_par_leak = eta(nu) * exp(j (delta0 + 2*pi*nu*tau_leak))
Test whether a free tau_leak improves the fit over HN2 (constant delta0), and a
model-independent check: peak arrival time t_peak(theta) of the raw sample
waveform vs angle (a shift toward theta->90 would signal a distinct leakage delay).

Because HN7 showed eta is degenerate with D_eff, the identifiability test is run
with D_eff PINNED (fixed at the M3 value) so tau_leak is cleanly isolated, plus a
free-D run for completeness. Ablation ladder on the DENSE 356att (leakage regime):
    M3  (no leakage)
    HN2 (leakage, delta const, tau_leak=0)
    HN8 (leakage + free tau_leak)
Read-only on originals; output under research/results/HN8/.
"""
import sys, json
from pathlib import Path
import numpy as np
import lmfit

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import fit_lib, model_core
sys.path.insert(0, str(HERE.parents[1]))
from unified_optimizer import config, model_blanco

OUT = HERE.parents[0] / "results" / "HN8"
OUT.mkdir(parents=True, exist_ok=True)
W_AMP, W_PHASE = fit_lib.W_AMP, fit_lib.W_PHASE
M3_D = {"att-11-16-356": 4.683, "test_grid_40_20": 11.447}


def compute_grid_hn8(angles_deg, freqs_thz, p, d, loss_factor, angle_offset, tau_ps,
                     gamma, tau_par_ps, eta0, eta_exp, delta0, tau_leak,
                     N=15, use_drude=True):
    """HN8 grid (leakage with group delay tau_leak) = model_core.compute_grid.

    Was an inline per-frequency Blanco loop; now a thin alias over the cached
    core (model_m5 and t20 import this name). Bit-identical to the old copy --
    asserted by verify_model_core.py.
    """
    return model_core.compute_grid(
        angles_deg, freqs_thz, p, d, loss_factor, angle_offset, tau_ps,
        gamma=gamma, N=N, use_drude=use_drude, tau_par_ps=tau_par_ps,
        eta0=eta0, eta_exp=eta_exp, delta0=delta0, tau_leak=tau_leak)


def residual(params, angles_val, freqs, exp, valid):
    return model_core.residual_from_params(params, angles_val, freqs, exp, valid,
                                           w_amp=W_AMP, w_phase=W_PHASE)


def fit(dataset, mode, fix_D=True):
    """mode in {m3, hn2, hn8}"""
    geo = fit_lib.GEOMETRY[dataset]
    from unified_optimizer.data_manager import DataManager
    data = DataManager(config.DATA_DIR).get_data_for_dataset(dataset)
    angles_val, freqs, exp, valid = fit_lib.build_experiment(data)
    p_um = geo["P_um"]
    P = lmfit.Parameters()
    P.add("P_um", value=p_um, vary=False)
    if fix_D:
        P.add("D_um", value=M3_D[dataset], vary=False)
    else:
        P.add("D_um", value=model_blanco.estimate_deff_initial(p_um, geo["D_phys_um"]),
              min=0.1, max=p_um-0.5)
    P.add("loss_factor", value=0.3, min=0.0, max=5.0)
    P.add("gamma", value=2.0, vary=False)
    P.add("angle_offset", value=0.0, min=-10, max=10)
    P.add("tau_ps", value=0.0, min=-10, max=10)
    P.add("tau_par_ps", value=0.0, min=-5, max=5)
    leak = mode in ("hn2", "hn8")
    P.add("eta0", value=0.02 if leak else 0.0, min=0.0, max=1.0, vary=leak)
    P.add("eta_exp", value=1.0, min=-2, max=6, vary=leak)
    P.add("delta0", value=0.0, min=-np.pi, max=np.pi, vary=leak)
    P.add("tau_leak", value=0.0, min=-5, max=5, vary=(mode == "hn8"))
    res = lmfit.Minimizer(residual, P, fcn_args=(angles_val, freqs, exp, valid)).minimize(method="leastsq")
    return {"dataset": dataset, "mode": mode, "fix_D": fix_D,
            "aic": float(res.aic), "bic": float(res.bic), "redchi": float(res.redchi),
            "nvarys": int(res.nvarys),
            "params": {k: float(v.value) for k, v in res.params.items()}}


def t_peak_profile(dataset):
    """Model-independent: peak arrival time of the raw sample waveform vs angle."""
    from unified_optimizer.data_manager import DataManager
    data = DataManager(config.DATA_DIR).get_data_for_dataset(dataset)
    angles = sorted(data.keys())
    if config.ANGLES_LIMIT_2D:
        lo, hi = config.ANGLES_LIMIT_2D; angles = [a for a in angles if lo <= a <= hi]
    prof = []
    for a in angles:
        t_s, E_s, t_b, E_b = data[a]
        t_s = np.asarray(t_s); E_s = np.asarray(E_s)
        # centroid of energy |E|^2 (robust) and argmax
        w = E_s**2
        t_centroid = float(np.sum(t_s*w)/np.sum(w)) if np.sum(w) > 0 else None
        t_argmax = float(t_s[int(np.argmax(np.abs(E_s)))])
        prof.append({"angle": float(a), "t_centroid_ps": t_centroid, "t_peak_ps": t_argmax})
    return prof


def main():
    rows = []
    for ds in ("att-11-16-356", "test_grid_40_20"):
        for fixD in (True, False):
            m3 = fit(ds, "m3", fixD); hn2 = fit(ds, "hn2", fixD); hn8 = fit(ds, "hn8", fixD)
            base = m3["aic"]
            for r in (m3, hn2, hn8):
                r["dAIC_vs_M3"] = r["aic"]-base
            hn8["dAIC_vs_HN2"] = hn8["aic"]-hn2["aic"]
            rows.extend([m3, hn2, hn8])
    tpk = {ds: t_peak_profile(ds) for ds in ("att-11-16-356", "test_grid_40_20")}
    (OUT / "hn8_ablation.json").write_text(json.dumps({"ablation": rows, "t_peak": tpk},
                                           ensure_ascii=False, indent=2), encoding="utf-8")

    print("=== T6 HN8: leakage group-delay tau_leak ===")
    print(f"{'dataset':<16}{'fixD':<6}{'mode':<5}{'AIC':>10}{'dAIC_M3':>9}{'dAIC_HN2':>10}{'tau_leak':>9}")
    for r in rows:
        print(f"{r['dataset']:<16}{str(r['fix_D']):<6}{r['mode']:<5}{r['aic']:>10.1f}"
              f"{r['dAIC_vs_M3']:>9.1f}{r.get('dAIC_vs_HN2', 0.0):>10.1f}"
              f"{r['params']['tau_leak']:>9.4f}")
    print("\nModel-independent peak-time shift (deep shadow vs theta~0):")
    for ds, pr in tpk.items():
        p0 = pr[0]; p90 = pr[-1]
        print(f"  {ds:<16} theta={p0['angle']:.0f}: t_peak={p0['t_peak_ps']:.3f}ps  "
              f"theta={p90['angle']:.0f}: t_peak={p90['t_peak_ps']:.3f}ps  "
              f"shift={p90['t_peak_ps']-p0['t_peak_ps']:+.3f}ps")
    print(f"\nartifacts -> {OUT}")


if __name__ == "__main__":
    main()
