"""T8 — consolidated model M5: resolve the D_eff<->eta degeneracy.

Phase-1 left an ambiguity: the deep-shadow parallel channel can be fit EITHER by a
small (collapsed) D_eff OR by explicit leakage eta -- same physics, degenerate.
M5 resolves it by PINNING D_eff (external info: lattice-sum/physical, T5b) and
giving the deep shadow to an explicit small leakage eta(nu)*e^{j(delta0+2pi nu tau_leak)}
(HN7 scale, HN8 delay). Also adds the dynamic-range NOISE MASK to the experiment
builder (was missing; critical in deep shadow).

Decisive sub-question tested here: does pinning D_eff at the PHYSICAL value already
remove the T2 deep-shadow deficit WITHOUT leakage (i.e. was the residual just an
under-estimated D_eff)? Variants fit on both datasets:
  M3ref        : D free,             no leakage   (reference, masked builder)
  M5_physD     : D pinned = D_phys,  no leakage
  M5_physD_leak: D pinned = D_phys,  + leakage(eta,tau_leak)
  M5_lowD_leak : D pinned = M3 value + leakage    (= HN2/HN8 fixed-D)

Read-only on originals; grid built in the research copy (model_core.compute_grid).
"""
import sys, json
from pathlib import Path
import numpy as np
import lmfit
from scipy import stats

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import fit_lib, t6_hn6, model_core
sys.path.insert(0, str(HERE.parents[1]))
from unified_optimizer import config, model_blanco

W_AMP, W_PHASE = fit_lib.W_AMP, fit_lib.W_PHASE
M3_D = {"356att": 4.683, "test_grid_40_20": 11.447}
PHYS_D = {"356att": 11.0, "test_grid_40_20": 18.0}   # test_grid: micrograph GT D=18 (was nominal 20)


def build_experiment_masked(dataset):
    """fit_lib.build_experiment + dynamic-range noise mask (|exp|>1.5 sqrt(t_noise))."""
    from unified_optimizer.data_manager import DataManager
    data = DataManager(config.DATA_DIR).get_data_for_dataset(dataset)
    angles, freqs, exp, valid = fit_lib.build_experiment(data)
    tnoise = t6_hn6.noise_floor(dataset)
    dr = np.abs(exp) > 1.5*np.sqrt(np.broadcast_to(tnoise[None, :], exp.shape))
    return angles, freqs, exp, (valid & dr)


def residual(params, angles_val, freqs, exp, valid):
    return model_core.residual_from_params(params, angles_val, freqs, exp, valid,
                                           w_amp=W_AMP, w_phase=W_PHASE)


def fit_variant(dataset, D_um, D_vary, leakage):
    angles_val, freqs, exp, valid = build_experiment_masked(dataset)
    p_um = fit_lib.GEOMETRY[dataset]["P_um"]
    P = lmfit.Parameters()
    P.add("P_um", value=p_um, vary=False)
    if D_vary:
        P.add("D_um", value=D_um, min=0.1, max=p_um-0.5)
    else:
        P.add("D_um", value=D_um, vary=False)
    P.add("loss_factor", value=0.3, min=0.0, max=5.0)
    P.add("gamma", value=2.0, vary=False)
    P.add("angle_offset", value=0.0, min=-10, max=10)
    P.add("tau_ps", value=0.0, min=-10, max=10)
    P.add("tau_par_ps", value=0.0, min=-5, max=5)
    P.add("eta0", value=0.02 if leakage else 0.0, min=0.0, max=1.0, vary=leakage)
    P.add("eta_exp", value=1.0, min=-2, max=6, vary=leakage)
    P.add("delta0", value=0.0, min=-np.pi, max=np.pi, vary=leakage)
    P.add("tau_leak", value=0.0, min=-5, max=5, vary=leakage)
    res = lmfit.Minimizer(residual, P, fcn_args=(angles_val, freqs, exp, valid)).minimize(method="leastsq")

    theo = model_core.grid_from_params(res.params, angles_val, freqs)
    ampdb_grid = np.where(valid,
        20*np.log10(np.maximum(np.abs(exp),1e-12))-20*np.log10(np.maximum(np.abs(theo),1e-12)), np.nan)
    hi = angles_val >= np.percentile(angles_val, 75)
    ds_bias = float(np.nanmean(ampdb_grid[hi, :])); ds_rmse = float(np.sqrt(np.nanmean(ampdb_grid[hi, :]**2)))
    amp = (np.abs(exp)-np.abs(theo))[valid]; x = amp[np.isfinite(amp)]
    dw = float(np.sum(np.diff(x)**2)/np.sum(x**2)) if np.sum(x**2) > 0 else None
    rmsedb = float(np.sqrt(np.nanmean(ampdb_grid[np.isfinite(ampdb_grid)]**2)))
    return {"aic": float(res.aic), "bic": float(res.bic), "redchi": float(res.redchi),
            "npts": int(np.sum(valid)), "rmse_amp_db": rmsedb,
            "deepshadow_bias_db": ds_bias, "deepshadow_rmse_db": ds_rmse, "durbin_watson": dw,
            "D_um": float(res.params["D_um"].value), "eta0": float(res.params["eta0"].value),
            "eta_exp": float(res.params["eta_exp"].value), "tau_leak": float(res.params["tau_leak"].value),
            "tau_par_ps": float(res.params["tau_par_ps"].value)}
