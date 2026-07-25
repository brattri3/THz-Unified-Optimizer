"""T6 test HN6 — noise model / weighted likelihood + INTEGRITY CHECK of the cluster.

Two goals:
 (1) Naftaly 2008 (10.1364/josab.25.001059): THz-TDS noise is largely
     multiplicative (drift), so an unweighted chi^2 is mis-specified. Test a
     weighted residual (multiplicative: w ~ 1/|E|; additive: w ~ 1/sqrt(noise))
     and check whether residual normality/autocorrelation improves.
 (2) INTEGRITY: fit_lib.build_experiment applies only the water mask, NOT the
     dynamic-range noise mask that optimizer_2d uses (|exp| > 1.5*sqrt(t_noise)).
     HN7 showed the dense deep shadow is only ~3.6x above noise. Re-fit M3 vs
     M3+HN2 WITH the noise mask to check the leakage improvement is not a
     noise-floor artifact.

Dense 356att, D_eff pinned (clean regime). Read-only; output results/HN6/.
"""
import sys, json
from pathlib import Path
import numpy as np
import lmfit
from scipy import stats

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import fit_lib, t6_hn2
sys.path.insert(0, str(HERE.parents[1]))
from unified_optimizer import config, model_blanco

OUT = HERE.parents[0] / "results" / "HN6"
OUT.mkdir(parents=True, exist_ok=True)
M3_D = {"356att": 4.683, "test_grid_40_20": 11.447}


def noise_floor(dataset):
    from unified_optimizer.data_manager import DataManager
    from unified_optimizer.optimizer_2d import get_transmission_spectra
    data = DataManager(config.DATA_DIR).get_data_for_dataset(dataset)
    a_sorted = sorted(data.keys())
    if config.ANGLES_LIMIT_2D:
        lo, hi = config.ANGLES_LIMIT_2D; a_sorted = [a for a in a_sorted if lo <= a <= hi]
    bg = []; fc = None
    for a in a_sorted:
        t_s, E_s, t_b, E_b = data[a]; f, s_s, s_b, tr = get_transmission_spectra(t_s, E_s, t_b, E_b)
        bg.append(s_b); fc = f if fc is None else fc
    bg_avg = np.mean(bg, axis=0)
    tail = np.where(fc >= 3.0)[0]
    if len(tail) == 0: tail = np.where(fc >= 2.5)[0]
    nfa = np.mean(bg_avg[tail]); tnoise = 1.0/((bg_avg/nfa)**2)
    fmax = getattr(config, "F_MAX_2D", config.F_MAX)
    idx = np.where((fc >= config.F_MIN) & (fc <= fmax))[0]
    return tnoise[idx]   # per-frequency transmission-power noise floor (analysis band)


def make_residual(weight_mode, tnoise_band):
    def residual(params, angles_val, freqs, exp, valid):
        p = params["P_um"].value*1e-6; d = params["D_um"].value*1e-6
        if d >= p:
            return np.ones(int(np.sum(valid))*2)*1e6
        theo = t6_hn2.compute_grid_hn2(angles_val, freqs, p, d, params["loss_factor"].value,
            params["angle_offset"].value, params["tau_ps"].value, params["gamma"].value,
            params["tau_par_ps"].value, params["eta0"].value, params["eta_exp"].value,
            params["delta0"].value)
        em, tm = exp[valid], theo[valid]
        amp = np.abs(em)-np.abs(tm)
        ph = np.angle(em)-np.angle(tm); ph = np.arctan2(np.sin(ph), np.cos(ph))
        if weight_mode == "unweighted":
            w = np.ones_like(amp)
        elif weight_mode == "multiplicative":     # relative error ~ const => w ~ 1/|E|
            w = 1.0/np.maximum(np.abs(em), 1e-6)
        elif weight_mode == "additive_noise":     # w ~ 1/sigma, sigma=sqrt(noise power)
            sig = np.sqrt(np.broadcast_to(tnoise_band[None, :], exp.shape)[valid])
            w = 1.0/np.maximum(sig, 1e-9)
        w = w/np.mean(w)                            # normalize so chi^2 scale comparable
        return np.concatenate([amp*w*fit_lib.W_AMP, ph*w*fit_lib.W_PHASE])
    return residual


def build_params(dataset, leakage):
    p_um = fit_lib.GEOMETRY[dataset]["P_um"]
    P = lmfit.Parameters()
    P.add("P_um", value=p_um, vary=False)
    P.add("D_um", value=M3_D[dataset], vary=False)     # pinned (clean regime)
    P.add("loss_factor", value=0.3, min=0.0, max=5.0)
    P.add("gamma", value=2.0, vary=False)
    P.add("angle_offset", value=0.0, min=-10, max=10)
    P.add("tau_ps", value=0.0, min=-10, max=10)
    P.add("tau_par_ps", value=0.0, min=-5, max=5)
    P.add("eta0", value=0.02 if leakage else 0.0, min=0.0, max=1.0, vary=leakage)
    P.add("eta_exp", value=1.0, min=-2, max=6, vary=leakage)
    P.add("delta0", value=0.0, min=-np.pi, max=np.pi, vary=leakage)
    return P


def run(dataset, weight_mode, leakage, apply_noise_mask):
    from unified_optimizer.data_manager import DataManager
    data = DataManager(config.DATA_DIR).get_data_for_dataset(dataset)
    angles_val, freqs, exp, valid = fit_lib.build_experiment(data)
    tnoise = noise_floor(dataset)
    if apply_noise_mask:                                # optimizer_2d dynamic-range cut
        dr = np.abs(exp) > 1.5*np.sqrt(np.broadcast_to(tnoise[None, :], exp.shape))
        valid = valid & dr
    P = build_params(dataset, leakage)
    res = lmfit.Minimizer(make_residual(weight_mode, tnoise), P,
                          fcn_args=(angles_val, freqs, exp, valid)).minimize(method="leastsq")
    # unweighted residual normality (for comparison, always same metric)
    theo = t6_hn2.compute_grid_hn2(angles_val, freqs, res.params["P_um"].value*1e-6,
        res.params["D_um"].value*1e-6, res.params["loss_factor"].value,
        res.params["angle_offset"].value, res.params["tau_ps"].value, res.params["gamma"].value,
        res.params["tau_par_ps"].value, res.params["eta0"].value, res.params["eta_exp"].value,
        res.params["delta0"].value)
    amp = (np.abs(exp)-np.abs(theo))[valid]
    x = amp[np.isfinite(amp)]
    dw = float(np.sum(np.diff(x)**2)/np.sum(x**2)) if np.sum(x**2) > 0 else None
    sh = float(stats.shapiro(x[:5000]).pvalue) if len(x) >= 8 else None
    return {"dataset": dataset, "weight": weight_mode, "leakage": leakage,
            "noise_mask": apply_noise_mask, "npts": int(np.sum(valid)),
            "aic": float(res.aic), "bic": float(res.bic), "redchi": float(res.redchi),
            "rmse_amp": float(np.sqrt(np.mean(x**2))), "durbin_watson": dw, "shapiro_p": sh,
            "eta0": float(res.params["eta0"].value)}


def main():
    ds = "356att"
    rows = []
    # A) weighting effect on residual normality (no noise mask), M3
    for wm in ("unweighted", "multiplicative", "additive_noise"):
        rows.append(run(ds, wm, leakage=False, apply_noise_mask=False))
    # B) integrity: does HN2 leakage survive the dynamic-range noise mask? (unweighted)
    for mask in (False, True):
        m3 = run(ds, "unweighted", leakage=False, apply_noise_mask=mask)
        hn2 = run(ds, "unweighted", leakage=True, apply_noise_mask=mask)
        m3["dAIC_HN2"] = 0.0; hn2["dAIC_HN2"] = hn2["aic"]-m3["aic"]
        rows.extend([m3, hn2])
    (OUT / "hn6_noise.json").write_text(json.dumps(rows, ensure_ascii=False, indent=2), encoding="utf-8")

    print("=== T6 HN6: noise model + cluster integrity check (356att, D pinned) ===")
    print("\nA) Weighting effect on residual normality (M3, no mask):")
    print(f"{'weight':>16}{'redchi':>12}{'DW':>8}{'shapiro_p':>12}")
    for r in rows[:3]:
        print(f"{r['weight']:>16}{r['redchi']:>12.5f}{r['durbin_watson']:>8.3f}{r['shapiro_p']:>12.2e}")
    print("\nB) HN2 leakage survival vs dynamic-range noise mask (unweighted):")
    print(f"{'noise_mask':>12}{'variant':>9}{'npts':>7}{'AIC':>11}{'dAIC_HN2':>10}{'eta0':>9}")
    for r in rows[3:]:
        v = "M3+HN2" if r["leakage"] else "M3"
        print(f"{str(r['noise_mask']):>12}{v:>9}{r['npts']:>7}{r['aic']:>11.1f}"
              f"{r.get('dAIC_HN2',0.0):>10.1f}{r['eta0']:>9.4f}")
    print(f"\nartifacts -> {OUT}")


if __name__ == "__main__":
    main()
