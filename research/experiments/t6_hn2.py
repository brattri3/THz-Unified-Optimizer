"""T6 test HN2 — explicit frequency-dependent leakage of the parallel channel.

Physical idea (Karimi 2026, 10.21203/rs.3.rs-8474834/v1): a real WGP has a
finite, frequency-dependent parallel-transmission floor eta(nu) with phase
retardance delta -- the "leaky" coefficient t_par = eta * e^{j delta}. In the
sandwiched Malus model the deep-shadow (theta->90) signal IS the parallel
channel, so an explicit leakage floor should absorb the T2 deep-shadow deficit
that a near-zero Blanco t_par (driven by the low fitted D_eff) cannot.

Model (research COPY of compute_theoretical_grid_2d; original untouched):
    t_par_total = t_par_blanco + eta(nu) * e^{j delta},   eta(nu)=eta0 * nu^a
    E_out = cos^2(th-off) t_perp_eff + sin^2(th-off) t_par_eff
Ablation: M3 (eta0=0 fixed) vs M3+HN2 (eta0,a,delta free), both datasets, same
pipeline. Verdict by dAIC/dBIC + deep-shadow residual + parameter sanity.

Read-only on originals; output under research/results/HN2/.
"""
import sys, json
from pathlib import Path
import numpy as np
import lmfit
from scipy import stats

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import fit_lib
sys.path.insert(0, str(HERE.parents[1]))
from unified_optimizer import config, model_blanco

OUT = HERE.parents[0] / "results" / "HN2"
OUT.mkdir(parents=True, exist_ok=True)

W_AMP, W_PHASE = fit_lib.W_AMP, fit_lib.W_PHASE


def compute_grid_hn2(angles_deg, freqs_thz, p, d, loss_factor, angle_offset, tau_ps,
                     gamma, tau_par_ps, eta0, eta_exp, delta0, N=15, use_drude=True):
    """Faithful copy of compute_theoretical_grid_2d + additive leakage on t_par."""
    d_over_p = d / p
    t_perp_arr, t_par_arr = [], []
    for f in freqs_thz:
        lam = model_blanco.C_LIGHT / (f*1e12)
        pol = p/lam
        t_perp_arr.append(model_blanco.compute_t_perp(pol, d_over_p, N, freq_thz=f, use_drude=use_drude))
        t_par_arr.append(model_blanco.compute_t_par(pol, d_over_p, N, freq_thz=f, use_drude=use_drude))
    t_perp_arr = np.array(t_perp_arr); t_par_arr = np.array(t_par_arr)

    # HN2 leakage floor added to parallel channel
    eta = eta0 * (freqs_thz ** eta_exp)
    t_par_arr = t_par_arr + eta * np.exp(1j*delta0)

    angles_deg = np.array(angles_deg)
    adj = np.deg2rad(angles_deg - angle_offset)
    cos_a = np.cos(adj)[:, None]; sin_a = np.sin(adj)[:, None]
    t_perp = t_perp_arr[None, :]; t_par = t_par_arr[None, :]

    if config.USE_POWER_LAW:
        lf = loss_factor/4.343
        loss_amp = np.exp(-0.5*lf*(freqs_thz**gamma))[None, :]
    else:
        loss_amp = np.exp(-0.5*loss_factor*freqs_thz)[None, :]
    phase_perp = np.exp(-1j*2*np.pi*freqs_thz*tau_ps)[None, :]
    phase_par = np.exp(-1j*2*np.pi*freqs_thz*(tau_ps+tau_par_ps))[None, :]

    E_out = cos_a**2 * (t_perp*loss_amp*phase_perp) + sin_a**2 * (t_par*loss_amp*phase_par)
    return E_out


def residual_hn2(params, angles_val, freqs, exp, valid):
    p = params["P_um"].value*1e-6; d = params["D_um"].value*1e-6
    if d >= p:
        return np.ones(int(np.sum(valid))*2)*1e6
    theo = compute_grid_hn2(
        angles_val, freqs, p, d, params["loss_factor"].value, params["angle_offset"].value,
        params["tau_ps"].value, params["gamma"].value, params["tau_par_ps"].value,
        params["eta0"].value, params["eta_exp"].value, params["delta0"].value)
    em, tm = exp[valid], theo[valid]
    amp = np.abs(em)-np.abs(tm)
    ph = np.angle(em)-np.angle(tm); ph = np.arctan2(np.sin(ph), np.cos(ph))
    return np.concatenate([amp*W_AMP, ph*W_PHASE])


def fit_variant(dataset, leakage, fix_D=None):
    geo = fit_lib.GEOMETRY[dataset]
    from unified_optimizer.data_manager import DataManager
    dm = DataManager(config.DATA_DIR)
    data = dm.get_data_for_dataset(dataset)
    angles_val, freqs, exp, valid = fit_lib.build_experiment(data)

    p_um = geo["P_um"]
    P = lmfit.Parameters()
    P.add("P_um", value=p_um, vary=False)
    d_init = model_blanco.estimate_deff_initial(p_um, geo["D_phys_um"])
    if fix_D is not None:
        P.add("D_um", value=fix_D, vary=False)
    else:
        P.add("D_um", value=d_init, min=0.1, max=p_um-0.5)
    P.add("loss_factor", value=0.3, min=0.0, max=5.0)
    P.add("gamma", value=2.0, min=-5.0, max=10.0, vary=False)   # M3 uses fixed gamma=2
    P.add("angle_offset", value=0.0, min=-10.0, max=10.0)
    P.add("tau_ps", value=0.0, min=-10.0, max=10.0)
    P.add("tau_par_ps", value=0.0, min=-5.0, max=5.0)
    # HN2 leakage params
    P.add("eta0", value=0.02 if leakage else 0.0, min=0.0, max=1.0, vary=leakage)
    P.add("eta_exp", value=1.0, min=-2.0, max=6.0, vary=leakage)
    P.add("delta0", value=0.0, min=-np.pi, max=np.pi, vary=leakage)

    mini = lmfit.Minimizer(residual_hn2, P, fcn_args=(angles_val, freqs, exp, valid))
    res = mini.minimize(method="leastsq")

    # residual decomposition
    theo = compute_grid_hn2(angles_val, freqs, res.params["P_um"].value*1e-6,
        res.params["D_um"].value*1e-6, res.params["loss_factor"].value,
        res.params["angle_offset"].value, res.params["tau_ps"].value,
        res.params["gamma"].value, res.params["tau_par_ps"].value,
        res.params["eta0"].value, res.params["eta_exp"].value, res.params["delta0"].value)
    em, tm = exp[valid], theo[valid]
    amp_res = np.abs(em)-np.abs(tm)
    ampdb = 20*np.log10(np.maximum(np.abs(em),1e-12))-20*np.log10(np.maximum(np.abs(tm),1e-12))
    # deep-shadow (top 25% angles) amp residual in dB, full grid
    ampdb_grid = np.where(valid,
        20*np.log10(np.maximum(np.abs(exp),1e-12))-20*np.log10(np.maximum(np.abs(theo),1e-12)), np.nan)
    hi = angles_val >= np.percentile(angles_val, 75)
    ds_bias = float(np.nanmean(ampdb_grid[hi, :])); ds_rmse = float(np.sqrt(np.nanmean(ampdb_grid[hi, :]**2)))
    x = amp_res[np.isfinite(amp_res)]
    dw = float(np.sum(np.diff(x)**2)/np.sum(x**2)) if np.sum(x**2)>0 else None
    vname = ("M3+HN2" if leakage else "M3") + ("_fixD" if fix_D is not None else "")
    return {
        "dataset": dataset, "variant": vname,
        "success": bool(res.success), "ndata": int(res.ndata), "nvarys": int(res.nvarys),
        "chisqr": float(res.chisqr), "redchi": float(res.redchi),
        "aic": float(res.aic), "bic": float(res.bic),
        "rmse_amp_db": float(np.sqrt(np.mean(ampdb[np.isfinite(ampdb)]**2))),
        "deepshadow_bias_db": ds_bias, "deepshadow_rmse_db": ds_rmse,
        "durbin_watson_amp": dw,
        "params": {k: float(v.value) for k, v in res.params.items()},
    }


def main():
    # M3-fitted D_eff (from T1 baseline) to isolate leakage from the eta0-Deff degeneracy
    M3_D = {"356att": 4.683, "test_grid_40_20": 11.447}
    rows = []
    for ds in fit_lib.GEOMETRY:
        m3 = fit_variant(ds, leakage=False)
        hn2 = fit_variant(ds, leakage=True)
        hn2f = fit_variant(ds, leakage=True, fix_D=M3_D[ds])  # leakage-only (D fixed at M3)
        for r in (m3, hn2, hn2f):
            r["dAIC_vs_M3"] = r["aic"]-m3["aic"]
            r["dBIC_vs_M3"] = r["bic"]-m3["bic"]
        rows.extend([m3, hn2, hn2f])
    (OUT / "hn2_ablation.json").write_text(json.dumps(rows, ensure_ascii=False, indent=2), encoding="utf-8")

    print("=== T6 HN2: explicit parallel-channel leakage ablation ===")
    hdr = f"{'dataset':<16}{'variant':<8}{'AIC':>10}{'dAIC':>9}{'rmse_dB':>9}{'DS_bias':>9}{'DS_rmse':>9}{'DW':>7}"
    print(hdr)
    for r in rows:
        print(f"{r['dataset']:<16}{r['variant']:<8}{r['aic']:>10.1f}{r['dAIC_vs_M3']:>9.1f}"
              f"{r['rmse_amp_db']:>9.3f}{r['deepshadow_bias_db']:>9.3f}{r['deepshadow_rmse_db']:>9.3f}"
              f"{(r['durbin_watson_amp'] or float('nan')):>7.3f}")
    print("\nHN2 fitted leakage params:")
    for r in rows:
        if r["variant"].startswith("M3+HN2"):
            pp = r["params"]
            print(f"  {r['dataset']:<16} eta0={pp['eta0']:.4f} eta_exp={pp['eta_exp']:.3f} "
                  f"delta0={pp['delta0']:+.3f} rad | D_um={pp['D_um']:.3f}")
    print(f"\nartifacts -> {OUT}")


if __name__ == "__main__":
    main()
