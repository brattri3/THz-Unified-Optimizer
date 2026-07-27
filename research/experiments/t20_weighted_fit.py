"""T15/Phase-3 — WEIGHTED multi-rep fit with per-point errors (Monday-ready).

Ingests MULTI-REP angular data (e.g. Monday's PureWave run), builds per-point
statistics across reps -> proper inverse-variance weighted fit, and reports
error bars on D_eff / eta0 / tau_leak / gamma (lmfit covariance + optional
bootstrap over reps). Adds the sample to the cross-sample D_eff-law table.

WHY: single-rep fits can't weight points or give errors; the leakage cluster
(eta, tau_leak) and gamma are degenerate, so real per-point sigma (from reps) is
what tightens them — especially near the crossed minimum where SNR is worst.

Robust to the pre-measurement state: if the dataset has no real (non-empty) data
yet, it prints a friendly notice and exits 0 (safe to keep in the repo).

Usage:
    .venv/Scripts/python.exe research/experiments/t20_weighted_fit.py --dataset purewave
    ... --nboot 40      # add bootstrap CIs over reps (slower)
    ... --leak/--no-leak, --dphys 10.6
Accesses DataManager.raw_data directly (per-rep), pairing each rep's signal with
the angle's available background (reps may share one bg). Read-only on originals.
"""
import sys, json, argparse
from pathlib import Path
import numpy as np
import lmfit

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import fit_lib, t6_hn6, model_core
sys.path.insert(0, str(HERE.parents[1]))
from unified_optimizer import config, model_blanco
from unified_optimizer.data_manager import DataManager
from unified_optimizer.optimizer_2d import get_transmission_spectra
from unified_optimizer.utils import find_auto_water_mask

OUT = HERE.parents[0] / "results" / "validation"
OUT.mkdir(parents=True, exist_ok=True)
W_AMP, W_PHASE = fit_lib.W_AMP, fit_lib.W_PHASE
# known effective diameters for the cross-sample D_eff law table
KNOWN = {"356att": (4.683, 11.0, 15.5), "test_grid_40_20": (11.447, 18.0, 38.8),
         "specac": (6.240, 14.0, 24.9)}   # (D_eff, D_phys, P)


def collect_reps(dataset):
    """Return {angle: {'trans':[complex spectra per rep], 'freq':f}} using per-rep
    signals paired with the angle's background (reps may share one bg)."""
    dm = DataManager(config.DATA_DIR)
    per_angle = {}
    bg_of_angle = {}
    for (ds, angle, rep), d in dm.raw_data.items():
        if ds != dataset:
            continue
        if 'bg' in d and angle not in bg_of_angle:
            bg_of_angle[angle] = d['bg']
    for (ds, angle, rep), d in dm.raw_data.items():
        if ds != dataset or 'sig' not in d:
            continue
        bg = d.get('bg') or bg_of_angle.get(angle)
        if bg is None:
            continue
        t_s, E_s = d['sig']; t_b, E_b = bg
        if len(E_s) != len(E_b):
            n = min(len(E_s), len(E_b)); t_s, E_s, t_b, E_b = t_s[:n], E_s[:n], t_b[:n], E_b[:n]
        f, _, _, trans = get_transmission_spectra(t_s, E_s, t_b, E_b)
        per_angle.setdefault(angle, {"trans": [], "freq": f})["trans"].append(trans)
    return per_angle


def build_weighted(dataset):
    per = collect_reps(dataset)
    angles = sorted(a for a in per if per[a]["trans"])
    if config.ANGLES_LIMIT_2D is not None:
        lo, hi = config.ANGLES_LIMIT_2D
        angles = [a for a in angles if lo <= a <= hi]
    if not angles:
        return None
    freq = per[angles[0]]["freq"]
    fmax = getattr(config, "F_MAX_2D", config.F_MAX)
    idx = np.where((freq >= config.F_MIN) & (freq <= fmax))[0]
    freqs = freq[idx]
    A = len(angles); F = len(freqs)
    exp_mean = np.zeros((A, F), complex); amp_std = np.zeros((A, F)); nrep = np.zeros(A, int)
    reps_by_angle = []
    for i, a in enumerate(angles):
        stack = np.array([tr[idx] for tr in per[a]["trans"]])   # [reps, F]
        reps_by_angle.append(stack)
        nrep[i] = stack.shape[0]
        exp_mean[i] = stack.mean(axis=0)
        amp_std[i] = np.std(np.abs(stack), axis=0) if stack.shape[0] > 1 else np.nan
    # water + noise-floor + valid masks
    Ymean_db = np.mean(20*np.log10(np.maximum(np.abs(exp_mean), 1e-12)), axis=0)
    water, _ = find_auto_water_mask(freqs, Ymean_db)
    try:
        tnoise = t6_hn6.noise_floor(dataset)
        dr = np.abs(exp_mean) > 1.5*np.sqrt(np.broadcast_to(tnoise[None, :], exp_mean.shape))
    except Exception:
        dr = np.ones_like(exp_mean, bool)
    valid = np.broadcast_to(water[None, :], exp_mean.shape) & dr
    # amplitude sigma floor: median of finite stds (avoid div0 / single-rep angles)
    finite = amp_std[np.isfinite(amp_std) & (amp_std > 0)]
    sfloor = float(np.median(finite)) if finite.size else 1e-3
    sig = np.where(np.isfinite(amp_std) & (amp_std > 0), amp_std, sfloor)
    sig = np.maximum(sig, 0.05*sfloor)
    return {"angles": np.array(angles, float), "freqs": freqs, "exp": exp_mean,
            "sig": sig, "valid": valid, "nrep": nrep, "reps": reps_by_angle,
            "sigma_floor": sfloor}


def make_resid(G):
    angles, freqs, exp, sig, valid = G["angles"], G["freqs"], G["exp"], G["sig"], G["valid"]
    w = (1.0/sig); w = w/np.mean(w[valid])
    def resid(params):
        return model_core.residual_from_params(params, angles, freqs, exp, valid,
                                               w=w, w_amp=W_AMP, w_phase=W_PHASE)
    return resid


def params0(P_um, D_um, D_vary, leak):
    P = lmfit.Parameters()
    P.add("P_um", value=P_um, vary=False)
    P.add("D_um", value=D_um, min=0.1, max=P_um-0.5, vary=D_vary)
    P.add("loss_factor", value=0.3, min=0.0, max=5.0)
    P.add("gamma", value=2.0, min=0.2, max=5.0)
    P.add("angle_offset", value=0.0, min=-15, max=15)
    P.add("tau_ps", value=0.0, min=-10, max=10)
    P.add("tau_par_ps", value=0.0, min=-5, max=5)
    P.add("eta0", value=0.02 if leak else 0.0, min=0.0, max=1.0, vary=leak)
    P.add("eta_exp", value=1.0, min=-2, max=6, vary=leak)
    P.add("delta0", value=0.0, min=-np.pi, max=np.pi, vary=leak)
    P.add("tau_leak", value=0.1, min=-5, max=5, vary=leak)
    return P


def fit(G, P_um, D_um, D_vary, leak):
    P = params0(P_um, D_um, D_vary, leak)
    res = lmfit.Minimizer(lambda pp: make_resid(G)(pp), P).minimize(method="leastsq")
    return res


def bootstrap(G, P_um, dataset, nboot, seed=1):
    rng = np.random.default_rng(seed)
    keys = ["D_um", "eta0", "eta_exp", "tau_leak", "gamma", "angle_offset"]
    samples = {k: [] for k in keys}
    for _ in range(nboot):
        exp_b = G["exp"].copy()
        for i, stack in enumerate(G["reps"]):
            nr = stack.shape[0]
            if nr > 1:
                sel = rng.integers(0, nr, nr)
                exp_b[i] = stack[sel].mean(axis=0)
        Gb = dict(G); Gb["exp"] = exp_b
        r = lmfit.Minimizer(lambda pp: make_resid(Gb)(pp),
                            params0(P_um, D_um=G["D_eff_pt"], D_vary=True, leak=True)).minimize(method="leastsq")
        for k in keys:
            samples[k].append(float(r.params[k].value))
    ci = {k: [float(np.percentile(v, 16)), float(np.percentile(v, 50)), float(np.percentile(v, 84))]
          for k, v in samples.items()}
    return ci


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dataset", default="purewave")
    ap.add_argument("--dphys", type=float, default=None, help="physical D (um) for D_eff/D_phys")
    ap.add_argument("--nboot", type=int, default=0, help=">0 to add bootstrap CIs over reps")
    a = ap.parse_args()
    ds = a.dataset

    G = build_weighted(ds)
    if G is None:
        print(f"[{ds}] no usable multi-rep data yet (empty placeholders?). "
              f"Fill data_pool/{ds}/ per MEASUREMENT_ORDER.md, then re-run. Nothing to do.")
        return
    P_um = fit_lib.GEOMETRY.get(ds, {}).get("P_um")
    D_phys = a.dphys or fit_lib.GEOMETRY.get(ds, {}).get("D_phys_um")
    print(f"=== weighted multi-rep fit: {ds} (P={P_um}, D_phys={D_phys}) ===")
    print(f"angles={len(G['angles'])} reps/angle={G['nrep'].tolist()} valid_pts={int(np.sum(G['valid']))} "
          f"sigma_floor={G['sigma_floor']:.3e}")

    # M3-style (D free, no leakage) -> effective D; then weighted M5 (D free + leakage)
    m3 = fit(G, P_um, D_um=min(D_phys or P_um*0.6, P_um-1.0), D_vary=True, leak=False)
    Deff = float(m3.params["D_um"].value); G["D_eff_pt"] = Deff
    m5 = fit(G, P_um, D_um=Deff, D_vary=True, leak=True)

    def val(res, k):
        v = res.params[k]; return (float(v.value), (float(v.stderr) if v.stderr else None))
    out = {"dataset": ds, "P_um": P_um, "D_phys_um": D_phys,
           "n_angles": len(G["angles"]), "reps_per_angle": G["nrep"].tolist(),
           "valid_points": int(np.sum(G["valid"])),
           "M3_free": {k: val(m3, k) for k in ("D_um", "gamma", "loss_factor")},
           "M5_weighted": {k: val(m5, k) for k in ("D_um", "eta0", "eta_exp", "tau_leak", "gamma", "angle_offset")},
           "M5_aic": float(m5.aic), "M3_aic": float(m3.aic)}
    if D_phys:
        out["Deff_over_Dphys"] = float(m5.params["D_um"].value)/D_phys

    print("\nM5 weighted fit (value ± covariance-stderr):")
    for k in ("D_um", "eta0", "eta_exp", "tau_leak", "gamma", "angle_offset"):
        v, e = val(m5, k); print(f"  {k:<12} = {v:+.4f} ± {e if e else 'NA'}")
    if D_phys:
        print(f"  D_eff/D_phys = {out['Deff_over_Dphys']:.3f}  (D/P_phys={D_phys/P_um:.2f})")

    if a.nboot > 0:
        print(f"\nbootstrap over reps (nboot={a.nboot}, slower)...")
        out["bootstrap_ci_16_50_84"] = bootstrap(G, P_um, ds, a.nboot)
        for k, ci in out["bootstrap_ci_16_50_84"].items():
            print(f"  {k:<12} CI = [{ci[0]:+.4f}, {ci[1]:+.4f}, {ci[2]:+.4f}]")

    # cross-sample D_eff law table (append this sample)
    tbl = {}
    for s, (de, dp, pp) in KNOWN.items():
        tbl[f"{s}(D/P{dp/pp:.2f})"] = round(de/dp, 3)
    if D_phys:
        tbl[f"{ds}(D/P{D_phys/P_um:.2f})"] = round(float(m5.params["D_um"].value)/D_phys, 3)
    out["cross_sample_Deff_over_Dphys"] = tbl
    print("\nD_eff/D_phys across samples (law check):")
    for k, v in sorted(tbl.items()):
        print(f"  {k:<22} {v}")

    (OUT / f"{ds}_weighted_fit.json").write_text(json.dumps(out, ensure_ascii=False, indent=2), encoding="utf-8")
    print(f"\nartifacts -> {OUT/ (ds+'_weighted_fit.json')}")


if __name__ == "__main__":
    main()
