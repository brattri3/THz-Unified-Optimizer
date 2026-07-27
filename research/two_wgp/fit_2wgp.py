"""A2 — refit 356att with the TWO-WGP model (reframed after A1).

A1 already falsified hypothesis A analytically (alignment theorem: with WGP2's
transmission axis along the detector axis, d^T J_WGP2 == t_perp2(nu) * d^T, so a
real WGP2 is indistinguishable from an ideal analyser times a scalar).  This
script converts that into a MEASURED statement on the real dataset and tests the
one channel through which WGP2 could still matter — MISALIGNMENT.

Ladder (all variants share the masked experiment builder, the amp/phase weights
and the free-parameter set of M5, so AIC is directly comparable):

  M5_1wgp        production one-WGP M5 grid (compute_grid_hn8)     [baseline]
  W2_shared      two-WGP, WGP2 IDENTICAL to WGP1, phi2=det=E_in=0,
                 loss applied once ("global")  -> isolates the pure WGP2 effect.
                 PREDICTION from the theorem: dAIC ~ 0, D_eff unchanged.
  W2_shared_pe   same but loss applied per element (physical for two grids);
                 differs from W2_shared only by loss_factor renormalisation.
  W2_mis         + free phi2 (misalignment delta = phi2 - det).  At E_in = 0 the
                 theta=90 channel is still proportional to t_par (A1/T9), so this
                 should mostly trade against angle_offset -> watch the correlation.
  W2_mis_psi     + free emitter tilt psi.  Only the PRODUCT sin(delta)*sin(psi)
                 creates a genuine angle-independent floor -> the physical,
                 1-parameter-in-effect replacement of the phenomenological
                 eps_cross.  Expect a strong delta<->psi correlation (ridge).
  M5_eps         one-WGP M5 + angle-independent complex eps_cross (T12 form,
                 2 free params) -> the incumbent phenomenology to beat.

Read-only on originals; artefacts under research/results/two_wgp/.
"""
from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import numpy as np
import lmfit

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(REPO / "research" / "experiments"))
sys.path.insert(0, str(REPO))

import model_2wgp as m2                                    # noqa: E402
import model_m5 as m5                                      # noqa: E402
from t6_hn8 import compute_grid_hn8                        # noqa: E402

OUT = REPO / "research" / "results" / "two_wgp"
OUT.mkdir(parents=True, exist_ok=True)
W_AMP, W_PHASE = m5.W_AMP, m5.W_PHASE
PHYS_D = {"356att": 11.0}


# --------------------------------------------------------------------------
def base_params(p_um, *, phi2=False, psi=False, eps=False, leak=True, seed=None):
    P = lmfit.Parameters()
    P.add("P_um", value=p_um, vary=False)
    P.add("D_um", value=m5.M3_D.get("356att", 4.683), min=0.5, max=p_um - 0.5)
    P.add("loss_factor", value=0.3, min=0.0, max=5.0)
    P.add("gamma", value=2.0, vary=False)          # HN11: degenerate, pinned
    P.add("angle_offset", value=0.0, min=-10, max=10)
    P.add("tau_ps", value=0.0, min=-10, max=10)
    P.add("tau_par_ps", value=0.0, min=-5, max=5)
    P.add("eta0", value=0.02 if leak else 0.0, min=0.0, max=1.0, vary=leak)
    P.add("eta_exp", value=1.0, min=-2, max=6, vary=leak)
    P.add("delta0", value=0.0, min=-np.pi, max=np.pi, vary=leak)
    P.add("tau_leak", value=0.0, min=-5, max=5, vary=leak)
    P.add("phi2_deg", value=0.0, min=-30, max=30, vary=phi2)
    P.add("psi_deg", value=0.0, min=-30, max=30, vary=psi)
    P.add("eps_re", value=0.0, min=-0.2, max=0.2, vary=eps)
    P.add("eps_im", value=0.0, min=-0.2, max=0.2, vary=eps)
    # warm start from another variant's optimum (guards against the local minima
    # that made the cold-started psi fit land ABOVE its own nested submodel)
    if seed:
        for k, v in seed.items():
            if k in P:
                P[k].set(value=min(max(v, P[k].min), P[k].max))
    return P


def model_grid(params, angles, freqs, kind, loss_mode="global"):
    p_um = params["P_um"].value
    d_um = params["D_um"].value
    if kind == "one":
        E = compute_grid_hn8(angles, freqs, p_um * 1e-6, d_um * 1e-6,
                             params["loss_factor"].value, params["angle_offset"].value,
                             params["tau_ps"].value, params["gamma"].value,
                             params["tau_par_ps"].value, params["eta0"].value,
                             params["eta_exp"].value, params["delta0"].value,
                             params["tau_leak"].value)
    else:
        E = m2.compute_grid_2wgp(
            angles, freqs, p_um, d_um,
            params["loss_factor"].value, params["angle_offset"].value,
            params["tau_ps"].value, params["gamma"].value, params["tau_par_ps"].value,
            phi2_deg=params["phi2_deg"].value, e_in_deg=params["psi_deg"].value,
            det_deg=0.0, eta0=params["eta0"].value, eta_exp=params["eta_exp"].value,
            delta0=params["delta0"].value, tau_leak_ps=params["tau_leak"].value,
            loss_mode=loss_mode)
    if params["eps_re"].vary or params["eps_im"].vary:
        E = E + (params["eps_re"].value + 1j * params["eps_im"].value)
    return E


def residual(params, angles, freqs, exp, valid, kind, loss_mode):
    if params["D_um"].value >= params["P_um"].value:
        return np.ones(int(np.sum(valid)) * 2) * 1e6
    theo = model_grid(params, angles, freqs, kind, loss_mode)
    em, tm = exp[valid], theo[valid]
    amp = np.abs(em) - np.abs(tm)
    ph = np.angle(em) - np.angle(tm)
    ph = np.arctan2(np.sin(ph), np.cos(ph))
    return np.concatenate([amp * W_AMP, ph * W_PHASE])


def _correl(res, a, b):
    p = res.params.get(a)
    if p is None or not p.correl or b not in p.correl:
        return None
    return float(p.correl[b])


def run_variant(label, dataset, angles, freqs, exp, valid, *, kind,
                loss_mode="global", phi2=False, psi=False, eps=False,
                leak=True, seed=None):
    p_um = m5.fit_lib.GEOMETRY[dataset]["P_um"]
    P = base_params(p_um, phi2=phi2, psi=psi, eps=eps, leak=leak, seed=seed)
    t0 = time.time()
    res = lmfit.Minimizer(residual, P,
                          fcn_args=(angles, freqs, exp, valid, kind, loss_mode)
                          ).minimize(method="leastsq")
    theo = model_grid(res.params, angles, freqs, kind, loss_mode)
    ampdb = np.where(valid, 20 * np.log10(np.maximum(np.abs(exp), 1e-12))
                     - 20 * np.log10(np.maximum(np.abs(theo), 1e-12)), np.nan)
    hi = angles >= np.percentile(angles, 75)
    amp = (np.abs(exp) - np.abs(theo))[valid]
    x = amp[np.isfinite(amp)]
    dw = float(np.sum(np.diff(x) ** 2) / np.sum(x ** 2)) if np.sum(x ** 2) > 0 else None
    D = float(res.params["D_um"].value)
    out = {
        "label": label, "kind": kind, "loss_mode": loss_mode,
        "nvarys": int(res.nvarys), "npts": int(np.sum(valid)),
        "aic": float(res.aic), "bic": float(res.bic), "redchi": float(res.redchi),
        "rmse_amp_db": float(np.sqrt(np.nanmean(ampdb[np.isfinite(ampdb)] ** 2))),
        "deepshadow_bias_db": float(np.nanmean(ampdb[hi, :])),
        "deepshadow_rmse_db": float(np.sqrt(np.nanmean(ampdb[hi, :] ** 2))),
        "durbin_watson": dw,
        "D_um": D, "D_over_Dphys": D / PHYS_D[dataset],
        "seconds": round(time.time() - t0, 1),
        "params": {k: {"v": float(v.value),
                       "err": (float(v.stderr) if v.stderr else None),
                       "vary": bool(v.vary)} for k, v in res.params.items()},
        "correl": {
            "phi2~angle_offset": _correl(res, "phi2_deg", "angle_offset"),
            "phi2~psi": _correl(res, "phi2_deg", "psi_deg"),
            "D~eta0": _correl(res, "D_um", "eta0"),
            "psi~angle_offset": _correl(res, "psi_deg", "angle_offset"),
        },
    }
    return out


def main():
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("--dataset", default="356att")
    args = ap.parse_args()
    ds = args.dataset

    angles, freqs, exp, valid = m5.build_experiment_masked(ds)
    print(f"=== A2: two-WGP refit of {ds} ===")
    print(f"angles {angles.min():.0f}..{angles.max():.0f} (n={len(angles)}), "
          f"freqs {freqs.min():.2f}..{freqs.max():.2f} THz (n={len(freqs)}), "
          f"valid pts {int(np.sum(valid))}")

    common = dict(dataset=ds, angles=angles, freqs=freqs, exp=exp, valid=valid)
    rows = []
    # --- no-leakage rung: this is where the published D_eff/D_phys=0.426 lives
    rows.append(run_variant("M3ref_1wgp", kind="one", leak=False, **common))
    rows.append(run_variant("W2_shared_noleak", kind="two", leak=False, **common))
    # --- full M5 rung (leakage on) — the production baseline
    rows.append(run_variant("M5_1wgp", kind="one", **common))
    rows.append(run_variant("W2_shared", kind="two", loss_mode="global", **common))
    rows.append(run_variant("W2_shared_pe", kind="two", loss_mode="per_element", **common))
    # --- the only channel through which WGP2 can matter: misalignment
    r_mis = run_variant("W2_mis", kind="two", loss_mode="global", phi2=True, **common)
    rows.append(r_mis)
    seed = {k: v["v"] for k, v in r_mis["params"].items()}
    seed["psi_deg"] = 5.0                       # off zero: only the PRODUCT matters
    rows.append(run_variant("W2_mis_psi", kind="two", loss_mode="global",
                            phi2=True, psi=True, seed=seed, **common))
    # --- incumbent phenomenology, and whether misalignment removes the need for it
    r_eps = run_variant("M5_eps", kind="one", eps=True, **common)
    rows.append(r_eps)
    # Same model, TWO starts: the misalignment basin and the eps_cross basin.
    # The combined model approximately contains M5_eps (phi2=psi=0 leaves only the
    # t_perp2 factor), so a start-dependent AIC gap is evidence of multimodality,
    # and the deeper basin is the data's actual preference.
    rows.append(run_variant("W2_mis_psi_eps", kind="two", loss_mode="global",
                            phi2=True, psi=True, eps=True, seed=seed, **common))
    seed_eps = {k: v["v"] for k, v in r_eps["params"].items()}
    seed_eps.update(phi2_deg=0.0, psi_deg=0.0)
    rows.append(run_variant("W2_mis_psi_eps_seedEps", kind="two", loss_mode="global",
                            phi2=True, psi=True, eps=True, seed=seed_eps, **common))

    by = {r["label"]: r for r in rows}
    base = by["M5_1wgp"]["aic"]
    for r in rows:
        r["dAIC_vs_M5_1wgp"] = r["aic"] - base
    by["W2_shared_noleak"]["dAIC_vs_M3ref"] = (by["W2_shared_noleak"]["aic"]
                                               - by["M3ref_1wgp"]["aic"])
    for lab in ("W2_mis", "W2_mis_psi", "W2_mis_psi_eps", "W2_mis_psi_eps_seedEps"):
        by[lab]["dAIC_vs_W2_shared"] = by[lab]["aic"] - by["W2_shared"]["aic"]
    by["W2_mis_psi"]["dAIC_vs_M5_eps"] = by["W2_mis_psi"]["aic"] - by["M5_eps"]["aic"]
    for lab in ("W2_mis_psi_eps", "W2_mis_psi_eps_seedEps"):
        by[lab]["dAIC_vs_M5_eps"] = by[lab]["aic"] - by["M5_eps"]["aic"]
    by["W2_mis_psi_eps"]["basin_gap_vs_seedEps"] = (by["W2_mis_psi_eps"]["aic"]
                                                    - by["W2_mis_psi_eps_seedEps"]["aic"])
    by["W2_mis_psi_eps"]["dAIC_vs_W2_mis_psi"] = (by["W2_mis_psi_eps"]["aic"]
                                                  - by["W2_mis_psi"]["aic"])

    payload = {"dataset": ds, "D_phys_um": PHYS_D[ds], "variants": rows,
               "cache": m2.cache_stats()}
    (OUT / f"a2_{ds}_ladder.json").write_text(
        json.dumps(payload, ensure_ascii=False, indent=2), encoding="utf-8")

    print(f"\n{'variant':<14}{'nv':>3}{'AIC':>11}{'dAIC':>9}{'D_eff':>8}{'D/Dphys':>9}"
          f"{'DS_bias':>9}{'DS_rmse':>9}{'DW':>7}{'sec':>7}")
    for r in rows:
        print(f"{r['label']:<14}{r['nvarys']:>3}{r['aic']:>11.1f}{r['dAIC_vs_M5_1wgp']:>+9.1f}"
              f"{r['D_um']:>8.3f}{r['D_over_Dphys']:>9.3f}{r['deepshadow_bias_db']:>+9.2f}"
              f"{r['deepshadow_rmse_db']:>9.3f}{r['durbin_watson']:>7.3f}{r['seconds']:>7.1f}")

    def _f(x, fmt="{:+.3f}"):
        return fmt.format(x) if x is not None else "n/a"

    print("\nconfirmatory (theorem T8 predicts ~0 effect):")
    print(f"  W2_shared_noleak vs M3ref_1wgp: dAIC = "
          f"{by['W2_shared_noleak']['dAIC_vs_M3ref']:+.2f}, "
          f"D_eff {by['M3ref_1wgp']['D_um']:.3f} -> {by['W2_shared_noleak']['D_um']:.3f} um "
          f"(D/D_phys {by['M3ref_1wgp']['D_over_Dphys']:.3f} -> "
          f"{by['W2_shared_noleak']['D_over_Dphys']:.3f})")
    print(f"  W2_shared        vs M5_1wgp   : dAIC = "
          f"{by['W2_shared']['dAIC_vs_M5_1wgp']:+.2f}, "
          f"D_eff {by['M5_1wgp']['D_um']:.3f} -> {by['W2_shared']['D_um']:.3f} um")

    print("\nmisalignment channel:")
    for lab in ("W2_mis", "W2_mis_psi", "W2_mis_psi_eps", "W2_mis_psi_eps_seedEps"):
        r = by[lab]
        ph, ps = r["params"]["phi2_deg"], r["params"]["psi_deg"]
        print(f"  {lab:<15} dAIC vs W2_shared = {r['dAIC_vs_W2_shared']:+.1f}   "
              f"phi2 = {_f(ph['v'])} +- {_f(ph['err'], '{:.3f}')}   "
              f"psi = {_f(ps['v'])} +- {_f(ps['err'], '{:.3f}')}")
        print(f"                  correl phi2~angle_offset="
              f"{_f(r['correl']['phi2~angle_offset'])}  "
              f"phi2~psi={_f(r['correl']['phi2~psi'])}")
    for lab in ("M5_eps", "W2_mis_psi_eps", "W2_mis_psi_eps_seedEps"):
        e = by[lab]["params"]
        print(f"  {lab:<15} |eps_cross| = "
              f"{np.hypot(e['eps_re']['v'], e['eps_im']['v']):.5f}   "
              f"AIC = {by[lab]['aic']:.1f}")
    print(f"  W2_mis_psi     vs M5_eps    : dAIC = {by['W2_mis_psi']['dAIC_vs_M5_eps']:+.1f}")
    print(f"  W2_mis_psi_eps vs M5_eps    : dAIC = {by['W2_mis_psi_eps']['dAIC_vs_M5_eps']:+.1f}")
    print(f"  W2_mis_psi_eps vs W2_mis_psi: dAIC = "
          f"{by['W2_mis_psi_eps']['dAIC_vs_W2_mis_psi']:+.1f}  "
          f"(is eps_cross STILL needed once misalignment is free?)")
    print(f"  SAME model, two starts: misalignment basin {by['W2_mis_psi_eps']['aic']:.1f} vs "
          f"eps basin {by['W2_mis_psi_eps_seedEps']['aic']:.1f}  -> gap "
          f"{by['W2_mis_psi_eps']['basin_gap_vs_seedEps']:+.1f}")
    print(f"\nartifact -> {OUT / f'a2_{ds}_ladder.json'}")


if __name__ == "__main__":
    main()
