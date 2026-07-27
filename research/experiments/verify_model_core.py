"""Regression guard for model_core: bit-exactness vs the pre-refactor code paths.

Checks, on the REAL analysis frequency grids of both datasets:
  1. model_core.blanco_t          == scalar loop over model_blanco.compute_t_*
  2. model_core.compute_grid      == optimizer_2d.compute_theoretical_grid_2d   (eta0=0)
  3. model_core.compute_grid      == t6_hn2.compute_grid_hn2                    (tau_leak=0)
  4. model_core.compute_grid      == t6_hn8.compute_grid_hn8                    (general)
  5. model_core.residual_from_params == t6_hn8.residual  (same vector)
and reports the speed-up. Tolerance is 0 (exact float equality); the vectorised
core performs the same operations in the same order as the scalar originals.

Run:  .venv/Scripts/python.exe research/experiments/verify_model_core.py
Exit code 1 if anything differs.
"""
import sys, time, warnings
from pathlib import Path

import numpy as np

warnings.filterwarnings("ignore", category=UserWarning)
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parents[1]))

import model_core, fit_lib                                       # noqa: E402
from unified_optimizer import config, model_blanco               # noqa: E402
from unified_optimizer.data_manager import DataManager           # noqa: E402
from unified_optimizer.optimizer_2d import compute_theoretical_grid_2d  # noqa: E402

DATASETS = ("356att", "test_grid_40_20")
FAILURES = []


def check(name, a, b):
    a = np.asarray(a); b = np.asarray(b)
    if a.shape != b.shape:
        FAILURES.append(f"{name}: shape {a.shape} != {b.shape}")
        print(f"  FAIL {name}: shape mismatch")
        return
    exact = bool(np.array_equal(a, b))
    denom = np.maximum(np.abs(b), 1e-300)
    rel = float(np.max(np.abs(a - b) / denom)) if a.size else 0.0
    tag = "EXACT" if exact else f"rel_max={rel:.3e}"
    print(f"  {'ok  ' if exact else 'DIFF'} {name:<46} {tag}")
    if not exact:
        FAILURES.append(f"{name}: not bit-exact (rel_max={rel:.3e})")


def legacy_grid_hn8(angles_deg, freqs_thz, p, d, loss_factor, angle_offset, tau_ps,
                    gamma, tau_par_ps, eta0, eta_exp, delta0, tau_leak,
                    N=15, use_drude=True):
    """FROZEN pre-refactor copy of t6_hn8.compute_grid_hn8 (verbatim).

    Kept here so the equivalence stays testable after t6_hn2/t6_hn8 became thin
    aliases over model_core. With tau_leak=0 this is also the old
    t6_hn2.compute_grid_hn2, and with eta0=0 it is compute_theoretical_grid_2d.
    """
    d_over_p = d/p
    t_perp_arr, t_par_arr = [], []
    for f in freqs_thz:
        lam = model_blanco.C_LIGHT/(f*1e12); pol = p/lam
        t_perp_arr.append(model_blanco.compute_t_perp(pol, d_over_p, N, freq_thz=f, use_drude=use_drude))
        t_par_arr.append(model_blanco.compute_t_par(pol, d_over_p, N, freq_thz=f, use_drude=use_drude))
    t_perp_arr = np.array(t_perp_arr); t_par_arr = np.array(t_par_arr)
    eta = eta0*(freqs_thz**eta_exp)
    delta = delta0 + 2*np.pi*freqs_thz*tau_leak
    t_par_arr = t_par_arr + eta*np.exp(1j*delta)

    adj = np.deg2rad(np.array(angles_deg)-angle_offset)
    cos_a = np.cos(adj)[:, None]; sin_a = np.sin(adj)[:, None]
    if config.USE_POWER_LAW:
        loss_amp = np.exp(-0.5*(loss_factor/4.343)*(freqs_thz**gamma))[None, :]
    else:
        loss_amp = np.exp(-0.5*loss_factor*freqs_thz)[None, :]
    phase_perp = np.exp(-1j*2*np.pi*freqs_thz*tau_ps)[None, :]
    phase_par = np.exp(-1j*2*np.pi*freqs_thz*(tau_ps+tau_par_ps))[None, :]
    return (cos_a**2*(t_perp_arr[None, :]*loss_amp*phase_perp)
            + sin_a**2*(t_par_arr[None, :]*loss_amp*phase_par))


def scalar_blanco(freqs, p, d, N=15, use_drude=True):
    d_over_p = d / p
    tp, ta = [], []
    for f in freqs:
        lam = model_blanco.C_LIGHT / (f * 1e12)
        pol = p / lam
        tp.append(model_blanco.compute_t_perp(pol, d_over_p, N, freq_thz=f, use_drude=use_drude))
        ta.append(model_blanco.compute_t_par(pol, d_over_p, N, freq_thz=f, use_drude=use_drude))
    return np.array(tp), np.array(ta)


def main():
    import t6_hn2, t6_hn8
    import lmfit

    for ds in DATASETS:
        geo = fit_lib.GEOMETRY[ds]
        data = DataManager(config.DATA_DIR).get_data_for_dataset(ds)
        angles, freqs, exp, valid = fit_lib.build_experiment(data)
        p = geo["P_um"] * 1e-6
        print(f"\n=== {ds}: {len(angles)} angles x {len(freqs)} freqs, "
              f"{int(valid.sum())} valid pts ===")

        for d_um in (geo["D_phys_um"], 4.683, 0.5 * geo["P_um"]):
            d = d_um * 1e-6
            if d >= p:
                continue
            for use_drude in (True, False):
                model_core.clear_cache()
                ref = scalar_blanco(freqs, p, d, use_drude=use_drude)
                got = model_core.blanco_t(freqs, p, d, use_drude=use_drude)
                tag = f"blanco_t D={d_um:g}um drude={int(use_drude)}"
                check(tag + " t_perp", got[0], ref[0])
                check(tag + " t_par ", got[1], ref[1])

        # --- grid equivalence, several parameter points -------------------
        d = 4.683e-6 if ds == "356att" else 11.447e-6
        pts = [
            dict(loss_factor=0.3, angle_offset=0.0, tau_ps=0.0, gamma=2.0,
                 tau_par_ps=0.0),
            dict(loss_factor=1.7, angle_offset=-2.3, tau_ps=0.41, gamma=1.85,
                 tau_par_ps=-0.12),
        ]
        for i, q in enumerate(pts):
            ref = compute_theoretical_grid_2d(
                angles, freqs, p, d, q["loss_factor"], q["angle_offset"],
                q["tau_ps"], gamma=q["gamma"], tau_par_ps=q["tau_par_ps"])
            got = model_core.compute_grid(
                angles, freqs, p, d, q["loss_factor"], q["angle_offset"],
                q["tau_ps"], gamma=q["gamma"], tau_par_ps=q["tau_par_ps"])
            check(f"compute_grid[{i}] vs compute_theoretical_grid_2d", got, ref)

            # leakage variants against the FROZEN legacy grid, incl. eta0=0
            # (model_core short-circuits the leakage term when eta0 == 0)
            for lk, tag in (((0.021, 1.3, 0.7, 0.0), "legacy hn2 (tau_leak=0)"),
                            ((0.021, 1.3, 0.7, 0.19), "legacy hn8"),
                            ((0.0, 1.3, 0.7, 0.19), "legacy hn8 eta0=0")):
                ref_l = legacy_grid_hn8(
                    angles, freqs, p, d, q["loss_factor"], q["angle_offset"],
                    q["tau_ps"], q["gamma"], q["tau_par_ps"], *lk)
                got_l = model_core.compute_grid(
                    angles, freqs, p, d, q["loss_factor"], q["angle_offset"],
                    q["tau_ps"], gamma=q["gamma"], tau_par_ps=q["tau_par_ps"],
                    eta0=lk[0], eta_exp=lk[1], delta0=lk[2], tau_leak=lk[3])
                check(f"compute_grid[{i}] vs {tag}", got_l, ref_l)

            # the live aliases must still agree with the core
            check(f"t6_hn2.compute_grid_hn2 alias[{i}]",
                  t6_hn2.compute_grid_hn2(
                      angles, freqs, p, d, q["loss_factor"], q["angle_offset"],
                      q["tau_ps"], q["gamma"], q["tau_par_ps"], 0.021, 1.3, 0.7),
                  legacy_grid_hn8(
                      angles, freqs, p, d, q["loss_factor"], q["angle_offset"],
                      q["tau_ps"], q["gamma"], q["tau_par_ps"], 0.021, 1.3, 0.7, 0.0))
            check(f"t6_hn8.compute_grid_hn8 alias[{i}]",
                  t6_hn8.compute_grid_hn8(
                      angles, freqs, p, d, q["loss_factor"], q["angle_offset"],
                      q["tau_ps"], q["gamma"], q["tau_par_ps"], 0.021, 1.3, 0.7, 0.19),
                  legacy_grid_hn8(
                      angles, freqs, p, d, q["loss_factor"], q["angle_offset"],
                      q["tau_ps"], q["gamma"], q["tau_par_ps"], 0.021, 1.3, 0.7, 0.19))

        # --- residual equivalence -----------------------------------------
        P = lmfit.Parameters()
        P.add("P_um", value=geo["P_um"], vary=False)
        P.add("D_um", value=d * 1e6, vary=False)
        P.add("loss_factor", value=1.7)
        P.add("gamma", value=1.85, vary=False)
        P.add("angle_offset", value=-2.3)
        P.add("tau_ps", value=0.41)
        P.add("tau_par_ps", value=-0.12)
        P.add("eta0", value=0.021)
        P.add("eta_exp", value=1.3)
        P.add("delta0", value=0.7)
        P.add("tau_leak", value=0.19)
        theo_l = legacy_grid_hn8(angles, freqs, p, d, 1.7, -2.3, 0.41, 1.85,
                                 -0.12, 0.021, 1.3, 0.7, 0.19)
        em, tm = exp[valid], theo_l[valid]
        amp = np.abs(em) - np.abs(tm)
        ph = np.angle(em) - np.angle(tm)
        ph = np.arctan2(np.sin(ph), np.cos(ph))
        ref_r = np.concatenate([amp * fit_lib.W_AMP, ph * fit_lib.W_PHASE])
        got_r = model_core.residual_from_params(P, angles, freqs, exp, valid)
        check("residual_from_params vs legacy inline residual", got_r, ref_r)
        check("t6_hn8.residual (alias) vs legacy inline residual",
              t6_hn8.residual(P, angles, freqs, exp, valid), ref_r)

        # fit_lib.residual (no leakage; historical gamma fallback = 2.0)
        Pf = lmfit.Parameters()
        Pf.add("P_um", value=geo["P_um"], vary=False)
        Pf.add("D_um", value=d * 1e6, vary=False)
        Pf.add("loss_factor", value=1.7)
        Pf.add("angle_offset", value=-2.3)
        Pf.add("tau_ps", value=0.41)
        Pf.add("tau_par_ps", value=-0.12)
        for use_drude in (True, False):
            for gam, has_gamma in ((1.85, True), (2.0, False)):   # 2.0 = fallback
                if has_gamma:
                    Pf.add("gamma", value=gam, vary=False)
                elif "gamma" in Pf:
                    del Pf["gamma"]
                theo_f = legacy_grid_hn8(angles, freqs, p, d, 1.7, -2.3, 0.41, gam,
                                         -0.12, 0.0, 1.0, 0.0, 0.0, use_drude=use_drude)
                em, tm = exp[valid], theo_f[valid]
                amp = np.abs(em) - np.abs(tm)
                ph = np.angle(em) - np.angle(tm)
                ph = np.arctan2(np.sin(ph), np.cos(ph))
                ref_f = np.concatenate([amp * fit_lib.W_AMP, ph * fit_lib.W_PHASE])
                check(f"fit_lib.residual drude={int(use_drude)} gamma_param={has_gamma}",
                      fit_lib.residual(Pf, angles, freqs, exp, valid, use_drude), ref_f)
        Pf.add("gamma", value=1.85, vary=False)

        # weighted variant (t20 style): w normalised to mean 1 on `valid`
        rng = np.random.default_rng(0)
        sig = 1.0 + 0.1 * rng.random(exp.shape)
        w = 1.0 / sig
        w = w / np.mean(w[valid])
        theo = model_core.grid_from_params(P, angles, freqs)
        em, tm, ww = exp[valid], theo[valid], w[valid]
        amp = (np.abs(em) - np.abs(tm)) * ww
        ph = np.angle(em) - np.angle(tm)
        ph = np.arctan2(np.sin(ph), np.cos(ph)) * ww
        ref_w = np.concatenate([amp * fit_lib.W_AMP, ph * fit_lib.W_PHASE])
        got_w = model_core.residual_from_params(P, angles, freqs, exp, valid, w=w)
        check("weighted residual vs inline t20 form", got_w, ref_w)

        # --- speed ---------------------------------------------------------
        model_core.clear_cache()
        n = 20
        t0 = time.perf_counter()
        for k in range(n):
            compute_theoretical_grid_2d(angles, freqs, p, d * (1 + 1e-6 * k),
                                        0.3, 0.0, 0.0, gamma=2.0)
        t_old = (time.perf_counter() - t0) / n
        t0 = time.perf_counter()
        for k in range(n):
            model_core.compute_grid(angles, freqs, p, d * (1 + 1e-6 * k),
                                    0.3, 0.0, 0.0, gamma=2.0)
        t_new_miss = (time.perf_counter() - t0) / n
        t0 = time.perf_counter()
        for k in range(n):
            model_core.compute_grid(angles, freqs, p, d, 0.3 + 0.01 * k,
                                    0.0, 0.0, gamma=2.0)
        t_new_hit = (time.perf_counter() - t0) / n
        print(f"  grid eval: old {t_old*1e3:7.2f} ms | new(D changes) "
              f"{t_new_miss*1e3:7.2f} ms ({t_old/t_new_miss:5.1f}x) | "
              f"new(D cached) {t_new_hit*1e3:7.2f} ms ({t_old/t_new_hit:5.1f}x)")
        print(f"  cache: {model_core.cache_stats()}")

    print()
    if FAILURES:
        print(f"FAILED ({len(FAILURES)}):")
        for f in FAILURES:
            print("  -", f)
        return 1
    print("ALL CHECKS EXACT")
    return 0


if __name__ == "__main__":
    sys.exit(main())
