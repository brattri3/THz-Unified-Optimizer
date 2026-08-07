"""T6 test HN7 — Parseval energy U(theta) vs cos^4, INDEPENDENT leakage cross-check.

Kaltenecker 2019 (10.1007/s10762-019-00608-x): an ideal WGP attenuator gives
U_out = U_in cos^4(theta). A real WGP with nonzero parallel transmission t_par
deviates -- a deep-shadow floor. From E_out = t_perp cos^2 + t_par sin^2,
    U(theta) = A cos^4 + B sin^4 + C cos^2 sin^2,
with A=<|t_perp|^2>, B=<|t_par|^2>, C=2 Re<t_perp t_par^*>. The ratio B/A is the
deep-shadow floor = |t_par/t_perp|^2 = 1/ER = leakage^2 -- measured directly from
the angular ENERGY profile, with NO spectral-coefficient fit. This independently
cross-checks HN2 (which fitted eta ~ 0.1 for the dense set => B/A ~ 0.01).

U(theta) uses Parseval: U = sum_{valid nu} |exp(theta,nu)|^2 (band-integrated,
water lines masked). Read-only on originals; output under research/results/HN7/.
"""
import sys, json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import fit_lib
sys.path.insert(0, str(HERE.parents[1]))
from unified_optimizer import config

OUT = HERE.parents[0] / "results" / "HN7"
OUT.mkdir(parents=True, exist_ok=True)


def energy_profile(dataset):
    from unified_optimizer.data_manager import DataManager
    from unified_optimizer.optimizer_2d import get_transmission_spectra
    dm = DataManager(config.DATA_DIR)
    data = dm.get_data_for_dataset(dataset)
    angles, freqs, exp, valid = fit_lib.build_experiment(data)
    # Parseval energy per angle over valid (non-water) frequencies
    U = np.array([np.sum(np.abs(exp[i, valid[i]])**2) for i in range(len(angles))])
    # measurement energy noise floor (replicate optimizer_2d dynamic-range floor)
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
    U_noise = float(np.sum(tnoise[idx][valid[0]]))
    return np.asarray(angles, float), U, U_noise


def fit_models(angles_deg, U, U_noise):
    th = np.deg2rad(angles_deg)
    c2, s2 = np.cos(th)**2, np.sin(th)**2
    # normalize by theta~0 energy for interpretability
    U0 = U[np.argmin(angles_deg)]
    Un = U / U0
    # Model A: pure cos^4 (1 param scale)
    basis1 = (c2**2)[:, None]
    coef1, *_ = np.linalg.lstsq(basis1, Un, rcond=None)
    pred1 = basis1 @ coef1
    rss1 = float(np.sum((Un-pred1)**2))
    # Model B: leaky A cos^4 + B sin^4 + C cos^2 sin^2 (3 params)
    basis3 = np.stack([c2**2, s2**2, c2*s2], axis=1)
    coef3, *_ = np.linalg.lstsq(basis3, Un, rcond=None)
    pred3 = basis3 @ coef3
    rss3 = float(np.sum((Un-pred3)**2))
    A, B, Cc = [float(x) for x in coef3]
    tss = float(np.sum((Un-np.mean(Un))**2))
    # model-independent direct floor at the largest angle (cleaner than lstsq B)
    imax = int(np.argmax(angles_deg))
    floor_direct = float(Un[imax])
    return {
        "U0_raw": float(U0),
        "floor_direct_U90_over_U0": floor_direct,
        "eta_direct": float(np.sqrt(max(floor_direct, 0.0))),
        "cos4_only": {"scale": float(coef1[0]), "rss": rss1, "R2": 1-rss1/tss},
        "leaky": {"A": A, "B": B, "C": Cc, "rss": rss3, "R2": 1-rss3/tss,
                  "floor_BA": B/A if A != 0 else None,
                  "leakage_eta": float(np.sqrt(max(B/A, 0.0))) if A != 0 else None},
        "F_improve_rss_ratio": rss1/rss3 if rss3 > 0 else None,
        "noise_floor_U_over_U0": U_noise/U0, "floor_over_noise_SNR": (floor_direct*U0)/U_noise if U_noise > 0 else None,
        "Un": Un.tolist(), "angles_deg": angles_deg.tolist(),
        "pred_cos4": pred1.tolist(), "pred_leaky": pred3.tolist(),
    }


def main():
    results = {}
    fig, ax = plt.subplots(figsize=(7.5, 5))
    for ds, color in (("att-11-16-356", "tab:red"), ("test_grid_40_20", "tab:blue")):
        ang, U, U_noise = energy_profile(ds)
        r = fit_models(ang, U, U_noise)
        results[ds] = r
        ax.plot(ang, r["Un"], "o", color=color, label=f"{ds} data")
        thf = np.linspace(0, 90, 200); c2f = np.cos(np.deg2rad(thf))**2
        ax.plot(thf, r["cos4_only"]["scale"]*c2f**2, "--", color=color, alpha=0.6,
                label=f"{ds} cos^4 (R2={r['cos4_only']['R2']:.3f})")
        ax.plot(thf, r["leaky"]["A"]*c2f**2 + r["leaky"]["B"]*(1-c2f)**2
                + r["leaky"]["C"]*c2f*(1-c2f), "-", color=color,
                label=f"{ds} leaky (floor B/A={r['leaky']['floor_BA']:.4f})")
    ax.set_yscale("log"); ax.set_xlabel("WGP angle theta (deg)")
    ax.set_ylabel("U(theta)/U(0)  (Parseval energy)")
    ax.set_title("HN7: energy transmission vs cos^4 (deep-shadow floor = leakage^2)")
    ax.legend(fontsize=7); ax.grid(alpha=0.3, which="both"); fig.tight_layout()
    fig.savefig(OUT / "U_theta_vs_cos4.png", dpi=120); plt.close(fig)

    (OUT / "hn7_energy.json").write_text(json.dumps(results, ensure_ascii=False, indent=2), encoding="utf-8")

    print("=== T6 HN7: U(theta) vs cos^4 (independent leakage) ===")
    for ds, r in results.items():
        print(f"\n{ds}:")
        print(f"  cos^4-only : R2={r['cos4_only']['R2']:.4f}  rss={r['cos4_only']['rss']:.3e}")
        print(f"  leaky 3par : R2={r['leaky']['R2']:.4f}  rss={r['leaky']['rss']:.3e}")
        print(f"  DIRECT floor U(90)/U(0) = {r['floor_direct_U90_over_U0']:.3e}  "
              f"=> leakage eta = {r['eta_direct']:.4f}")
        print(f"  noise floor U/U0 = {r['noise_floor_U_over_U0']:.3e}  "
              f"=> deep-shadow SNR = {r['floor_over_noise_SNR']:.1f}x above noise")
        print(f"  (lstsq B/A={r['leaky']['floor_BA']:.5f} -- biased by mid-angles, use DIRECT)")
    print(f"\nartifacts -> {OUT}")


if __name__ == "__main__":
    main()
