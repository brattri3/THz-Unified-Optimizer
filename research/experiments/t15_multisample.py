"""T15 (step 1) — cross-sample generalization of the leakage floor (model-independent).

Runs the HN7 Parseval energy diagnostic U(theta) on ALL data_pool datasets
(356att, test_grid_40_20, series1..5) to test whether the deep-shadow leakage
floor eta = sqrt(U(90)/U(0)) generalizes beyond the 2 analyzed samples, and how
it sits vs the measurement noise floor. NO geometry needed (energy is
model-independent). Read-only; output results/validation/.
"""
import sys, json
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import fit_lib, t6_hn6
sys.path.insert(0, str(HERE.parents[1]))
from unified_optimizer import config
from unified_optimizer.data_manager import DataManager

OUT = HERE.parents[0] / "results" / "validation"
OUT.mkdir(parents=True, exist_ok=True)
DATASETS = ["att-11-16-356", "test_grid_40_20", "att-11-16-s1", "att-11-16-s2", "att-11-16-s3", "att-11-16-s4-artifact", "att-11-16-s5-artifact"]


def leak_floor(dataset):
    dm = DataManager(config.DATA_DIR)
    data = dm.get_data_for_dataset(dataset)
    if not data:
        return {"dataset": dataset, "error": "no data"}
    angles, freqs, exp, valid = fit_lib.build_experiment(data)
    U = np.array([np.sum(np.abs(exp[i, valid[i]])**2) for i in range(len(angles))])
    i0 = int(np.argmin(angles)); i90 = int(np.argmax(angles))
    U0 = U[i0]; U90 = U[i90]
    floor = float(U90/U0); eta = float(np.sqrt(max(floor, 0.0)))
    try:
        tnoise = t6_hn6.noise_floor(dataset)
        Unoise = float(np.sum(tnoise[valid[i0]]))
        snr = U90/Unoise if Unoise > 0 else None
    except Exception:
        snr = None
    return {"dataset": dataset, "n_angles": int(len(angles)),
            "angle_max": float(angles[i90]), "U90_over_U0": floor, "eta_direct": eta,
            "deepshadow_SNR": (float(snr) if snr else None)}


def main():
    rows = [leak_floor(ds) for ds in DATASETS]
    (OUT / "multisample_leakage.json").write_text(json.dumps(rows, ensure_ascii=False, indent=2), encoding="utf-8")
    print("=== T15: cross-sample leakage floor (model-independent, HN7 method) ===")
    print(f"{'dataset':<18}{'n_ang':>6}{'ang_max':>8}{'U90/U0':>12}{'eta':>9}{'SNR':>8}")
    for r in rows:
        if "error" in r:
            print(f"{r['dataset']:<18} {r['error']}"); continue
        snr = r['deepshadow_SNR']
        print(f"{r['dataset']:<18}{r['n_angles']:>6}{r['angle_max']:>8.0f}{r['U90_over_U0']:>12.3e}"
              f"{r['eta_direct']:>9.4f}{(snr if snr else float('nan')):>8.1f}")
    print(f"\nartifacts -> {OUT}")


if __name__ == "__main__":
    main()
