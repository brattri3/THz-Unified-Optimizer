"""HN12 (owner micrograph-seeded) — multi-plane depth Delta_z -> delayed leakage.

Micrograph of the P40 grid (test_grid_40_20) shows wires at two focal depths
Delta_z ~ 200 um along the beam axis. Prediction: the crossed (theta~90) residual
leakage traverses the DEEPER wire plane, arriving Delta_z/c ~ 0.67 ps (single) or
2 Delta_z/c ~ 1.33 ps (inter-plane etalon) LATER than the co-pol (theta~0) signal;
a DETECTOR cross-pol would arrive with ZERO delay. Model-independent test:
cross-correlate the raw time-domain sample waveform E_s(t) at theta~90 vs theta~0,
read the lag; also look for etalon ripple ~0.75-1.5 THz in the crossed spectrum.
This directly separates SAMPLE leakage (delayed) from DETECTOR cross-pol (T12).
Read-only; output results/HN12/.
"""
import sys, json
from pathlib import Path
import numpy as np
from scipy.signal import correlate, correlation_lags

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parents[1]))
from unified_optimizer import config
from unified_optimizer.data_manager import DataManager

OUT = HERE.parents[0] / "results" / "HN12"
OUT.mkdir(parents=True, exist_ok=True)


def nearest(data, target):
    return min(data.keys(), key=lambda a: abs(a-target))


def analyze(dataset):
    data = DataManager(config.DATA_DIR).get_data_for_dataset(dataset)
    a0 = nearest(data, 0.0); a90 = nearest(data, 90.0)
    t0, E0, tb0, Eb0 = data[a0]
    t9, E9, tb9, Eb9 = data[a90]
    t0 = np.asarray(t0); E0 = np.asarray(E0); E9 = np.asarray(E9); t9 = np.asarray(t9)
    dt = float(np.median(np.diff(t0)))   # ps assumed (project time unit)
    # cross-correlate crossed(90) vs co(0); lag>0 => 90 arrives later
    x0 = (E0-E0.mean())/(E0.std()+1e-12)
    x9 = (E9-E9.mean())/(E9.std()+1e-12)
    xc = correlate(x9, x0, mode="full")
    lags = correlation_lags(len(x9), len(x0), mode="full")
    lag_samp = int(lags[np.argmax(xc)])
    lag_ps = lag_samp*dt
    # peak-time shift (argmax|E|) as a second estimate
    dpeak_ps = float(t9[np.argmax(np.abs(E9))] - t0[np.argmax(np.abs(E0))])
    # echo search in crossed trace: second-largest local max after main peak
    aE9 = np.abs(E9); main = int(np.argmax(aE9))
    tail = aE9[main+3:]
    echo_ps = float((main+3+np.argmax(tail) - main)*dt) if len(tail) else None
    echo_rel = float(np.max(tail)/aE9[main]) if len(tail) else None
    # Delta_z predictions (P40 micrograph: 200 um)
    c = 3e8
    pred_single = 200e-6/c*1e12   # ps
    pred_etalon = 2*200e-6/c*1e12
    return {"dataset": dataset, "angle0": float(a0), "angle90": float(a90), "dt_ps": dt,
            "xcorr_lag_ps": float(lag_ps), "peak_shift_ps": dpeak_ps,
            "echo_delay_ps": echo_ps, "echo_rel_amp": echo_rel,
            "pred_single_pass_ps": pred_single, "pred_etalon_ps": pred_etalon}


def main():
    rows = [analyze(ds) for ds in ("356att", "test_grid_40_20")]
    (OUT / "hn12_depth_delay.json").write_text(json.dumps(rows, ensure_ascii=False, indent=2), encoding="utf-8")
    print("=== HN12: multi-plane depth -> delayed leakage (time-domain, model-independent) ===")
    print("Micrograph prediction (P40, Delta_z=200um): single-pass 0.67 ps / etalon 1.33 ps LATER at 90deg")
    for r in rows:
        print(f"\n{r['dataset']} (theta {r['angle0']:.0f} vs {r['angle90']:.0f}, dt={r['dt_ps']:.4f}):")
        print(f"  xcorr lag (90 vs 0)  : {r['xcorr_lag_ps']:+.3f} ps  (>0 = 90 arrives later)")
        print(f"  peak-time shift       : {r['peak_shift_ps']:+.3f} ps")
        print(f"  echo in 90 trace      : delay {r['echo_delay_ps']} ps, rel-amp {r['echo_rel_amp']}")
    print(f"\nartifacts -> {OUT}")


if __name__ == "__main__":
    main()
