"""T8 runner — build/compare M5 variants on both datasets, deliver M5 baseline."""
import sys, json
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
import model_m5 as m5

OUT = Path(__file__).resolve().parents[1] / "results" / "consolidated"
OUT.mkdir(parents=True, exist_ok=True)


def main():
    variants = [
        ("M3ref",         lambda ds: dict(D_um=m5.M3_D[ds],   D_vary=True,  leakage=False)),
        ("M5_physD",      lambda ds: dict(D_um=m5.PHYS_D[ds], D_vary=False, leakage=False)),
        ("M5_physD_leak", lambda ds: dict(D_um=m5.PHYS_D[ds], D_vary=False, leakage=True)),
        ("M5_lowD_leak",  lambda ds: dict(D_um=m5.M3_D[ds],   D_vary=False, leakage=True)),
    ]
    rows = []
    for ds in ("356att", "test_grid_40_20"):
        base = None
        for name, cfg in variants:
            r = m5.fit_variant(ds, **cfg(ds))
            r["dataset"] = ds; r["variant"] = name
            if name == "M3ref": base = r["aic"]
            r["dAIC_vs_M3ref"] = r["aic"]-base
            rows.append(r)
    (OUT / "m5_variants.json").write_text(json.dumps(rows, ensure_ascii=False, indent=2), encoding="utf-8")

    print("=== T8: M5 consolidation (masked builder, dynamic-range noise mask ON) ===")
    h = f"{'dataset':<16}{'variant':<15}{'npts':>6}{'AIC':>10}{'dAIC':>8}{'rmseDB':>8}{'DS_bias':>8}{'DS_rmse':>8}{'DW':>6}"
    print(h)
    for r in rows:
        print(f"{r['dataset']:<16}{r['variant']:<15}{r['npts']:>6}{r['aic']:>10.1f}{r['dAIC_vs_M3ref']:>8.1f}"
              f"{r['rmse_amp_db']:>8.3f}{r['deepshadow_bias_db']:>8.3f}{r['deepshadow_rmse_db']:>8.3f}"
              f"{(r['durbin_watson'] or float('nan')):>6.2f}")
    print("\nfitted params (D_um / eta0 / eta_exp / tau_leak / tau_par):")
    for r in rows:
        print(f"  {r['dataset']:<16}{r['variant']:<15} D={r['D_um']:.3f} eta0={r['eta0']:.4f} "
              f"a={r['eta_exp']:+.2f} tau_leak={r['tau_leak']:+.3f} tau_par={r['tau_par_ps']:+.3f}")
    print(f"\nartifacts -> {OUT}")


if __name__ == "__main__":
    main()
