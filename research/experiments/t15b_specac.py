"""T15b — full M5 fit on the Specac external sample; D_eff/D_phys + leakage transfer.

Specac micrograph geometry: P=24.9 um (pinned), D_phys~14 um (wide prior), D/P~0.55.
Variants (reuse model_m5.fit_variant, which takes D_um directly):
  M3ref         : D free, no leakage       -> effective D_eff (no degeneracy w/o leak)
  M5_physD      : D pinned = D_phys(14)     -> does physical pin over-predict (as 356att)?
  M5_physD_leak : D pinned = 14 + leakage
  M5_lowD_leak  : D pinned = M3 D_eff + leakage  -> does the leakage cluster transfer?
Reports D_eff/D_phys for cross-sample D_eff-law check and eta0/tau_leak transfer.
Read-only originals; output results/validation/.
"""
import sys, json
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
import model_m5 as m5

OUT = Path(__file__).resolve().parents[1] / "results" / "validation"
OUT.mkdir(parents=True, exist_ok=True)
DS = "specac"
D_PHYS = m5.fit_lib.GEOMETRY[DS]["D_phys_um"]   # 14.0
P = m5.fit_lib.GEOMETRY[DS]["P_um"]             # 24.9


def main():
    m3 = m5.fit_variant(DS, D_um=D_PHYS, D_vary=True, leakage=False)
    Deff = m3["D_um"]
    physD = m5.fit_variant(DS, D_um=D_PHYS, D_vary=False, leakage=False)
    physD_leak = m5.fit_variant(DS, D_um=D_PHYS, D_vary=False, leakage=True)
    lowD_leak = m5.fit_variant(DS, D_um=Deff, D_vary=False, leakage=True)
    for r, name in ((m3, "M3ref"), (physD, "M5_physD"), (physD_leak, "M5_physD_leak"), (lowD_leak, "M5_lowD_leak")):
        r["variant"] = name
    base = m3["aic"]
    rows = [m3, physD, physD_leak, lowD_leak]
    for r in rows:
        r["dAIC_vs_M3ref"] = r["aic"]-base
        r["Deff_over_Dphys"] = r["D_um"]/D_PHYS

    summary = {"dataset": DS, "P_um": P, "D_phys_um": D_PHYS, "D_over_P_phys": D_PHYS/P,
               "M3_Deff_um": Deff, "M3_Deff_over_Dphys": Deff/D_PHYS,
               "cross_sample_Deff_over_Dphys": {
                   "356att(D/P0.71)": round(4.683/11.0, 3),
                   "test_grid(D/P0.46)": round(11.447/18.0, 3),
                   f"specac(D/P{D_PHYS/P:.2f})": round(Deff/D_PHYS, 3)},
               "variants": rows}
    (OUT / "specac_m5.json").write_text(json.dumps(summary, ensure_ascii=False, indent=2), encoding="utf-8")

    print(f"=== T15b: M5 on Specac (P={P} pinned, D_phys={D_PHYS}, D/P={D_PHYS/P:.2f}) ===")
    print(f"{'variant':<15}{'AIC':>10}{'dAIC':>8}{'D_um':>8}{'Deff/Dphys':>11}{'DS_bias':>8}{'DS_rmse':>8}{'eta0':>8}{'tau_lk':>8}")
    for r in rows:
        print(f"{r['variant']:<15}{r['aic']:>10.1f}{r['dAIC_vs_M3ref']:>8.1f}{r['D_um']:>8.3f}"
              f"{r['Deff_over_Dphys']:>11.3f}{r['deepshadow_bias_db']:>8.3f}{r['deepshadow_rmse_db']:>8.3f}"
              f"{r['eta0']:>8.4f}{r['tau_leak']:>8.3f}")
    print(f"\nD_eff/D_phys across samples (effective diameter law):")
    for k, v in summary["cross_sample_Deff_over_Dphys"].items():
        print(f"  {k:<22} {v}")
    print(f"\nartifacts -> {OUT}")


if __name__ == "__main__":
    main()
