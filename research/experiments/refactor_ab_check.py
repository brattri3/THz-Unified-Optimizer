"""A/B harness for the model_core refactor: AIC before vs after.

Runs the SAME public fit entry points that the tN_ scripts expose
(t6_hn2.fit_variant, t6_hn8.fit, model_m5.fit_variant, t20 weighted fit) and
dumps AIC/BIC/redchi/nfev + every fitted parameter to JSON. Run it once on the
pre-refactor tree, once after switching the scripts to model_core, then diff:

    .venv/Scripts/python.exe research/experiments/refactor_ab_check.py --tag pre
    ... (switch the scripts) ...
    .venv/Scripts/python.exe research/experiments/refactor_ab_check.py --tag post
    .venv/Scripts/python.exe research/experiments/refactor_ab_check.py --diff pre post

Writes to research/results/refactor/ -- never touches the real artifacts.
Read-only on originals (guardrails RESEARCH_PLAN.md sec.6).
"""
import sys, json, time, argparse, traceback, warnings
from pathlib import Path

import numpy as np

warnings.filterwarnings("ignore", category=UserWarning)
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parents[1]))

OUT = HERE.parents[0] / "results" / "refactor"
OUT.mkdir(parents=True, exist_ok=True)

DATASETS = ("att-11-16-356", "test_grid_40_20")
M3_D = {"att-11-16-356": 4.683, "test_grid_40_20": 11.447}


def _rec(res):
    """lmfit MinimizerResult -> comparable record."""
    return {"aic": float(res.aic), "bic": float(res.bic),
            "redchi": float(res.redchi), "chisqr": float(res.chisqr),
            "nfev": int(res.nfev), "nvarys": int(res.nvarys),
            "ndata": int(res.ndata), "success": bool(res.success),
            "params": {k: float(v.value) for k, v in res.params.items()}}


# fit_lib.fit_model ablation ladder (t1/t2 entry point)
FITLIB_MODELS = {
    "M0": dict(use_drude=False, use_scattering=False, free_gamma=False, use_tau_par=False),
    "M1": dict(use_drude=True,  use_scattering=False, free_gamma=False, use_tau_par=False),
    "M2": dict(use_drude=True,  use_scattering=True,  free_gamma=False, use_tau_par=False),
    "M3": dict(use_drude=True,  use_scattering=True,  free_gamma=False, use_tau_par=True),
    "M4": dict(use_drude=True,  use_scattering=True,  free_gamma=True,  use_tau_par=True),
}


def _runner(out):
    def run(name, fn):
        t0 = time.perf_counter()
        try:
            out[name] = fn()
            out[name]["_sec"] = round(time.perf_counter() - t0, 3)
            print(f"  {name:<44} AIC={out[name].get('aic', float('nan')):12.6f}"
                  f"  {out[name]['_sec']:7.2f}s")
        except Exception as exc:
            out[name] = {"error": f"{type(exc).__name__}: {exc}",
                         "trace": traceback.format_exc()[-800:]}
            print(f"  {name:<44} ERROR {type(exc).__name__}: {exc}")
    return run


def collect_fitlib():
    """Iteration-2 target: everything routed through fit_lib.fit_model."""
    import fit_lib
    out, t_start = {}, time.perf_counter()
    run = _runner(out)
    print("fit_lib.fit_model (M0..M4 ladder):")
    for ds in list(DATASETS) + ["specac"]:
        for label, cfg in FITLIB_MODELS.items():
            def one(ds=ds, label=label, cfg=cfg):
                m, arr = fit_lib.fit_model(ds, label=label, **cfg)
                # drop the heavy arrays; keep a checksum of the residual map so a
                # silent change in the post-fit grid cannot slip through
                fin = arr["amp_res_db"][np.isfinite(arr["amp_res_db"])]
                m["_ampres_sum"] = float(np.sum(fin))
                m["_ampres_sumsq"] = float(np.sum(fin**2))
                m["_theo_absum"] = float(np.sum(np.abs(arr["theo"])))
                return m
            run(f"fit_lib/{ds}/{label}", one)
    print(f"\ntotal {time.perf_counter()-t_start:.1f}s")
    return out


def collect():
    import t6_hn2, t6_hn8, model_m5, t20_weighted_fit as t20
    out, t_start = {}, time.perf_counter()
    run = _runner(out)

    print("t6_hn2.fit_variant:")
    for ds in DATASETS:
        run(f"t6_hn2/{ds}/M3", lambda ds=ds: t6_hn2.fit_variant(ds, leakage=False))
        run(f"t6_hn2/{ds}/HN2", lambda ds=ds: t6_hn2.fit_variant(ds, leakage=True))
        run(f"t6_hn2/{ds}/HN2_fixD",
            lambda ds=ds: t6_hn2.fit_variant(ds, leakage=True, fix_D=M3_D[ds]))

    print("t6_hn8.fit:")
    for ds in DATASETS:
        for fixD in (True, False):
            for mode in ("m3", "hn2", "hn8"):
                run(f"t6_hn8/{ds}/{mode}/fixD={fixD}",
                    lambda ds=ds, m=mode, f=fixD: t6_hn8.fit(ds, m, f))

    print("model_m5.fit_variant (via t8 variant table):")
    variants = [("M3ref", lambda ds: dict(D_um=model_m5.M3_D[ds], D_vary=True, leakage=False)),
                ("M5_physD", lambda ds: dict(D_um=model_m5.PHYS_D[ds], D_vary=False, leakage=False)),
                ("M5_physD_leak", lambda ds: dict(D_um=model_m5.PHYS_D[ds], D_vary=False, leakage=True)),
                ("M5_lowD_leak", lambda ds: dict(D_um=model_m5.M3_D[ds], D_vary=False, leakage=True))]
    for ds in DATASETS:
        for nm, cfg in variants:
            run(f"model_m5/{ds}/{nm}", lambda ds=ds, c=cfg: model_m5.fit_variant(ds, **c(ds)))

    print("t20_weighted_fit (multi-rep):")
    for ds in ("specac", "purewave", "att-11-16-356"):
        def t20_fit(ds=ds):
            G = t20.build_weighted(ds)
            if G is None:
                return {"skipped": "no usable multi-rep data"}
            import fit_lib
            P_um = fit_lib.GEOMETRY[ds]["P_um"]
            D_phys = fit_lib.GEOMETRY[ds]["D_phys_um"]
            m3 = t20.fit(G, P_um, D_um=min(D_phys, P_um - 1.0), D_vary=True, leak=False)
            G["D_eff_pt"] = float(m3.params["D_um"].value)
            m5 = t20.fit(G, P_um, D_um=G["D_eff_pt"], D_vary=True, leak=True)
            r = {"M3": _rec(m3), "M5": _rec(m5), "aic": float(m5.aic),
                 "valid_points": int(np.sum(G["valid"])),
                 "sigma_floor": float(G["sigma_floor"])}
            return r
        run(f"t20/{ds}", t20_fit)

    print(f"\ntotal {time.perf_counter()-t_start:.1f}s")
    return out


def flatten(d, prefix=""):
    """dict-of-dicts -> {'a/b/c': scalar} for numeric leaves."""
    flat = {}
    for k, v in d.items():
        key = f"{prefix}{k}"
        if isinstance(v, dict):
            flat.update(flatten(v, key + "/"))
        elif isinstance(v, (int, float)) and not isinstance(v, bool):
            flat[key] = float(v)
    return flat


def diff(tag_a, tag_b):
    A = json.loads((OUT / f"{tag_a}.json").read_text(encoding="utf-8"))
    B = json.loads((OUT / f"{tag_b}.json").read_text(encoding="utf-8"))
    fa, fb = flatten(A), flatten(B)
    keys = sorted(set(fa) | set(fb))
    worst, nbit, ndiff, nmiss = [], 0, 0, 0
    for k in keys:
        if k.endswith("/_sec") or k == "_sec":
            continue
        if k not in fa or k not in fb:
            nmiss += 1
            print(f"  MISSING {k} ({'only ' + tag_b if k not in fa else 'only ' + tag_a})")
            continue
        x, y = fa[k], fb[k]
        if x == y:
            nbit += 1
            continue
        ndiff += 1
        rel = abs(x - y) / max(abs(x), 1e-300)
        worst.append((rel, k, x, y))
    worst.sort(reverse=True)
    print(f"\n=== diff {tag_a} -> {tag_b} ===")
    print(f"  identical (bit-for-bit): {nbit}")
    print(f"  differing              : {ndiff}")
    print(f"  missing on one side    : {nmiss}")
    for rel, k, x, y in worst[:25]:
        print(f"    {k:<52} {x!r} -> {y!r}  rel={rel:.3e}")
    # AIC-only summary
    aics = [(k, fa[k], fb[k]) for k in keys
            if k.endswith("aic") and k in fa and k in fb]
    bad = [(k, x, y) for k, x, y in aics if x != y]
    print(f"\n  AIC entries: {len(aics)}, differing: {len(bad)}")
    for k, x, y in bad:
        print(f"    AIC {k}: {x!r} -> {y!r}  (d={y-x:+.6e})")
    return 0 if (ndiff == 0 and nmiss == 0) else 1


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--tag", help="collect and save as <tag>.json")
    ap.add_argument("--diff", nargs=2, metavar=("A", "B"))
    ap.add_argument("--set", choices=("core", "fitlib"), default="core",
                    help="core = t6/t8/model_m5/t20 (iter 1); fitlib = fit_lib.fit_model (iter 2)")
    a = ap.parse_args()
    if a.diff:
        return diff(*a.diff)
    if not a.tag:
        ap.error("need --tag or --diff")
    res = collect_fitlib() if a.set == "fitlib" else collect()
    (OUT / f"{a.tag}.json").write_text(json.dumps(res, ensure_ascii=False, indent=2),
                                       encoding="utf-8")
    print(f"artifacts -> {OUT / (a.tag + '.json')}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
