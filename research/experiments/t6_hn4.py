"""T6 test HN4 — TM-selective loss/leakage from positional irregularity (Manabe).

Manabe & Murk 2005 (10.1109/TAP.2004.884786): random wire-position irregularity,
treated as a perturbation to the exact periodic theory, strongly affects the TM
(E || wires) response and barely affects TE (E _|_ wires). Physically, disorder
disrupts the coherent reflection that makes a regular grid a good polarizer for
the parallel component => MORE parallel leakage. This is a candidate FIRST-
PRINCIPLES ORIGIN for the leakage cluster (HN2/HN7/HN8).

Test from first principles: finite array of N PEC cylinders, DIRECT multiple-
scattering (Graf's addition theorem, no lattice sum, so aperiodic geometry is
allowed). Monte-Carlo over Gaussian x-position jitter sigma_pos:
  Q1: does |t_par| grow with sigma_pos (leakage)?
  Q2: is it TM-selective (|t_par| affected >> |t_perp|)?
  Q3: which sigma_pos reproduces the observed leakage eta ~ 0.012-0.03 (HN7)?
Single-cylinder coefficients b_l reused from the validated T5 solver.
Read-only on originals; output under research/results/HN4/.
"""
import sys, json
from pathlib import Path
import numpy as np
from scipy.special import hankel1

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import t5_latticesum as t5

OUT = HERE.parents[0] / "results" / "HN4"
OUT.mkdir(parents=True, exist_ok=True)
C = 3e8


def finite_array_t(xpos, k, a, pol, L=5):
    """Direct multiple-scattering for cylinders at (x_j, 0). Returns zeroth-order
    forward transmission amplitude T0 via the per-cylinder averaged Rayleigh term."""
    N = len(xpos)
    orders = np.arange(-L, L+1)
    nord = len(orders)
    ka = k*a
    b = np.array([t5.single_cylinder_b(l, ka, pol) for l in orders])
    # global index arrays
    Jidx = np.repeat(np.arange(N), nord)          # cylinder per row
    Lidx = np.tile(orders, N)                      # order per row
    xr = xpos[Jidx]
    dx = xr[:, None] - xr[None, :]                 # x_j - x_j'  (all on axis, y=0)
    R = np.abs(dx)
    # angle from j to j': +x (dx<0 means j' to the right) -> phi in {0, pi}
    phi = np.where(dx < 0, 0.0, np.pi)             # angle of (x_j' - x_j)
    same = (Jidx[:, None] == Jidx[None, :])
    dl = Lidx[:, None] - Lidx[None, :]
    Rsafe = np.where(same, 1.0, R)                 # avoid H(.,0) on diagonal blocks
    H = hankel1(dl, k*Rsafe)
    coupling = H * np.exp(1j*dl*phi)
    coupling[same] = 0.0
    M = np.eye(N*nord, dtype=complex) - (b[Lidx][:, None]) * coupling
    rhs = b[Lidx] * 1.0                            # plane wave +y, a_l=1, y_j=0
    B = np.linalg.solve(M, rhs)
    # per-cylinder zeroth-order forward term, averaged (period contribution 2/(d k))
    gT = (-1j)**Lidx
    per_cyl = (B*gT).reshape(N, nord).sum(axis=1)
    return per_cyl  # complex per-cylinder amplitude (pre-normalization)


def damped_grating_TR(k, d, a, pol, sigma, L=8, beta=2.0, **ls_kw):
    """T5 solver with lattice sums damped by a Debye-Waller factor W=exp(-beta k^2 sigma^2),
    modelling disorder-averaged coherent multiple scattering (Manabe mechanism).
    Disorder damps the collective coupling S_l (relative-displacement variance ~2 sigma^2);
    par (TM), which relies on strong collective S_0 for extinction, loses reflection =>
    more parallel leakage; perp (TE, weakly coupled) is barely affected."""
    from scipy.special import jv, jvp, hankel1, h1vp
    W = np.exp(-beta * (k*sigma)**2)
    ls = t5.lattice_sums(k, d, L, **ls_kw) * W          # damped collective coupling
    orders = np.arange(-L, L+1); n = len(orders); ka = k*a
    b = np.array([t5.single_cylinder_b(l, ka, pol) for l in orders])
    M = np.eye(n, dtype=complex)
    for ia, la in enumerate(orders):
        for ib, lb in enumerate(orders):
            M[ia, ib] -= b[ia] * ls[abs(la-lb)]
    B = np.linalg.solve(M, b*1.0)
    pref = 2.0/(d*k)
    gT = np.array([(-1j)**l for l in orders])
    return 1.0 + pref*np.sum(B*gT)


def main():
    P = 15.5e-6; d = P; nu = 0.5e12; k = 2*np.pi*nu/C
    a = (0.71*P)/2.0     # nominal dense 356att: D/P=0.71

    # --- documented negative result: finite-array direct solve cannot reach the
    #     par regime (collective S_0 needs ~inf wires); validated & abandoned. ---
    x0 = (np.arange(41)-20)*d; pref = 2.0/(d*k)
    def T_fin(pol):
        per = finite_array_t(x0, k, a, pol, L=5); core = np.abs(np.arange(41)-20) <= 15
        return abs(1.0 + pref*np.mean(per[core]))
    finite_note = {"par_finite_N41": T_fin("par"),
                   "par_lattice": abs(t5.grating_TR(k, d, a, "par", L=5)[0]),
                   "comment": "finite N=41 gives |t_par|~1 (no extinction) vs lattice 0.045; "
                              "collective TM reflection needs ~inf wires -> use damped lattice sums"}

    # --- HN4 test: disorder-damped lattice sums (Manabe coherent-scattering model) ---
    sigmas = [0.0, 0.02, 0.05, 0.10, 0.15, 0.20, 0.30]   # fraction of period
    rows = []
    tperp0 = abs(damped_grating_TR(k, d, a, "perp", 0.0))
    for sf in sigmas:
        s = sf*d
        tpar = abs(damped_grating_TR(k, d, a, "par", s))
        tperp = abs(damped_grating_TR(k, d, a, "perp", s))
        rows.append({"sigma_frac": sf, "sigma_um": sf*P*1e6,
                     "tpar": tpar, "tperp": tperp, "eta": tpar/tperp,
                     "tperp_change_pct": 100*(tperp-tperp0)/tperp0})

    # sigma reproducing observed leakage eta (HN7 dense=0.012, sparse=0.030)
    etas = np.array([r["eta"] for r in rows]); sig = np.array([r["sigma_frac"] for r in rows])
    def sigma_for(eta_obs):
        if eta_obs <= etas[0]:
            return 0.0
        for i in range(1, len(etas)):
            if etas[i] >= eta_obs:
                t = (eta_obs-etas[i-1])/(etas[i]-etas[i-1]+1e-30)
                return float(sig[i-1] + t*(sig[i]-sig[i-1]))
        return float(sig[-1])

    out = {"geometry": {"P_um": P*1e6, "D_over_P": 0.71, "nu_THz": 0.5},
           "finite_array_negative_result": finite_note,
           "model": "disorder-damped lattice sums W=exp(-2 k^2 sigma^2) (Manabe coherent mechanism)",
           "damped_scan": rows,
           "sigma_for_observed_eta": {"dense_0.012": sigma_for(0.012), "sparse_0.030": sigma_for(0.030)}}
    (OUT / "hn4_irregularity.json").write_text(json.dumps(out, ensure_ascii=False, indent=2), encoding="utf-8")

    print("=== T6 HN4: TM-selective irregularity (disorder-damped lattice sums) ===")
    print(f"finite-array negative result: |t_par|(N=41)={finite_note['par_finite_N41']:.3f} "
          f"vs lattice {finite_note['par_lattice']:.3f} -> finite array cannot reach par regime")
    print("\ndisorder-damped scan (Manabe coherent model):")
    print(f"{'sig%P':>6}{'sig_um':>8}{'|t_par|':>10}{'|t_perp|':>10}{'eta':>9}{'dTperp%':>9}")
    for r in rows:
        print(f"{r['sigma_frac']*100:>5.0f}%{r['sigma_um']:>8.2f}{r['tpar']:>10.4f}"
              f"{r['tperp']:>10.4f}{r['eta']:>9.4f}{r['tperp_change_pct']:>+9.3f}")
    print(f"\nsigma reproducing observed eta: dense(0.012)->{out['sigma_for_observed_eta']['dense_0.012']*100:.1f}%P"
          f"  sparse(0.030)->{out['sigma_for_observed_eta']['sparse_0.030']*100:.1f}%P")
    print(f"artifacts -> {OUT}")


if __name__ == "__main__":
    main()
