"""model_core — cached, vectorised Blanco core shared by the research fits.

WHY: every fit script (fit_lib, t6_hn2, t6_hn8, model_m5, t20) carried its own
near-identical copy of `compute_theoretical_grid_2d`, and each copy re-ran the
per-frequency Python loop over `model_blanco.compute_t_perp/compute_t_par` on
EVERY `nfev`. Blanco's t_perp/t_par depend only on (nu-grid, P, D, N, drude) --
never on loss_factor / gamma / angle_offset / tau_ps / tau_par_ps / leakage --
so during lmfit's numerical Jacobian they are recomputed identically for every
column that does not touch D.

This module fixes both:
  1. `blanco_t()`  -- vectorised over frequency + memoised on (nu, P, D, N, drude);
  2. `compute_grid()` -- ONE grid function that is a strict superset of
     `compute_theoretical_grid_2d` (leakage off) / `compute_grid_hn2`
     (tau_leak=0) / `compute_grid_hn8`;
  3. `complex_residual()` / `residual_from_params()` -- one weighted-residual.

NUMERICAL CONTRACT: the vectorised code performs *exactly* the same float64
operations, in the same order, as the scalar originals in
`unified_optimizer/model_blanco.py` (the m-sums are still accumulated in a
Python loop m=1..N so the summation order is preserved; every scalar `if`
guard became an `np.where` on the same predicate). Equality is asserted
bit-for-bit by `verify_model_core.py`. Originals are NOT edited (guardrails
RESEARCH_PLAN.md sec.6) -- this is a research-side copy.

Cache keys use exact float bits (no rounding): lmfit perturbs D by ~1e-8
relative for the Jacobian, and rounding the key would silently merge those
neighbours and corrupt the derivative.
"""
from __future__ import annotations

import sys
from collections import OrderedDict
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[2]      # THz-Unified-Optimizer
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))

from unified_optimizer import config, model_blanco   # noqa: E402

EPS = model_blanco.EPS
C_LIGHT = model_blanco.C_LIGHT
Z0 = model_blanco.Z0

# weights of the amplitude / phase halves of the residual vector (kept in sync
# with fit_lib.W_AMP / W_PHASE, which every tN_ script imports)
W_AMP = 1.0
W_PHASE = float(np.sqrt(0.1))

# ---------------------------------------------------------------------------
# cache
# ---------------------------------------------------------------------------
_T_CACHE: "OrderedDict[tuple, tuple]" = OrderedDict()
_CACHE_MAX = 512            # ~5.4 kB/entry at 169 freqs -> a few MB worst case
_CACHE_STATS = {"hits": 0, "misses": 0}


def cache_stats() -> dict:
    """{'hits', 'misses', 'size'} for the Blanco t_perp/t_par memo."""
    return {**_CACHE_STATS, "size": len(_T_CACHE)}


def clear_cache() -> None:
    _T_CACHE.clear()
    _CACHE_STATS.update(hits=0, misses=0)


# ---------------------------------------------------------------------------
# vectorised Blanco
# ---------------------------------------------------------------------------
def _cpython_cdiv(a, b) -> np.ndarray:
    """Complex division bit-identical to CPython's `complex.__truediv__`.

    NumPy and CPython both use Smith's method but differ in the last step
    (NumPy scales by the reciprocal of the denominator, CPython divides), so
    they disagree by ~1 ULP. NumPy *scalar* and *array* complex division do
    agree with each other -- verified over 4000 random pairs -- so this helper
    is only needed where the scalar reference code operates on plain Python
    complex objects, i.e. inside get_drude_impedance_normalized. Everything
    else in model_blanco already runs on np.complex128 (compute_C returns
    np.sqrt(...)), and matches the vectorised path exactly.
    """
    a = np.asarray(a, dtype=complex)
    b = np.asarray(b, dtype=complex)
    ar, ai, br, bi = a.real, a.imag, b.real, b.imag
    big = np.abs(br) >= np.abs(bi)
    with np.errstate(divide="ignore", invalid="ignore"):
        # branch |b.real| >= |b.imag|: divide top and bottom by b.real
        r1 = bi / np.where(big, br, 1.0)
        d1 = br + bi * r1
        re1 = (ar + ai * r1) / d1
        im1 = (ai - ar * r1) / d1
        # branch |b.imag| > |b.real|: divide top and bottom by b.imag
        r2 = br / np.where(big, 1.0, bi)
        d2 = br * r2 + bi
        re2 = (ar * r2 + ai) / d2
        im2 = (ai * r2 - ar) / d2
    out = np.empty(np.broadcast(a, b).shape, dtype=complex)
    out.real = np.where(big, re1, re2)
    out.imag = np.where(big, im1, im2)
    return out


def _cmul(a, b) -> np.ndarray:
    """complex*complex bit-identical to NumPy's SCALAR multiplication.

    The complex128 multiply *ufunc* contracts `ar*br - ai*bi` into an FMA on
    x86-64; the scalar (np.complex128) path does not, so the two disagree by
    ~1 ULP on ~45% of random inputs. Computing the four real products as
    separate float64 ufunc calls reproduces the scalar result exactly
    (0 mismatches over 20k random pairs). Complex *division*, addition and
    real*complex already agree scalar-vs-array, so only products need this.
    """
    a = np.asarray(a, dtype=complex)
    b = np.asarray(b, dtype=complex)
    ar, ai, br, bi = a.real, a.imag, b.real, b.imag
    out = np.empty(np.broadcast(a, b).shape, dtype=complex)
    out.real = ar * br - ai * bi
    out.imag = ar * bi + ai * br
    return out


def _csq(a) -> np.ndarray:
    """a**2 -- matches NumPy's scalar `z**2` exactly (array `**2` does not)."""
    return _cmul(a, a)


def _drude_z(freqs_thz: np.ndarray, use_drude: bool) -> np.ndarray:
    """Vectorised model_blanco.get_drude_impedance_normalized (bit-exact)."""
    if not use_drude:
        return np.zeros(freqs_thz.shape, dtype=complex)
    omega = freqs_thz * 2 * np.pi * 1e12
    sigma0 = 1.8e7      # S/m (tungsten)
    tau = 8.0e-15       # s
    mu0 = 4 * np.pi * 1e-7
    # these two divisions are Python-complex in the scalar reference -> _cpython_cdiv
    sigma_omega = _cpython_cdiv(sigma0, 1.0 - 1j * omega * tau)
    with np.errstate(divide="ignore", invalid="ignore"):
        Z_s = np.sqrt(_cpython_cdiv(1j * omega * mu0, sigma_omega))
    return np.where(freqs_thz > 0, Z_s / Z0, 0j)   # np scalar/array div agree


def _blanco_tt(freqs_thz: np.ndarray, p: float, d: float, N: int,
               use_drude: bool):
    """(t_perp, t_par) over the whole frequency grid. p, d in METRES."""
    d_over_p = d / p
    lam = C_LIGHT / (freqs_thz * 1e12)
    pol = p / lam                                    # P/lambda per frequency
    Z_drude = _drude_z(freqs_thz, use_drude)
    ones = np.ones(freqs_thz.shape)

    if d_over_p <= 0 or d_over_p >= 1:
        # degenerate branch of compute_fa/fb/fc/fd (scalar early returns)
        fa = ones * 1e6
        fb = ones * -1e6
        fc = ones * 0.0
        fd = ones * 0.0
    else:
        pidl = np.pi * d_over_p * pol                # pi * d/p * p/lambda
        log_arg = 1.0 / (np.pi * d_over_p)
        slog = model_blanco.safe_log(log_arg)        # scalar, freq-independent

        # m-sums: loop preserves the accumulation order of the scalar code.
        # sum_inv -> shared by compute_A1 and compute_fc
        # sum_a2  -> compute_A2
        sum_inv = np.zeros(freqs_thz.shape, dtype=complex)
        sum_a2 = np.zeros(freqs_thz.shape, dtype=complex)
        for m in range(1, N + 1):
            C_m = np.sqrt((m**2 - pol**2) + 1j * 1e-9)   # compute_C
            sum_inv = sum_inv + (1.0 / C_m - 1.0 / m)
            sum_a2 = sum_a2 + (m - 0.5 / m * (pol)**2 - C_m)

        # compute_A1 / compute_A2 (evaluated once each, unlike the scalar path
        # where compute_fa and compute_fb both recomputed A2)
        A1 = 1.0 + 0.5 * (pidl)**2 * (slog + 0.75) + 0.5 * (pidl)**2 * sum_inv
        A2 = (1.0
              + 0.5 * (pidl)**2 * (11.0 / 4.0 - slog)
              + (1.0 / 24.0) * (pidl)**2
              - ((pidl)**2 * sum_a2))

        ok_a2 = np.abs(A2) >= EPS
        A2s = np.where(ok_a2, A2, 1.0)               # safe denominator
        ok_pol = np.abs(pol) >= EPS
        pols = np.where(ok_pol, pol, 1.0)

        # compute_fa
        fa = np.where(ok_a2, 0.5 * pol * (np.pi * d_over_p)**2 / A2s, 1e6)
        # compute_fb
        fb_term1 = 2.0 / pols * (1.0 / (np.pi * d_over_p))**2 * A1
        fb_term2 = - pol / 4.0 * (np.pi * d_over_p)**2 / A2s
        fb = np.where(ok_a2 & ok_pol, fb_term1 + fb_term2, -1e6)
        # compute_fc / compute_fd
        fc = pol * (slog + sum_inv)
        fd = pol * (np.pi * d_over_p)**2

    # ---- compute_t_perp -----------------------------------------------------
    ok_fa = np.abs(fa) > EPS
    ok_fb = np.abs(fb) > EPS
    Za = np.where(ok_fa, -1j / np.where(ok_fa, fa, 1.0) + Z_drude, -1e9j)
    Zb = np.where(ok_fb, -1j / np.where(ok_fb, fb, 1.0) + Z_drude, -1e9j)

    za2 = _csq(Za)
    zazb = _cmul(Za, Zb)
    num1 = za2 + zazb + _cmul(za2, Zb)
    den1 = 2.0 * Za + za2 + Zb + zazb
    ok_den1 = np.abs(den1) > EPS
    Z1 = np.where(ok_den1, num1 / np.where(ok_den1, den1, 1.0), 0j)
    ok_za1 = np.abs(1.0 + Za) > EPS
    Z2 = np.where(ok_za1, Za / np.where(ok_za1, 1.0 + Za, 1.0), 1.0 + 0j)

    num_t = _cmul(2.0 * Z1, Z2)
    den_t = _cmul(1.0 + Z1, Zb + Z2)
    ok_dt = np.abs(den_t) > EPS
    t_perp = np.where(ok_dt, num_t / np.where(ok_dt, den_t, 1.0), 0j)

    # ---- compute_t_par ------------------------------------------------------
    Zc = _cmul(1j, fc) + Z_drude
    Zd = _cmul(-1j, fd) + Z_drude
    den34 = 1.0 + Zc + Zd
    ok34 = np.abs(den34) >= EPS
    den34s = np.where(ok34, den34, 1.0)
    Z3 = (Zd + _cmul(2.0 * Zc, Zd) + _csq(Zd) + Zc) / den34s
    Z4 = (Zc + _cmul(Zc, Zd)) / den34s
    num_p = _cmul(2.0 * Z3, Z4)
    den_p = _cmul(_cmul(1.0 + Z3, Zd + Z4), 1.0 + Zd)
    ok_dp = np.abs(den_p) > EPS
    t_par = np.where(ok34 & ok_dp, num_p / np.where(ok_dp, den_p, 1.0), 0j)

    return t_perp, t_par


def blanco_t(freqs_thz, p: float, d: float, N: int = 15,
             use_drude: bool = True):
    """Memoised (t_perp, t_par) for one wire grid. p, d in METRES.

    Returns READ-ONLY arrays -- they are the cached objects; copy before
    mutating. Cache key uses exact float bits (see module docstring).
    """
    freqs_thz = np.ascontiguousarray(freqs_thz, dtype=float)
    key = (freqs_thz.tobytes(), freqs_thz.shape,
           float(p), float(d), int(N), bool(use_drude))
    hit = _T_CACHE.get(key)
    if hit is not None:
        _CACHE_STATS["hits"] += 1
        _T_CACHE.move_to_end(key)
        return hit
    _CACHE_STATS["misses"] += 1

    t_perp, t_par = _blanco_tt(freqs_thz, p, d, N, use_drude)
    t_perp.flags.writeable = False
    t_par.flags.writeable = False
    out = (t_perp, t_par)
    if len(_T_CACHE) >= _CACHE_MAX:
        _T_CACHE.popitem(last=False)      # LRU eviction
    _T_CACHE[key] = out
    return out


# ---------------------------------------------------------------------------
# unified theoretical grid
# ---------------------------------------------------------------------------
def compute_grid(angles_deg, freqs_thz, p, d, loss_factor, angle_offset, tau_ps,
                 gamma=1.0, N=15, use_drude=True, tau_par_ps=0.0,
                 eta0=0.0, eta_exp=1.0, delta0=0.0, tau_leak=0.0):
    """Complex model grid [angle, freq] -- superset of all research copies.

    eta0 == 0            -> identical to optimizer_2d.compute_theoretical_grid_2d
    tau_leak == 0        -> identical to t6_hn2.compute_grid_hn2
    (general)            -> identical to t6_hn8.compute_grid_hn8

    Leakage channel added to the parallel component (HN2/HN8):
        t_par += eta0 * nu**eta_exp * exp(1j*(delta0 + 2*pi*nu*tau_leak))
    """
    freqs_thz = np.ascontiguousarray(freqs_thz, dtype=float)
    t_perp_arr, t_par_arr = blanco_t(freqs_thz, p, d, N, use_drude)

    if eta0 != 0.0:
        eta = eta0 * (freqs_thz ** eta_exp)
        delta = delta0 + 2 * np.pi * freqs_thz * tau_leak
        t_par_arr = t_par_arr + eta * np.exp(1j * delta)

    adj = np.deg2rad(np.asarray(angles_deg) - angle_offset)
    cos_a = np.cos(adj)[:, None]
    sin_a = np.sin(adj)[:, None]

    if config.USE_POWER_LAW:
        loss_amp = np.exp(-0.5 * (loss_factor / 4.343) * (freqs_thz ** gamma))[None, :]
    else:
        loss_amp = np.exp(-0.5 * loss_factor * freqs_thz)[None, :]
    phase_perp = np.exp(-1j * 2 * np.pi * freqs_thz * tau_ps)[None, :]
    phase_par = np.exp(-1j * 2 * np.pi * freqs_thz * (tau_ps + tau_par_ps))[None, :]

    return (cos_a**2 * (t_perp_arr[None, :] * loss_amp * phase_perp)
            + sin_a**2 * (t_par_arr[None, :] * loss_amp * phase_par))


# parameter names understood by grid_from_params, with the defaults that make
# compute_grid collapse to the plain Blanco model
_GRID_DEFAULTS = {"loss_factor": 0.0, "angle_offset": 0.0, "tau_ps": 0.0,
                  "gamma": 1.0, "tau_par_ps": 0.0,
                  "eta0": 0.0, "eta_exp": 1.0, "delta0": 0.0, "tau_leak": 0.0}


def grid_from_params(params, angles_deg, freqs_thz, N=15, use_drude=True):
    """compute_grid driven by an lmfit.Parameters (missing names -> defaults)."""
    def g(name):
        return params[name].value if name in params else _GRID_DEFAULTS[name]
    return compute_grid(
        angles_deg, freqs_thz,
        params["P_um"].value * 1e-6, params["D_um"].value * 1e-6,
        g("loss_factor"), g("angle_offset"), g("tau_ps"),
        gamma=g("gamma"), N=N, use_drude=use_drude, tau_par_ps=g("tau_par_ps"),
        eta0=g("eta0"), eta_exp=g("eta_exp"), delta0=g("delta0"),
        tau_leak=g("tau_leak"))


# ---------------------------------------------------------------------------
# unified residual
# ---------------------------------------------------------------------------
def complex_residual(exp, theo, valid, w=None, w_amp=W_AMP, w_phase=W_PHASE):
    """[amp_residual*w*w_amp, wrapped_phase_residual*w*w_phase] over `valid`.

    `w` (optional) is a full-grid array of per-point weights (t20: 1/sigma
    normalised to mean 1 over the valid mask).
    """
    em, tm = exp[valid], theo[valid]
    amp = np.abs(em) - np.abs(tm)
    ph = np.angle(em) - np.angle(tm)
    ph = np.arctan2(np.sin(ph), np.cos(ph))
    if w is not None:
        ww = w[valid]
        amp = amp * ww
        ph = ph * ww
    return np.concatenate([amp * w_amp, ph * w_phase])


def residual_from_params(params, angles_deg, freqs_thz, exp, valid, w=None,
                         N=15, use_drude=True, w_amp=W_AMP, w_phase=W_PHASE):
    """lmfit residual: guards D>=P, builds the grid, returns the weighted vector."""
    p = params["P_um"].value * 1e-6
    d = params["D_um"].value * 1e-6
    if d >= p:
        return np.ones(int(np.sum(valid)) * 2) * 1e6
    theo = grid_from_params(params, angles_deg, freqs_thz, N=N, use_drude=use_drude)
    return complex_residual(exp, theo, valid, w=w, w_amp=w_amp, w_phase=w_phase)
