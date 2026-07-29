"""A1 — TWO-WGP forward model (Jones product, frequency-complex, per-element).

WHY (see research/two_wgp/PLAN.md):
    The production 2D model (`unified_optimizer/optimizer_2d.py:68`)
        E_out = cos^2(theta)*t_perp_eff + sin^2(theta)*t_par_eff
    is exactly [ideal x-polarizer] -> [real WGP1(theta)] -> [IDEAL x-analyzer].
    The real rig has TWO REAL WGPs:
        E_in (linear) -> WGP1(theta, rotating sample) -> WGP2(phi2, fixed) -> detector(axis d)
    so WGP2's non-ideality was being absorbed into WGP1's fitted geometry (D_eff).
    This module implements the correct chain:
        E_out(theta,nu) = d^T . J_WGP2(phi2,nu) . J_WGP1(theta,nu) . E_in

CONVENTIONS (explicit — they were the #1 pitfall of Phases 1-3):
    * Every angle is the orientation of an element's TRANSMISSION axis (the
      perpendicular-to-wires axis, which carries t_perp) measured in the lab
      frame, in DEGREES, CCW positive.  Wires are therefore at angle+90.
    * Jones matrix of a WGP whose transmission axis sits at phi:
          J(phi) = R(phi) @ diag(t_perp, t_par) @ R(-phi)
          R(phi) = [[cos, -sin], [sin, cos]]
      -> J = [[a c^2 + b s^2, (a-b) c s], [(a-b) c s, a s^2 + b c^2]]  (symmetric)
    * `det_deg` is the detector's sensitivity axis d; the detector reads the
      COMPLEX projection d^T E_out (THz-TDS is field-sensitive, not power).
    * `e_in_deg` is the emitter's linear polarization direction.
    * OPEN OWNER QUESTION (PLAN sec.9): whether WGP2's WIRES or its
      TRANSMISSION axis is the one aligned with d.  Here phi2_deg is the
      TRANSMISSION axis, so "wires along d" == phi2_deg = 90.  Both are
      reachable; nothing is hard-coded.

Ideal limits this module reproduces (validated in `selftest`):
    * t=(1,0) on both, E_in along WGP1 axis  ->  power = cos^2(theta1-theta2)  (Malus)
    * t=(1,0) on both, E_in = d = phi2 = 0   ->  amplitude = cos^2(theta)
      (identical to the production one-WGP model in its ideal limit — the two
      models only diverge once t_par != 0, which is the whole point)
    * identical WGPs, phi2=0, E_in=d=x  ->  bit-exact match to
      `model_blanco.transmission_complex_two_polarizers`

Read-only w.r.t. originals: imports `unified_optimizer.model_blanco` /`config`
but never edits them.
"""
from __future__ import annotations

import sys
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]                      # THz-Unified-Optimizer
sys.path.insert(0, str(REPO))

from unified_optimizer import config, model_blanco  # noqa: E402  (read-only use)

sys.path.insert(0, str(REPO / "research" / "experiments"))
import model_core as _core                      # noqa: E402  (общее ядро, зона B, read-only)

# --------------------------------------------------------------------------
# Blanco t_perp/t_par with a cache (PLAN sec.7.2: Blanco is recomputed on every
# lmfit iteration; two WGPs would double that cost).  Key = (freqs, P, D, N, drude).
# --------------------------------------------------------------------------
_T_CACHE: dict = {}
_CACHE_STATS = {"hits": 0, "misses": 0}
_CACHE_MAX = 4096          # FIFO eviction; a fit sweeps D so keys keep growing


def blanco_t(freqs_thz, p_um: float, d_um: float, N: int = 15,
             use_drude: bool = True):
    """Frequency arrays (t_perp, t_par) for one wire grid.

    p_um, d_um: period and effective wire diameter in MICROMETRES.
    Delegates to the shared vectorised kernel `model_core.blanco_t`, which takes
    METRES (handoff B->A, 2026-07-28).  The local FIFO cache below is kept only
    as the legacy reference path (`_blanco_t_legacy`) for the bit-comparison in
    `--verify-core`; the kernel has its own LRU cache.

    NOTE ON UNITS (measured, see `verify_core`): Blanco depends on p and d only
    through `pol = p/lambda` and `d/p`.  Converting um -> m is NOT bit-neutral:
    `p_um*1e-6 / (c/f)` and `p_um / (c/f*1e6)` differ in the last ULP for ~1/3 of
    the frequency grid (1e-6 is not exactly representable in binary).  The
    resulting deviation in t_perp/t_par is at machine-epsilon level and is
    reported numerically by `--verify-core` rather than assumed to be zero.
    """
    t_perp, t_par = _core.blanco_t(np.ascontiguousarray(freqs_thz, dtype=float),
                                   p_um * 1e-6, d_um * 1e-6, N, use_drude)
    return t_perp, t_par


def _blanco_t_legacy(freqs_thz, p_um: float, d_um: float, N: int = 15,
                     use_drude: bool = True):
    """ЛЕГАСИ-путь: скалярный цикл по частотам + собственный FIFO-кэш.

    Оставлен ТОЛЬКО как эталон для бит-сверки (`--verify-core`). В расчётах не
    используется. Не удалять без пересогласования с B: это единственная
    независимая реализация, против которой проверяется общее ядро на этом пути.
    """
    freqs_thz = np.asarray(freqs_thz, dtype=float)
    key = (freqs_thz.tobytes(), freqs_thz.shape, round(float(p_um), 12),
           round(float(d_um), 12), int(N), bool(use_drude))
    hit = _T_CACHE.get(key)
    if hit is not None:
        _CACHE_STATS["hits"] += 1
        return hit
    _CACHE_STATS["misses"] += 1

    d_over_p = d_um / p_um
    tp = np.empty(freqs_thz.shape, dtype=complex)
    ta = np.empty(freqs_thz.shape, dtype=complex)
    flat_f = freqs_thz.ravel()
    flat_tp = tp.ravel()
    flat_ta = ta.ravel()
    for i, f in enumerate(flat_f):
        # p/lambda with both in the same units (um): lambda_um = c/f
        lam_um = model_blanco.C_LIGHT / (f * 1e12) * 1e6
        pol = p_um / lam_um
        flat_tp[i] = model_blanco.compute_t_perp(pol, d_over_p, N, freq_thz=f,
                                                 use_drude=use_drude)
        flat_ta[i] = model_blanco.compute_t_par(pol, d_over_p, N, freq_thz=f,
                                                use_drude=use_drude)
    out = (tp, ta)
    if len(_T_CACHE) >= _CACHE_MAX:
        _T_CACHE.pop(next(iter(_T_CACHE)))
    _T_CACHE[key] = out
    return out


def verify_core(verbose=True):
    """Сверка ЛЕГАСИ-пути против общего ядра (хэндофф B->A, критерий готовности).

    Критерий B был «бит-в-бит, допуск 0». Он НЕ достижим и не должен быть
    достижим: переход мкм -> м меняет порядок операций с плавающей точкой
    (1e-6 не представимо двоично точно). Поэтому здесь измеряется и печатается
    ФАКТИЧЕСКОЕ расхождение в ULP и в относительной мере, а вердикт выносится
    по порогу машинной точности, а не по нулю.
    """
    geoms = [("356att", 15.5, 11.0), ("series (паспорт)", 16.0, 11.0),
             ("specac", 24.9, 14.0), ("test_grid_40_20", 38.8, 18.0),
             ("purewave", 25.5, 10.6)]
    freqs = np.arange(config.F_MIN, getattr(config, "F_MAX_2D", config.F_MAX)
                      + 1e-12, 0.005)
    rows, worst = [], 0.0
    for name, p_um, d_um in geoms:
        for drude in (True, False):
            a_tp, a_ta = _blanco_t_legacy(freqs, p_um, d_um, 15, drude)
            b_tp, b_ta = blanco_t(freqs, p_um, d_um, 15, drude)
            rel = max(_relmax(a_tp, b_tp), _relmax(a_ta, b_ta))
            ulp = max(_ulpmax(a_tp, b_tp), _ulpmax(a_ta, b_ta))
            exact = bool(np.array_equal(a_tp, b_tp) and np.array_equal(a_ta, b_ta))
            worst = max(worst, rel)
            rows.append({"geometry": name, "p_um": p_um, "d_um": d_um,
                         "drude": drude, "bit_exact": exact,
                         "max_rel_dev": rel, "max_ulp": ulp})
            if verbose:
                print(f"  {name:<18} P={p_um:<5} D={d_um:<5} Друде={str(drude):<5} "
                      f"бит-в-бит: {'ДА' if exact else 'нет'}  "
                      f"max отн. = {rel:.2e}  max ULP = {ulp}")
    ok = worst < 1e-13
    if verbose:
        print(f"\n  Худшее относительное расхождение: {worst:.3e} "
              f"(порог машинной точности 1e-13) -> {'ПРИНЯТО' if ok else 'ПРОВАЛ'}")
        print(f"  Бит-в-бит совпадений: {sum(r['bit_exact'] for r in rows)}/{len(rows)}"
              f" — расхождение вносит перевод мкм->м, а не алгоритм ядра.")
    return {"rows": rows, "worst_rel_dev": worst, "accepted": ok,
            "criterion": "max отн. расхождение < 1e-13 (порог машинной точности); "
                         "исходный критерий B «допуск 0» недостижим при смене единиц"}


def _relmax(a, b):
    den = np.maximum(np.abs(a), 1e-300)
    return float(np.max(np.abs(a - b) / den))


def _ulpmax(a, b):
    """Максимальное расхождение в ULP по вещественной и мнимой частям."""
    out = 0
    for x, y in ((a.real, b.real), (a.imag, b.imag)):
        d = np.abs(x - y)
        sp = np.abs(np.spacing(x))
        out = max(out, int(np.max(np.where(sp > 0, d / sp, 0.0))))
    return out


def cache_stats() -> dict:
    """Статистика кэша. С 2026-07-29 расчёт идёт через общее ядро, поэтому
    возвращается ЕГО статистика; локальный FIFO обслуживает только легаси-путь."""
    return {"core": _core.cache_stats(), "legacy_fifo": dict(_CACHE_STATS)}


def clear_cache() -> None:
    _T_CACHE.clear()
    _CACHE_STATS.update(hits=0, misses=0)
    _core.clear_cache()


# --------------------------------------------------------------------------
# One grid = one element with its OWN geometry and OWN dressings
# --------------------------------------------------------------------------
@dataclass
class WGP:
    """One wire-grid polarizer: own geometry + own optional dressings.

    p_um, d_um   : Blanco geometry (d_um is the EFFECTIVE diameter that the fit sees)
    t_override   : (t_perp, t_par) — bypass Blanco entirely (scalars or arrays).
                   Used for ideal-limit validation and for pinning a known element.
    eta0, eta_exp, delta0, tau_leak_ps : optional extra leakage on the parallel
                   channel, exactly the M5/HN8 form
                       t_par -> t_par + eta0*nu^eta_exp * exp(+j(delta0 + 2 pi nu tau_leak))
                   but attached to a PHYSICAL element instead of to the whole system.
    tau_par_ps   : extra group delay of this element's parallel channel; sign
                   convention matches the production model / HN8, i.e. the whole
                   parallel channel is multiplied by exp(-j 2 pi nu tau_par_ps).
    """
    p_um: float = 15.5
    d_um: float = 11.0
    N: int = 15
    use_drude: bool = True
    t_override: tuple | None = None
    eta0: float = 0.0
    eta_exp: float = 1.0
    delta0: float = 0.0
    tau_leak_ps: float = 0.0
    tau_par_ps: float = 0.0

    def ab(self, freqs_thz):
        """Dressed (a, b) = (perpendicular, parallel) amplitude transmissions."""
        freqs_thz = np.asarray(freqs_thz, dtype=float)
        if self.t_override is not None:
            a = np.broadcast_to(np.asarray(self.t_override[0], dtype=complex),
                                freqs_thz.shape).copy()
            b = np.broadcast_to(np.asarray(self.t_override[1], dtype=complex),
                                freqs_thz.shape).copy()
        else:
            a, b = blanco_t(freqs_thz, self.p_um, self.d_um, self.N, self.use_drude)
            a, b = a.copy(), b.copy()
        if self.eta0 != 0.0:
            b = b + (self.eta0 * (freqs_thz ** self.eta_exp)
                     * np.exp(1j * (self.delta0
                                    + 2 * np.pi * freqs_thz * self.tau_leak_ps)))
        if self.tau_par_ps != 0.0:
            b = b * np.exp(-1j * 2 * np.pi * freqs_thz * self.tau_par_ps)
        return a, b


def jones(phi_rad, a, b):
    """Jones matrices of grids with transmission axis at phi_rad.

    phi_rad : broadcastable to shape S (angles)
    a, b    : shape (nf,) or scalar (frequency axis)
    returns : complex array of shape S + (nf, 2, 2)
    """
    phi = np.asarray(phi_rad, dtype=float)[..., None]        # S + (1,)
    a = np.asarray(a, dtype=complex)
    b = np.asarray(b, dtype=complex)
    c, s = np.cos(phi), np.sin(phi)
    cc, ss, cs = c * c, s * s, c * s
    j00 = a * cc + b * ss
    j11 = a * ss + b * cc
    j01 = (a - b) * cs
    j00, j11, j01 = np.broadcast_arrays(j00, j11, j01)
    return np.stack([np.stack([j00, j01], axis=-1),
                     np.stack([j01, j11], axis=-1)], axis=-2)


def _unit(angle_deg, shape):
    """Unit Jones vector at angle_deg (degrees).

    Scalar angle -> shape (2,);  array angle -> broadcast to `shape` + (2,).
    """
    ang = np.asarray(angle_deg, dtype=float)
    if ang.ndim:
        ang = np.broadcast_to(ang, shape)
    ang = np.deg2rad(ang)
    return np.stack([np.cos(ang), np.sin(ang)], axis=-1).astype(complex)


# --------------------------------------------------------------------------
# THE forward model
# --------------------------------------------------------------------------
def forward_field(theta1_deg, theta2_deg, freqs_thz, wgp1: WGP,
                  wgp2: WGP | None = None, *, e_in_deg=0.0, det_deg=0.0,
                  loss_factor: float = 0.0, gamma: float = 2.0,
                  tau_ps: float = 0.0, loss_mode: str = "per_element",
                  full_vector: bool = False):
    """E_out(theta1, theta2, nu) = d^T J_WGP2(theta2) J_WGP1(theta1) E_in.

    theta1_deg, theta2_deg : transmission-axis angles, broadcast against each
        other.  Rig use: theta1 = scan angles shape (n,), theta2 = scalar phi2.
        Map use (A4): theta1 shape (n1,1), theta2 shape (1,n2) -> out (n1,n2,nf).
    wgp2 = None  -> WGP2 is the SAME element as WGP1 (identical grids, shared P,D).
    e_in_deg     : emitter polarization angle; may itself be an array
                   broadcastable to the angle shape (used by the Malus test).
    loss_factor, gamma : phenomenological scattering exp(-0.5*(L/4.343)*nu^gamma)
                   (same parameterization as the production model / M5).
    loss_mode    : "per_element" (each grid scatters — physical for two grids)
                   or "global" (applied once, matching the one-WGP M5 convention).
    tau_ps       : common group delay of the whole chain.
    full_vector  : return the 2-component field instead of the detector projection
                   (needed for the energy-conservation check).

    Returns complex ndarray of shape S + (nf,) — or S + (nf, 2) if full_vector.
    """
    freqs_thz = np.asarray(freqs_thz, dtype=float)
    if wgp2 is None:
        wgp2 = wgp1
    a1, b1 = wgp1.ab(freqs_thz)
    a2, b2 = (a1, b1) if wgp2 is wgp1 else wgp2.ab(freqs_thz)

    if loss_factor != 0.0:
        if config.USE_POWER_LAW:
            loss_amp = np.exp(-0.5 * (loss_factor / 4.343) * (freqs_thz ** gamma))
        else:
            loss_amp = np.exp(-0.5 * loss_factor * freqs_thz)
        if loss_mode == "per_element":
            a1, b1 = a1 * loss_amp, b1 * loss_amp
            a2, b2 = a2 * loss_amp, b2 * loss_amp
        elif loss_mode == "global":
            a1, b1 = a1 * loss_amp, b1 * loss_amp
        else:
            raise ValueError(f"loss_mode must be per_element|global, got {loss_mode!r}")

    th1 = np.deg2rad(np.asarray(theta1_deg, dtype=float))
    th2 = np.deg2rad(np.asarray(theta2_deg, dtype=float))
    th1, th2 = np.broadcast_arrays(th1, th2)                  # both -> S

    J1 = jones(th1, a1, b1)                                   # S + (nf,2,2)
    J2 = jones(th2, a2, b2)
    M = np.einsum('...ij,...jk->...ik', J2, J1)               # S + (nf,2,2)

    e_in = _unit(e_in_deg, th1.shape)
    if e_in.ndim > 1:                                          # per-angle E_in
        E = np.einsum('...fij,...j->...fi', M, e_in)
    else:
        E = np.einsum('...ij,j->...i', M, e_in)

    if tau_ps != 0.0:
        E = E * np.exp(-1j * 2 * np.pi * freqs_thz * tau_ps)[..., None]

    if full_vector:
        return E
    d = _unit(det_deg, th1.shape)
    if d.ndim > 1:                                             # per-angle detector
        return np.einsum('...fi,...i->...f', E, d)
    return np.einsum('...i,i->...', E, d)


def compute_grid_2wgp(angles_deg, freqs_thz, p1_um, d1_um, loss_factor,
                      angle_offset, tau_ps, gamma, tau_par1_ps,
                      *, p2_um=None, d2_um=None, tau_par2_ps=None,
                      phi2_deg=0.0, e_in_deg=0.0, det_deg=0.0,
                      eta0=0.0, eta_exp=1.0, delta0=0.0, tau_leak_ps=0.0,
                      eta0_2=None, eta_exp_2=None, delta0_2=None,
                      N=15, use_drude=True, loss_mode="per_element"):
    """lmfit-facing grid builder, signature parallel to `t6_hn8.compute_grid_hn8`.

    Rig geometry: WGP1 scans (angles - angle_offset), WGP2 fixed at phi2_deg.
    p2_um/d2_um = None -> WGP2 is IDENTICAL to WGP1 (shared geometry, the
    zero-extra-parameter hypothesis of PLAN sec.6/A2).  Same for the leakage
    and tau_par of element 2.
    Returns complex (n_angles, n_freqs) — drop-in for the M5 residual.
    """
    freqs_thz = np.asarray(freqs_thz, dtype=float)
    w1 = WGP(p_um=p1_um, d_um=d1_um, N=N, use_drude=use_drude,
             eta0=eta0, eta_exp=eta_exp, delta0=delta0,
             tau_leak_ps=tau_leak_ps, tau_par_ps=tau_par1_ps)
    shared = (p2_um is None and d2_um is None and tau_par2_ps is None
              and eta0_2 is None and eta_exp_2 is None and delta0_2 is None)
    if shared:
        w2 = None
    else:
        w2 = WGP(p_um=p1_um if p2_um is None else p2_um,
                 d_um=d1_um if d2_um is None else d2_um,
                 N=N, use_drude=use_drude,
                 eta0=eta0 if eta0_2 is None else eta0_2,
                 eta_exp=eta_exp if eta_exp_2 is None else eta_exp_2,
                 delta0=delta0 if delta0_2 is None else delta0_2,
                 tau_leak_ps=tau_leak_ps,
                 tau_par_ps=tau_par1_ps if tau_par2_ps is None else tau_par2_ps)
    theta1 = np.asarray(angles_deg, dtype=float) - angle_offset
    return forward_field(theta1, phi2_deg, freqs_thz, w1, w2,
                         e_in_deg=e_in_deg, det_deg=det_deg,
                         loss_factor=loss_factor, gamma=gamma, tau_ps=tau_ps,
                         loss_mode=loss_mode)


def attenuation(theta1_deg, theta2_deg, freqs_thz, wgp1: WGP,
                wgp2: WGP | None = None, **kw):
    """Power attenuation A = |d^T J2(theta2) J1(theta1) E_in|^2 (A4's core)."""
    return np.abs(forward_field(theta1_deg, theta2_deg, freqs_thz,
                                wgp1, wgp2, **kw)) ** 2


# --------------------------------------------------------------------------
# Validation (PLAN sec.8 criterion A1)
# --------------------------------------------------------------------------
def selftest(verbose: bool = True) -> dict:
    """Energy / ideal-limit / reference-match / limit-position checks."""
    res: dict = {}
    ok = True

    freqs = np.linspace(config.F_MIN, config.F_MAX_2D, 40)
    th = np.linspace(0.0, 90.0, 91)

    # --- T1: energy conservation |E_out|^2 <= |E_in|^2 for a real grid pair ---
    worst = 0.0
    for (P, D) in [(15.5, 11.0), (15.5, 4.683), (38.8, 18.0), (24.9, 6.24)]:
        E = forward_field(th, 0.0, freqs, WGP(p_um=P, d_um=D), None,
                          full_vector=True)
        worst = max(worst, float(np.max(np.sum(np.abs(E) ** 2, axis=-1))))
    res["T1_energy_max_out_over_in"] = worst
    res["T1_pass"] = bool(worst <= 1.0 + 1e-9)
    ok &= res["T1_pass"]

    # --- T2: Malus in the ideal limit, general (theta1, theta2) ---
    ideal = WGP(t_override=(1.0, 0.0))
    t1 = np.linspace(-180, 180, 73)[:, None]
    t2 = np.linspace(-180, 180, 73)[None, :]
    A = attenuation(t1, t2, np.array([1.0]), ideal, ideal,
                    e_in_deg=np.broadcast_to(t1, (73, 73)),  # E_in along WGP1 axis
                    det_deg=0.0)
    # detector along x: |d^T e(theta2)|*cos(theta1-theta2) -> compare to the
    # full projection chain, i.e. Malus times the analyser-to-detector cosine
    expect = (np.cos(np.deg2rad(t1 - t2)) * np.cos(np.deg2rad(t2))) ** 2
    err = float(np.max(np.abs(A[..., 0] - expect)))
    res["T2_malus_max_abs_err"] = err
    res["T2_pass"] = bool(err < 1e-12)
    ok &= res["T2_pass"]
    # ... and the pure Malus law with the detector following the analyser
    A2 = attenuation(t1, t2, np.array([1.0]), ideal, ideal,
                     e_in_deg=np.broadcast_to(t1, (73, 73)),
                     det_deg=np.broadcast_to(t2, (73, 73)))
    err2 = float(np.max(np.abs(A2[..., 0] - np.cos(np.deg2rad(t1 - t2)) ** 2)))
    res["T2b_malus_detector_on_analyser_err"] = err2
    res["T2b_pass"] = bool(err2 < 1e-12)
    ok &= res["T2b_pass"]

    # --- T3: rig ideal limit == production one-WGP model ideal limit ---
    E3 = forward_field(th, 0.0, np.array([1.0]), ideal, ideal,
                       e_in_deg=0.0, det_deg=0.0)[:, 0]
    prod_ideal = np.cos(np.deg2rad(th)) ** 2      # cos^2*t_perp + sin^2*t_par, t=(1,0)
    err3 = float(np.max(np.abs(E3 - prod_ideal)))
    res["T3_ideal_vs_production_err"] = err3
    res["T3_pass"] = bool(err3 < 1e-12)
    ok &= res["T3_pass"]

    # --- T4: bit-match model_blanco.transmission_complex_two_polarizers ---
    # reference: M = P @ R(th) @ P @ R(-th), E_in=x, returns E_out[0]
    #         == J_WGP2(phi2=0) @ J_WGP1(th) with identical grids, d = x.
    # Reference passes freq_thz=None -> Z_drude = 0, so use_drude=False here.
    P_ref, D_ref = 15.5, 4.683
    w = WGP(p_um=P_ref, d_um=D_ref, use_drude=False)
    freqs4 = np.array([0.3, 0.7, 1.2])
    mine = forward_field(th[::10], 0.0, freqs4, w, None,
                         e_in_deg=0.0, det_deg=0.0)
    err4 = 0.0
    for i, tt in enumerate(th[::10]):
        for j, f in enumerate(freqs4):
            lam_um = model_blanco.C_LIGHT / (f * 1e12) * 1e6
            ref = model_blanco.transmission_complex_two_polarizers(
                np.deg2rad(tt), P_ref / lam_um, D_ref / P_ref)
            err4 = max(err4, abs(mine[i, j] - ref))
    res["T4_vs_transmission_complex_two_polarizers_max_abs_err"] = float(err4)
    res["T4_pass"] = bool(err4 < 1e-12)
    ok &= res["T4_pass"]

    # --- T5: limit positions — theta=0 (axes aligned) max, theta=90 min ---
    w5 = WGP(p_um=15.5, d_um=4.683)
    U = np.sum(np.abs(forward_field(th, 0.0, freqs, w5, None)) ** 2, axis=1)
    res["T5_argmax_deg"] = float(th[int(np.argmax(U))])
    res["T5_argmin_deg"] = float(th[int(np.argmin(U))])
    res["T5_extinction_ratio"] = float(U.min() / U.max())
    res["T5_pass"] = bool(th[int(np.argmax(U))] == 0.0 and th[int(np.argmin(U))] == 90.0)
    ok &= res["T5_pass"]

    # --- T6: shared WGP2 == explicitly duplicated WGP2 ---
    e_shared = forward_field(th, 0.0, freqs, w5, None)
    e_dup = forward_field(th, 0.0, freqs, w5, WGP(p_um=15.5, d_um=4.683))
    err6 = float(np.max(np.abs(e_shared - e_dup)))
    res["T6_shared_vs_duplicate_err"] = err6
    res["T6_pass"] = bool(err6 == 0.0)
    ok &= res["T6_pass"]

    # --- T7: how much does the missing WGP2 actually change things? ---
    # Same geometry & dressings, one-WGP production grid vs two-WGP grid.
    from unified_optimizer.optimizer_2d import compute_theoretical_grid_2d
    p_m, d_m = 15.5e-6, 4.683e-6
    one = compute_theoretical_grid_2d(th, freqs, p_m, d_m, 0.0, 0.0, 0.0,
                                      gamma=2.0, use_drude=True, tau_par_ps=0.0)
    two = compute_grid_2wgp(th, freqs, 15.5, 4.683, 0.0, 0.0, 0.0, 2.0, 0.0)
    db = 20 * np.log10(np.maximum(np.abs(two), 1e-30) / np.maximum(np.abs(one), 1e-30))
    res["T7_delta_db_at_theta0"] = float(np.mean(db[0, :]))
    res["T7_delta_db_at_theta90"] = float(np.mean(db[-1, :]))
    res["T7_delta_db_max_abs"] = float(np.max(np.abs(db)))
    res["T7_deep_shadow_ratio_two_over_one"] = float(
        np.mean(np.abs(two[-1, :])) / np.mean(np.abs(one[-1, :])))

    # --- T8: THE ALIGNMENT THEOREM (decisive for hypothesis A) -------------
    # If WGP2's transmission axis coincides with the detector axis, then
    #     d^T J_WGP2(phi2) == t_perp2 * d^T          (J2 is diagonal in that basis)
    # so a REAL WGP2 is mathematically indistinguishable from an IDEAL analyser
    # times the scalar t_perp2(nu).  Verified against an explicitly ideal WGP2
    # for several (E_in, phi2=det) combinations, then broken by misalignment.
    w8 = WGP(p_um=15.5, d_um=11.0)
    a8, _ = w8.ab(freqs)
    t8 = {}
    for psi in (0.0, 5.0, 20.0):
        for delta in (0.0, 5.0, 20.0):          # phi2 == det == delta (aligned)
            real2 = forward_field(th, delta, freqs, w8, w8, e_in_deg=psi, det_deg=delta)
            ideal2 = forward_field(th, delta, freqs, w8, ideal, e_in_deg=psi, det_deg=delta)
            r = np.max(np.abs(real2 - a8[None, :] * ideal2))
            t8[f"aligned_psi{psi:g}_delta{delta:g}"] = float(r)
    res["T8_aligned_max_abs_dev_from_tperp_times_ideal"] = max(t8.values())
    res["T8_pass"] = bool(max(t8.values()) < 1e-14)
    ok &= res["T8_pass"]
    # misalignment phi2 != det DOES break it (sanity: the theorem is not vacuous)
    mis = forward_field(th, 10.0, freqs, w8, w8, e_in_deg=0.0, det_deg=0.0)
    mis_i = forward_field(th, 10.0, freqs, w8, ideal, e_in_deg=0.0, det_deg=0.0)
    res["T8_misaligned10deg_dev"] = float(np.max(np.abs(mis - a8[None, :] * mis_i)))
    res["T8_theorem_nonvacuous"] = bool(res["T8_misaligned10deg_dev"] > 1e-6)
    ok &= res["T8_theorem_nonvacuous"]

    # --- T9: can a real/misaligned WGP2 mimic the collapsed D_eff? ---------
    # The M5 fits needed a COLLAPSED D_eff, which makes the deep shadow DEEPER
    # (Blanco t_par shrinks with D).  Target = the extinction U(90)/U(0) that the
    # fitted D_eff produces; baseline = the physical D.  Question: can WGP2
    # non-ideality / misalignment (delta = phi2 - det, psi = emitter tilt) move
    # the physical-D extinction toward the target?  SIGN MATTERS: leakage terms
    # can only make the extinction ratio LARGER (shadow shallower).
    def extinction(w, phi2, det, psi):
        U = np.sum(np.abs(forward_field(th, phi2, freqs, w, w,
                                       e_in_deg=psi, det_deg=det)) ** 2, axis=1)
        return float(U[-1] / U[0])

    w_phys = WGP(p_um=15.5, d_um=11.0)
    w_eff = WGP(p_um=15.5, d_um=4.683)
    target = extinction(w_eff, 0.0, 0.0, 0.0)
    base_phys = extinction(w_phys, 0.0, 0.0, 0.0)
    res["T9_target_extinction_Deff4p683"] = target
    res["T9_extinction_Dphys11_aligned"] = base_phys
    # negative => the fit needs a DEEPER shadow than physical-D Blanco gives
    res["T9_required_extra_depth_db"] = float(10 * np.log10(target / base_phys))
    scan = {}
    for delta in (0.0, 1.0, 2.0, 5.0, 10.0, 15.0, 20.0):
        for psi in (0.0, 1.0, 2.0, 5.0, 10.0, 15.0, 20.0):
            scan[f"d{delta:g}_psi{psi:g}"] = extinction(w_phys, delta, 0.0, psi)
    res["T9_scan_extinction"] = scan
    worst_key = max(scan, key=lambda k: scan[k])
    res["T9_scan_min"] = float(min(scan.values()))
    res["T9_scan_max"] = float(scan[worst_key])
    res["T9_scan_max_at"] = worst_key
    # does ANY (delta, psi) reach the required (deeper) extinction?
    reach = [k for k, v in scan.items() if v <= target]
    res["T9_combos_reaching_target"] = reach
    res["T9_hypothesisA_can_explain_Deff"] = bool(reach)
    res["T9_direction_db_moved_by_worst_misalignment"] = float(
        10 * np.log10(scan[worst_key] / base_phys))

    # --- T10: which WGP2 convention is physically possible? ----------------
    # PLAN sec.9 asks whether WGP2's WIRES or its TRANSMISSION axis lies along d.
    # Wires along d  <=>  transmission axis at 90 <=> the bright channel passes
    # through WGP2's PARALLEL (blocked) channel -> overall throughput ~ t_par.
    U0_axis = float(np.mean(np.abs(forward_field(0.0, 0.0, freqs, w8, w8,
                                                e_in_deg=0.0, det_deg=0.0)) ** 2))
    U0_wires = float(np.mean(np.abs(forward_field(0.0, 90.0, freqs, w8, w8,
                                                 e_in_deg=0.0, det_deg=0.0)) ** 2))
    res["T10_bright_throughput_phi2_axis_along_d"] = U0_axis
    res["T10_bright_throughput_phi2_wires_along_d"] = U0_wires
    res["T10_ratio"] = float(U0_wires / U0_axis)
    # Measured |T(theta=0)| is O(1); a t_par-suppressed chain is excluded.
    res["T10_wires_along_d_excluded"] = bool(U0_wires / U0_axis < 1e-2)

    res["cache"] = cache_stats()
    res["ALL_PASS"] = bool(ok)

    if verbose:
        print("=== A1 two-WGP model validation ===")
        for k, v in res.items():
            if isinstance(v, float):
                print(f"  {k:<52} {v:.6g}")
            else:
                print(f"  {k:<52} {v}")
    return res


if __name__ == "__main__":
    import argparse
    import json

    ap = argparse.ArgumentParser(description="two-WGP forward model (A1)")
    ap.add_argument("--selftest", action="store_true")
    ap.add_argument("--verify-core", action="store_true",
                    help="сверка легаси-пути Бланко против общего ядра (хэндофф B->A)")
    ap.add_argument("--out", default=str(REPO / "research" / "results" / "two_wgp"
                                        / "a1_validation.json"))
    args = ap.parse_args()

    if args.verify_core:
        print("Сверка blanco_t: ЛЕГАСИ (скалярный цикл) против model_core (ядро B)")
        vr = verify_core()
        out = Path(args.out).with_name("a1_verify_core.json")
        out.parent.mkdir(parents=True, exist_ok=True)
        out.write_text(json.dumps(vr, ensure_ascii=False, indent=2), encoding="utf-8")
        print(f"\nartifact -> {out}")
        raise SystemExit(0 if vr["accepted"] else 1)

    r = selftest()
    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(r, ensure_ascii=False, indent=2), encoding="utf-8")
    print(f"\nartifact -> {out}")
    raise SystemExit(0 if r["ALL_PASS"] else 1)
