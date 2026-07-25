"""
micrograph_period.py — geometric ground-truth from microscope images of a
wire-grid polarizer (WGP) plus a stage-micrometer / graticule calibration slide.

Pipeline
--------
1. Calibration: the graticule has a known smallest division (default 10 um =
   0.01 mm). Its fine comb of ticks is vertical -> variation is along x. We take
   the column-mean profile over a vertical band, and recover the tick spacing in
   pixels two independent ways (FFT fundamental + peak-to-peak spacing). The two
   must agree; their spread feeds the calibration uncertainty. Result: um/px.
2. WGP: wires are vertical -> column-mean profile. FFT fundamental gives the
   period; find_peaks gives the per-period spacing distribution (mean, std =
   departure from ideal periodicity, and N periods sampled). Wire diameter D is
   the FWHM of the bright wire stripes.

All lengths are converted to microns with the paired calibration (same
magnification). The WGP period, computed independently at each magnification,
must agree — that cross-magnification consistency is the built-in validation.

Usage
-----
    python micrograph_period.py --diagnostics <img> [<img> ...]
    python micrograph_period.py --calib <cal.bmp> --wgp <wgp.bmp> [--div-um 10]
"""
from __future__ import annotations
import argparse, os
import numpy as np
from PIL import Image, ImageFile
from scipy.signal import find_peaks
ImageFile.LOAD_TRUNCATED_IMAGES = True


def load_gray(path: str) -> np.ndarray:
    im = Image.open(path).convert("RGB")
    return np.asarray(im, dtype=np.float64).mean(axis=2)


def column_profile(gray: np.ndarray, band=(0.30, 0.70)) -> np.ndarray:
    """Mean intensity per column over a central horizontal band of rows."""
    h = gray.shape[0]
    r0, r1 = int(band[0] * h), int(band[1] * h)
    p = gray[r0:r1, :].mean(axis=0)
    return p - p.mean()


def fft_fundamental(profile: np.ndarray, pmin_px=3.0, pmax_px=None):
    """Dominant spatial period (px) of a 1D profile within [pmin,pmax]."""
    n = profile.size
    if pmax_px is None:
        pmax_px = n / 4.0
    w = np.hanning(n)
    spec = np.abs(np.fft.rfft(profile * w))
    freqs = np.fft.rfftfreq(n, d=1.0)          # cycles/px
    per = np.full_like(freqs, np.inf)
    per[1:] = 1.0 / freqs[1:]
    mask = (per >= pmin_px) & (per <= pmax_px)
    idx = np.where(mask)[0]
    k = idx[np.argmax(spec[idx])]
    # parabolic refinement around the peak bin
    if 1 <= k < len(spec) - 1:
        y0, y1, y2 = spec[k - 1], spec[k], spec[k + 1]
        denom = (y0 - 2 * y1 + y2)
        delta = 0.5 * (y0 - y2) / denom if denom != 0 else 0.0
    else:
        delta = 0.0
    f_ref = freqs[k] + delta * (freqs[1] - freqs[0])
    return 1.0 / f_ref, spec, per


def peak_spacings(profile: np.ndarray, approx_px: float):
    """Positions & spacings of bright stripes, min-distance = 0.6*approx."""
    dist = max(2, int(0.6 * approx_px))
    peaks, props = find_peaks(profile, distance=dist,
                              prominence=profile.std() * 0.5)
    spac = np.diff(peaks)
    return peaks, spac


def peak_widths_at(profile: np.ndarray, peaks: np.ndarray, approx_px: float,
                   frac: float = 0.10) -> np.ndarray:
    """Full width (px) of each bright peak at `frac` of its height above the
    local baseline. For ROUND wires in reflection the diameter is the full
    reflecting footprint, so a low frac (~0.10, base width) tracks the wire
    diameter; FWHM (frac=0.5) only spans the narrow specular glint and
    under-reads the diameter by ~2x. Search is capped at 0.6*period per side."""
    widths = []
    p = profile - profile.min()
    half = max(2, int(0.6 * approx_px))
    for pk in peaks:
        pk = int(pk)
        thr = frac * p[pk]
        i = pk
        while i > max(0, pk - half) and p[i] > thr:
            i -= 1
        left = i + (thr - p[i]) / (p[i + 1] - p[i] + 1e-12)
        j = pk
        while j < min(len(p) - 1, pk + half) and p[j] > thr:
            j += 1
        right = j - (thr - p[j]) / (p[j - 1] - p[j] + 1e-12)
        widths.append(right - left)
    return np.asarray(widths)


def in_focus_mask(profile: np.ndarray, peaks: np.ndarray, approx_px: float,
                  q_edge: float = 0.55, q_amp: float = 0.4) -> np.ndarray:
    """Boolean mask of in-focus wires: steep edges (top (1-q_edge) by max |grad|)
    AND bright (top (1-q_amp) by peak height). Use on samples whose wires sit at
    different depths so only sharp wires are used for the diameter."""
    p = profile - profile.min()
    g = np.gradient(p)
    h = int(0.5 * approx_px)
    idx = peaks.astype(int)
    steep = np.array([max(g[max(0, k - h):k].max(initial=0),
                          -g[k:k + h].min(initial=0)) for k in idx])
    amp = p[idx]
    return (steep >= np.quantile(steep, q_edge)) & (amp >= np.quantile(amp, q_amp))


def analyze_calibration(path, div_um=10.0, band=(0.30, 0.70)):
    g = load_gray(path)
    prof = column_profile(g, band)
    per_fft, spec, per = fft_fundamental(prof, pmin_px=4.0)
    peaks, spac = peak_spacings(prof, per_fft)
    per_pk = np.median(spac) if len(spac) else np.nan
    per_pk_std = np.std(spac) if len(spac) else np.nan
    # use FFT (robust) as central value; disagreement -> uncertainty
    px_per_div = per_fft
    umpx = div_um / px_per_div
    disagree = abs(per_fft - per_pk) / per_fft if np.isfinite(per_pk) else np.nan
    return dict(path=path, px_per_div_fft=per_fft, px_per_div_peak=per_pk,
                peak_std_px=per_pk_std, n_ticks=len(peaks), um_per_px=umpx,
                rel_disagree=disagree)


def lattice_fit(peaks_px, approx_px):
    """Assign each wire an integer lattice index and least-squares fit
    position = a + P*index. Handles missing wires (skipped indices) and
    returns the refined period (px) and RMS positional residual (px) =
    disorder relative to an ideal periodic grid."""
    x = np.sort(np.asarray(peaks_px, float))
    n = np.round((x - x[0]) / approx_px).astype(int)
    # drop accidental duplicate indices (split peaks)
    keep = np.concatenate(([True], np.diff(n) > 0))
    x, n = x[keep], n[keep]
    A = np.vstack([np.ones_like(n), n]).T
    (a, P), *_ = np.linalg.lstsq(A, x, rcond=None)
    resid = x - (a + P * n)
    n_missing = (n[-1] + 1) - len(n)          # skipped lattice sites
    return P, float(np.sqrt(np.mean(resid ** 2))), len(n), n_missing


def analyze_wgp(path, um_per_px, band=(0.25, 0.75), diam_frac=0.10,
                in_focus_only=False):
    """Period, wire diameter and lattice disorder from one WGP micrograph.
    diam_frac: height fraction for the base-width diameter (0.10 tracks the
    round-wire reflecting footprint; 0.5 = specular-glint FWHM, under-reads D by
    ~2x). in_focus_only: measure D only on sharp wires (needed when wires are
    non-coplanar, e.g. a bowed/wound grid where one focus captures only some)."""
    g = load_gray(path)
    prof = column_profile(g, band)
    per_fft, spec, per = fft_fundamental(prof, pmin_px=3.0)
    peaks, spac = peak_spacings(prof, per_fft)
    P_px_mean = np.mean(spac) if len(spac) else per_fft
    P_px_std = np.std(spac) if len(spac) else np.nan
    dpeaks = peaks[in_focus_mask(prof, peaks, per_fft)] if in_focus_only else peaks
    widths = peak_widths_at(prof, dpeaks, per_fft, frac=diam_frac)
    med_w = np.median(widths) if len(widths) else np.nan
    wok = widths[(widths > 0.4 * med_w) & (widths < 1.6 * med_w)] if len(widths) else widths
    D_px_mean, D_px_std = (np.mean(wok), np.std(wok)) if len(wok) else (np.nan, np.nan)
    P_lat, rms_px, n_wires, n_missing = lattice_fit(peaks, per_fft)
    return dict(
        path=path, um_per_px=um_per_px, n_periods=len(spac), n_wires=n_wires,
        n_diam_wires=len(wok),
        P_um_fft=per_fft * um_per_px, P_um_lattice=P_lat * um_per_px,
        disorder_um=rms_px * um_per_px, disorder_pct=100 * rms_px / P_lat,
        n_missing=n_missing,
        P_um_mean=P_px_mean * um_per_px, P_um_std=P_px_std * um_per_px,
        D_um_mean=D_px_mean * um_per_px, D_um_std=D_px_std * um_per_px,
        fill_factor=(D_px_mean / P_lat) if P_lat else np.nan,
    )


def diagnostics(paths):
    for pth in paths:
        g = load_gray(pth)
        for band in [(0.30, 0.70), (0.10, 0.45), (0.55, 0.90)]:
            prof = column_profile(g, band)
            per_fft, spec, per = fft_fundamental(prof, pmin_px=4.0)
            peaks, spac = peak_spacings(prof, per_fft)
            med = np.median(spac) if len(spac) else float("nan")
            print(f"{os.path.basename(pth):32s} band{band} "
                  f"FFT={per_fft:7.2f}px  peakmed={med:7.2f}px  Nteeth={len(peaks):3d}")
        print()


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--diagnostics", nargs="+")
    ap.add_argument("--calib")
    ap.add_argument("--wgp")
    ap.add_argument("--div-um", type=float, default=10.0)
    a = ap.parse_args()
    if a.diagnostics:
        diagnostics(a.diagnostics)
    elif a.calib and a.wgp:
        c = analyze_calibration(a.calib, a.div_um)
        w = analyze_wgp(a.wgp, c["um_per_px"])
        print("CAL:", c)
        print("WGP:", w)
