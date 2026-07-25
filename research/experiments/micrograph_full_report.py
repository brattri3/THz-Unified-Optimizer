"""
micrograph_full_report.py — comprehensive cross-method / cross-magnification
characterisation of the WGP micrographs. For every sample and magnification it
measures PERIOD, DIAMETER and PERIODICITY DISORDER by several independent
methods and writes a Markdown report with the full comparison tables.

Methods
  Period    : FFT fundamental | autocorrelation | median NN peak-spacing | lattice-fit slope
  Diameter  : FWHM(50%) glint | base-width@20% | base-width@10% | duty-cycle (Otsu)
  Disorder  : NN-spacing CV% | lattice-residual RMS% | missing/merged-wire fraction
Diameter is measured on in-focus wires at max mag (round-wire footprint = base
width; FWHM only spans the specular glint). Absolute low-mag scale is approximate
(fine graticule comb unresolved) so low-mag % metrics — which are calibration-
independent — are the reliable ones there.
"""
from __future__ import annotations
import numpy as np
from micrograph_period import (load_gray, column_profile, fft_fundamental,
                               peak_spacings, peak_widths_at, in_focus_mask,
                               lattice_fit)

MIC = "research/validation/micrographs/"
# PureWave calibrations (µm/px): high & mid from resolved 10 µm comb; low from
# magnification ratio (fine comb unresolved) -> flagged approximate.
UMPX = {"low": 10 / 15.52 * 2.0, "mid": 10 / 15.52, "high": 10 / 78.77}


def autocorr_period(p, approx_px, pmax=None):
    """Period = lag of the first strong autocorrelation peak. Search around the
    expected period so the monotonic near-zero-lag decay is not mistaken for it."""
    from scipy.signal import find_peaks
    n = p.size
    ac = np.correlate(p, p, "full")[n - 1:]
    ac = ac / ac[0]
    lo = max(2, int(0.5 * approx_px))
    hi = int(min(pmax or n / 4, 3.0 * approx_px))
    peaks, _ = find_peaks(ac[:hi], distance=max(2, int(0.5 * approx_px)))
    peaks = peaks[peaks >= lo]
    if len(peaks) == 0:
        return lo + int(np.argmax(ac[lo:hi]))
    return float(peaks[np.argmax(ac[peaks])])


def duty_cycle(p, umpx, per_px):
    p = p - p.min()
    # Otsu threshold
    hist, edges = np.histogram(p, 256)
    c = hist.cumsum(); m = (hist * (edges[:-1])).cumsum()
    tot = m[-1]; best_t, best_v = 0, -1
    for i in range(1, 256):
        w0 = c[i]; w1 = c[-1] - w0
        if w0 == 0 or w1 == 0:
            continue
        mu0 = m[i] / w0; mu1 = (tot - m[i]) / w1
        v = w0 * w1 * (mu0 - mu1) ** 2
        if v > best_v:
            best_v, best_t = v, edges[i]
    bright = (p > best_t).mean()
    return bright * per_px * umpx          # D = duty * P


def period_methods(p, umpx):
    per_fft, _, _ = fft_fundamental(p, pmin_px=3.0)
    pk, spac = peak_spacings(p, per_fft)
    spac = np.asarray(spac, float)
    med = np.median(spac) if len(spac) else np.nan
    inl = spac[(spac > 0.6 * med) & (spac < 1.4 * med)] if len(spac) else spac
    P_lat, rms, nw, nmiss = lattice_fit(pk, per_fft)
    return dict(
        fft=per_fft * umpx,
        acf=autocorr_period(p, per_fft) * umpx,
        nn=(np.mean(inl) * umpx) if len(inl) else np.nan,
        lat=P_lat * umpx,
        # disorder
        nn_cv=100 * np.std(inl) / np.mean(inl) if len(inl) else np.nan,
        lat_rms_pct=100 * rms / P_lat,
        missing_pct=100 * nmiss / (nw + nmiss) if (nw + nmiss) else np.nan,
        n_wires=nw, per_px=per_fft, peaks=pk,
    )


def diameter_methods(p, per_px, umpx, in_focus=False):
    pk, _ = peak_spacings(p, per_px)
    if in_focus:
        pk = pk[in_focus_mask(p, pk, per_px)]
    out = {}
    for lab, frac in [("fwhm50", 0.5), ("w20", 0.20), ("w10", 0.10)]:
        w = peak_widths_at(p, pk, per_px, frac=frac)
        w = w[np.isfinite(w) & (w > 0)]
        mw = np.median(w) if len(w) else np.nan
        wok = w[(w > 0.4 * mw) & (w < 1.6 * mw)] if len(w) else w
        out[lab] = np.mean(wok) * umpx if len(wok) else np.nan
    out["duty"] = duty_cycle(p, umpx, per_px)
    out["n"] = len(pk)
    return out


SAMPLES = {
    "PureWave": {
        "low":  [MIC + "purewave/purewave__wgp__lowmag.bmp"],
        "mid":  [MIC + "purewave/purewave__wgp__midmag.bmp"],
        "high": [MIC + f"purewave/purewave__wgp__highmag_0{i}.bmp" for i in (1, 2, 3)],
        "focus": False,
    },
    "Grid_D20_P40": {
        "low":  [MIC + "wgp_grid_D20_P40/wgp_grid_D20_P40__wgp__lowmag.bmp"],
        "mid":  [MIC + "wgp_grid_D20_P40/wgp_grid_D20_P40__wgp__midmag.bmp"],
        "high": [MIC + f"wgp_grid_D20_P40/wgp_grid_D20_P40__wgp__highmag_{f}.bmp"
                 for f in ("focus0um", "focus200um")],
        "focus": True,
    },
}


def run():
    lines = ["# Micrograph comprehensive analysis (auto-generated)\n",
             "Calibration µm/px: high=0.1270, mid=0.6443, low≈1.289 (approx — "
             "low fine comb unresolved; low absolute lengths indicative, % metrics valid).\n"]
    for sample, cfg in SAMPLES.items():
        lines += [f"\n## {sample}\n",
                  "### Period (µm) — 4 methods\n",
                  "| mag | N | FFT | autocorr | NN-median | lattice | spread |",
                  "|-----|---|-----|----------|-----------|---------|--------|"]
        per_store = {}
        for mag in ("low", "mid", "high"):
            vals = [period_methods(column_profile(load_gray(f)), UMPX[mag]) for f in cfg[mag]]
            m = {k: np.nanmean([v[k] for v in vals]) for k in ("fft", "acf", "nn", "lat")}
            per_store[mag] = vals
            nw = int(np.mean([v["n_wires"] for v in vals]))
            spread = np.nanstd([m["fft"], m["acf"], m["nn"], m["lat"]])
            lines.append(f"| {mag} | {nw} | {m['fft']:.2f} | {m['acf']:.2f} | "
                         f"{m['nn']:.2f} | {m['lat']:.2f} | ±{spread:.2f} |")
        lines += ["\n### Diameter (µm) — high mag, in-focus wires; 4 methods\n",
                  "| frame | FWHM(glint) | width@20% | width@10% | duty(Otsu) |",
                  "|-------|-------------|-----------|-----------|------------|"]
        for f in cfg["high"]:
            p = column_profile(load_gray(f))
            per, _, _ = fft_fundamental(p, pmin_px=3.0)
            d = diameter_methods(p, per, UMPX["high"], in_focus=cfg["focus"])
            name = f.split("__")[-1].replace(".bmp", "")
            lines.append(f"| {name} | {d['fwhm50']:.1f} | {d['w20']:.1f} | "
                         f"{d['w10']:.1f} | {d['duty']:.1f} |")
        lines += ["\n### Periodicity disorder — 3 methods\n",
                  "| mag | N | NN-spacing CV% | lattice-RMS% | missing/merged % |",
                  "|-----|---|----------------|--------------|------------------|"]
        for mag in ("low", "mid", "high"):
            vals = per_store[mag]
            cv = np.nanmean([v["nn_cv"] for v in vals])
            rms = np.nanmean([v["lat_rms_pct"] for v in vals])
            ms = np.nanmean([v["missing_pct"] for v in vals])
            nw = int(np.mean([v["n_wires"] for v in vals]))
            lines.append(f"| {mag} | {nw} | {cv:.1f} | {rms:.1f} | {ms:.1f} |")
    report = "\n".join(lines) + "\n"
    out = "research/validation/ANALYSIS_REPORT.md"
    with open(out, "w", encoding="utf-8") as fh:
        fh.write(report)
    print(report.encode("ascii", "replace").decode())
    print("written:", out)


if __name__ == "__main__":
    run()
