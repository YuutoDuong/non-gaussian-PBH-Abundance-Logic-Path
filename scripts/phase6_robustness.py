"""Phase 6 audit: how robust are the Phase 5 significances to the analysis choices?

For the cases that lie in the reliable range (sigma_P <= 4 for threshold statistics), recompute
sigma_P (posterior mass of the PBH-allowed region) and sigma_HPD for threshold statistics and
peaks 3b under one change at a time:
  * prior box: log10 k* in [6, 9] (2306.17149's plot range) instead of [5, 10]; log10 A down to -4;
  * likelihood: 5 or 30 frequency bins instead of 14; KDE tail smoothing 0, 0.03, 0.1 dex
    instead of 0.05;
  * posterior grid twice as fine in both directions;
  * peak-theory normalisation convention: the peaks line moved by log10 of
    R_convention / R_BBKS (the double ratio cancels conventions to < 0.01, so the shift is the
    same for every NG template): Gaussian peak volume (2 pi)^(3/2) r^3 (x0.955), PS factor 2
    (x1.03), pbhstat v0.0.6 code (x1.115), and the high-peak formula of 2406.20089 as used in
    2602.12252 (x0.80, indicative).
Base variant (2306.17149 conventions), f_PBH = 1. Writes results/phase6_robustness.md.
"""

from pathlib import Path

import numpy as np
from scipy.interpolate import CubicSpline

from ngpbh.pta import FreeSpectrum, credible_level, log_posterior_sigw, posterior_mass, sigma_1d, sigma_one_sided
from ngpbh.sigw import SIGWTemplate

from pta_tension import SHAPES, load_lines

CASES_SHAPES = ["BPL(4,3,1)", "LN D=0.5", "LN D=1"]
CASES_TEMPLATES = ["gaussian", "fNL=-2", "fNL=-1", "fNL=-0.5", "fNL=+0.5", "curvaton r=0.9", "USR mu=5/2"]
LK = np.linspace(5.0, 10.0, 501)
LA = np.linspace(-4.0, 0.0, 801)
CONVENTIONS = {"peak volume (2 pi)^3/2 r^3": 0.955, "PS factor 2": 1.03, "pbhstat v0.0.6 code": 1.115,
               "high-peak formula (2406.20089), indicative": 0.80}


def metrics(lk, la, logpost, line, shift=0.0):
    a_star = CubicSpline(line[:, 0], line[:, 1] + shift)(lk)
    allowed = la[None, :] < a_star[:, None]
    return sigma_1d(credible_level(logpost, allowed)), sigma_one_sided(posterior_mass(logpost, allowed))


def posterior(shape, n_bins=14, smooth=0.05, lk=LK, la=LA):
    return log_posterior_sigw(SIGWTemplate(SHAPES[shape]), FreeSpectrum(n_bins=n_bins, smooth=smooth), lk, la)


def box(lk, la, logpost, k_range=(5.0, 10.0), a_min=-3.0):
    ik = (lk >= k_range[0] - 1e-9) & (lk <= k_range[1] + 1e-9)
    ia = la >= a_min - 1e-9
    return lk[ik], la[ia], logpost[np.ix_(ik, ia)]


def main():
    lines = load_lines()
    rows = {}  # variation -> {(shape, template): (sig_thr, sig_3b)} with sig = (hpd, p)
    for shape in CASES_SHAPES:
        ref_full = posterior(shape)
        variations = {
            "reference (14 bins, smooth 0.05, k* in [5,10], A in [-3,0])": box(LK, LA, ref_full),
            "prior log10 k* in [6, 9]": box(LK, LA, ref_full, k_range=(6.0, 9.0)),
            "prior log10 A in [-4, 0]": box(LK, LA, ref_full, a_min=-4.0),
            "5 frequency bins": box(LK, LA, posterior(shape, n_bins=5)),
            "30 frequency bins": box(LK, LA, posterior(shape, n_bins=30)),
            "KDE tails as distributed (no smoothing)": box(LK, LA, posterior(shape, smooth=None)),
            "KDE smoothing 0.03 dex": box(LK, LA, posterior(shape, smooth=0.03)),
            "KDE smoothing 0.10 dex": box(LK, LA, posterior(shape, smooth=0.10)),
        }
        lk2, la2 = np.linspace(5.0, 10.0, 1001), np.linspace(-3.0, 0.0, 1201)
        variations["posterior grid 2x finer"] = (lk2, la2, posterior(shape, lk=lk2, la=la2))
        for name, (lk, la, lp) in variations.items():
            for t in CASES_TEMPLATES:
                thr = metrics(lk, la, lp, lines[(shape, t, "threshold", "base", 1.0)])
                p3b = metrics(lk, la, lp, lines[(shape, t, "peaks_3b", "base", 1.0)])
                rows.setdefault(name, {})[(shape, t)] = (thr, p3b)
        lk, la, lp = variations["reference (14 bins, smooth 0.05, k* in [5,10], A in [-3,0])"]
        for name, ratio in CONVENTIONS.items():
            for t in CASES_TEMPLATES:
                thr = metrics(lk, la, lp, lines[(shape, t, "threshold", "base", 1.0)])
                p3b = metrics(lk, la, lp, lines[(shape, t, "peaks_3b", "base", 1.0)], shift=np.log10(ratio))
                rows.setdefault(f"peaks normalisation: {name}", {})[(shape, t)] = (thr, p3b)
        print(f"done {shape}", flush=True)

    ref_name = next(iter(rows))
    md = ["# Phase 6 audit: robustness of the Phase 5 significances", "",
          "sigma_P (posterior mass of the PBH-allowed region, one-sided Gaussian equivalent), threshold "
          "statistics -> peaks 3b, base variant, f_PBH = 1. One change at a time from the reference. "
          "Last two columns: range over the cases of the statistics shift Delta sigma_P = 3b - threshold, "
          "and the largest change of any sigma_P relative to the reference.", ""]
    cols = [(s, t) for s in CASES_SHAPES for t in ("gaussian", "fNL=-2", "curvaton r=0.9")]
    md.append("| variation | " + " | ".join(f"{s}, {t}" for s, t in cols) + " | Δσ_P range | max |change| |")
    md.append("|---|" + "---|" * len(cols) + "---|---|")
    for name, res in rows.items():
        cells = []
        for c in cols:
            (h0, p0), (h1, p1) = res[c]
            cells.append(f"{p0:.2f} → {p1:.2f}")
        shifts = [res[c][1][1] - res[c][0][1] for c in res if res[c][1][1] <= 4.0]
        change = max(max(abs(res[c][0][1] - rows[ref_name][c][0][1]), abs(res[c][1][1] - rows[ref_name][c][1][1]))
                     for c in res if rows[ref_name][c][1][1] <= 4.0)
        md.append(f"| {name} | " + " | ".join(cells) + f" | {min(shifts):+.2f} to {max(shifts):+.2f} | {change:.2f} |")
    md += ["", "Cases in the ranges: shapes " + ", ".join(CASES_SHAPES) + "; templates " + ", ".join(CASES_TEMPLATES)
           + "; only cases with sigma_P(3b) <= 4 in the reference enter the last two columns.", ""]
    Path("results/phase6_robustness.md").write_text("\n".join(md) + "\n")
    print("\n".join(md))


if __name__ == "__main__":
    main()
