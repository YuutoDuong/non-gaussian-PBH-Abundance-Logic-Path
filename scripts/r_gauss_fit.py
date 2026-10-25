"""R_Gauss as a table and a fit, and the factorisation of the NG dependence for the BPL.

R = A*_peaks / A*_threshold (peaks of C_l; conventions of scan_amplitudes.py: Musco threshold with
the transfer function, k < 4.49/r, K = 3.3, gamma = 0.36, f_PBH = 1, k* = 1e7/Mpc):
  * Gaussian zeta on a dense grid of lognormal widths, fitted by R_Gauss = a + b ln Delta;
  * Gaussian zeta for BPL(4, beta, 1), beta = 2, 3, 4;
  * the NG templates of the PTA analysis for BPL(4,3,1), whose double ratio R/R_Gauss tests the
    factorisation A*_peaks(NG) ~ R_Gauss A*_threshold(NG) beyond the lognormal.
Writes results/r_gauss_fit.csv and results/r_gauss_fit.md.

Usage: python scripts/r_gauss_fit.py [--workers 4]   (~5 min on 4 cores)
"""

import argparse
import csv
import itertools
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np

from ngpbh.abundance import AbundanceModel
from ngpbh.moments import MomentsConfig
from ngpbh.nongaussian import Curvaton, Logarithmic, Quadratic
from ngpbh.spectra import BrokenPowerLaw, LogNormal
from ngpbh.stats import CollapseModel, PeaksConfig
from ngpbh.threshold import MuscoThreshold

WIDTHS = (0.1, 0.15, 0.2, 0.3, 0.4, 0.5, 0.7, 1.0, 1.4, 2.0)
BETAS = (2.0, 3.0, 4.0)
TEMPLATES = {"gaussian": None, "fNL=-2": Quadratic(-2.0), "fNL=-1": Quadratic(-1.0), "fNL=-0.5": Quadratic(-0.5),
             "fNL=+0.5": Quadratic(0.5), "fNL=+1": Quadratic(1.0), "fNL=+2": Quadratic(2.0),
             "USR mu=5/2": Logarithmic(2.5), "curvaton r=0.9": Curvaton(0.9)}
STATISTICS = {"threshold": ("threshold", PeaksConfig()), "peaks_3b": ("peaks", PeaksConfig(constraint="C_l"))}


def spectrum(name):
    if name.startswith("LN"):
        return LogNormal(1.0, 1e7, float(name.split("=")[1]))
    return BrokenPowerLaw(1.0, 1e7, 4.0, float(name.split(",")[1]), 1.0)


def job(args):
    shape, template, statistic = args
    spec = spectrum(shape)
    th = MuscoThreshold(transfer=True).evaluate(spec)
    stat, peaks = STATISTICS[statistic]
    model = AbundanceModel(spec, CollapseModel(C_th=th.C_th, K=3.3, gamma=0.36), ng=TEMPLATES[template],
                           statistic=stat, moments_cfg=MomentsConfig(transfer=True, x_cut=4.49), peaks=peaks, n_mh=50)
    try:
        a = model.solve_amplitude(1.0, bracket=(1e-4, 2.0))
    except RuntimeError:
        a = np.nan
    return {"shape": shape, "template": template, "statistic": statistic, "C_th": th.C_th, "A_star": a}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--workers", type=int, default=4)
    args = ap.parse_args()
    ln = [f"LN D={w:g}" for w in WIDTHS]
    bpl = [f"BPL(4,{b:g},1)" for b in BETAS]
    jobs = list(itertools.product(ln + bpl, ["gaussian"], STATISTICS))
    jobs += list(itertools.product(["BPL(4,3,1)"], [t for t in TEMPLATES if t != "gaussian"], STATISTICS))
    with ProcessPoolExecutor(args.workers) as pool:
        rows = list(pool.map(job, sorted(jobs, key=lambda j: (j[1] == "gaussian", j[2] != "peaks_3b"))))
    with open("results/r_gauss_fit.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)

    A = {(r["shape"], r["template"], r["statistic"]): r["A_star"] for r in rows}
    R = {(s, t): A[(s, t, "peaks_3b")] / A[(s, t, "threshold")] for s, t, _ in jobs}
    x = np.log(WIDTHS)
    y = np.array([R[(s, "gaussian")] for s in ln])
    b, a = np.polyfit(x, y, 1)
    res = y - (a + b * x)
    md = ["# R_Gauss: table, fit and the factorisation for the BPL", "",
          ("R = A*_peaks/A*_threshold, peaks of C_l, f_PBH = 1, k* = 1e7/Mpc, Musco threshold with the transfer "
           "function, k < 4.49/r, K = 3.3, gamma = 0.36, BBKS normalisation with V = 4 pi r^3/3."), "",
          "| spectrum | R_Gauss | fit | A*_threshold |", "|---|---|---|---|"]
    for s, xi, yi in zip(ln, x, y):
        md.append(f"| lognormal, Delta = {s.split('=')[1]} | {yi:.4f} | {a + b * xi:.4f} | {A[(s, 'gaussian', 'threshold')]:.4e} |")
    for s in bpl:
        md.append(f"| {s} | {R[(s, 'gaussian')]:.4f} | | {A[(s, 'gaussian', 'threshold')]:.4e} |")
    md += ["", (f"Fit for the lognormal, 0.1 <= Delta <= 2: R_Gauss = {a:.4f} + {b:.4f} ln Delta, largest residual "
                f"{np.abs(res).max():.4f}."), "", "Double ratio R/R_Gauss for BPL(4,3,1):", "",
           "| template | R | R/R_Gauss |", "|---|---|---|"]
    rg = R[("BPL(4,3,1)", "gaussian")]
    for t in TEMPLATES:
        md.append(f"| {t} | {R[('BPL(4,3,1)', t)]:.4f} | {R[('BPL(4,3,1)', t)] / rg:.4f} |")
    quad = [R[("BPL(4,3,1)", t)] / rg for t in TEMPLATES if t.startswith("fNL")]
    md += ["", f"Quadratic templates, |f_NL| <= 2: R/R_Gauss = {min(quad):.3f} to {max(quad):.3f}.", ""]
    Path("results/r_gauss_fit.md").write_text("\n".join(md))
    print("\n".join(md))


if __name__ == "__main__":
    main()
