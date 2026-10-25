"""R and the NANOGrav 15 significances when the collapse threshold changes by +-5%.

The paper uses the same Musco et al. threshold for every NG template, although NG also deforms
the typical profile and hence the threshold (Kehagias, Musco & Riotto 1906.07135; Escriva et al.
2202.01028). As in Iovino et al. 2406.20089, C_th is scaled by 0.95 and 1.05, here for the
templates that relieve the tension (f_NL < 0, the curvaton) and for the Gaussian case:
  * R = A*_peaks / A*_threshold on the lognormal widths of the R map (conventions of
    scan_amplitudes.py: threshold with the transfer function, k < 4.49/r);
  * sigma_P of the PBH-allowed region for the lines of pta_tension.py (variant base, f_PBH = 1).
Writes results/threshold_sensitivity_R.csv, results/threshold_sensitivity_lines.csv and
results/threshold_sensitivity.md.

Usage: python scripts/threshold_sensitivity.py [--workers 4]   (~20 min on 4 cores)
"""

import argparse
import csv
import itertools
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np
import pta_tension as pt

from ngpbh.abundance import AbundanceModel
from ngpbh.cosmology import HorizonMass
from ngpbh.moments import MomentsConfig
from ngpbh.nongaussian import Curvaton, Quadratic
from ngpbh.pta import sigma_one_sided
from ngpbh.spectra import LogNormal
from ngpbh.stats import CollapseModel, PeaksConfig
from ngpbh.threshold import MuscoThreshold

FACTORS = (0.95, 1.0, 1.05)
TEMPLATES = {"gaussian": None, "fNL=-2": Quadratic(-2.0), "fNL=-1": Quadratic(-1.0), "fNL=-0.5": Quadratic(-0.5),
             "curvaton r=0.9": Curvaton(0.9)}
WIDTHS = (0.1, 0.3, 0.5, 1.0, 2.0)
STATISTICS = {"threshold": ("threshold", PeaksConfig()), "peaks_3b": ("peaks", PeaksConfig(constraint="C_l"))}
LINE_SHAPES = ("BPL(4,3,1)", "LN D=0.5", "LN D=1")


def r_job(job):
    width, template, statistic, factor = job
    spec = LogNormal(1.0, 1e7, width)
    th = MuscoThreshold(transfer=True).evaluate(spec)
    stat, peaks = STATISTICS[statistic]
    model = AbundanceModel(spec, CollapseModel(C_th=factor * th.C_th, K=3.3, gamma=0.36), ng=TEMPLATES[template],
                           statistic=stat, moments_cfg=MomentsConfig(transfer=True, x_cut=4.49), peaks=peaks, n_mh=50)
    try:
        a = model.solve_amplitude(1.0, bracket=(1e-4, 2.0))
    except RuntimeError:
        a = np.nan
    return {"width": width, "template": template, "statistic": statistic, "factor": factor,
            "C_th": factor * th.C_th, "A_star": a}


def line_job(job):
    shape, template, statistic, factor = job
    spec = pt.SHAPES[shape]
    th = MuscoThreshold(transfer=False).evaluate(spec)  # variant base
    stat, peaks = STATISTICS[statistic]
    model = AbundanceModel(spec, CollapseModel(C_th=factor * th.C_th, K=3.3, gamma=0.36), ng=TEMPLATES[template],
                           statistic=stat, moments_cfg=MomentsConfig(transfer=True, x_cut=4.49),
                           horizon=HorizonMass(convention="fivt23"), peaks=peaks, n_mh=50)
    try:
        a = model.amplitude_curve(10**pt.LOG10_KSTAR_LINE / pt.K_REF, f_targets=(1.0,), bracket=(1e-4, 2.0))[0]
    except RuntimeError:
        a = np.full(pt.LOG10_KSTAR_LINE.size, np.nan)
    return [{"shape": shape, "template": template, "statistic": statistic, "factor": factor,
             "log10_kstar": round(lk, 3), "log10_A_star": np.log10(v)} for lk, v in zip(pt.LOG10_KSTAR_LINE, a)]


def write(rows, path):
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--workers", type=int, default=4)
    args = ap.parse_args()

    r_jobs = list(itertools.product(WIDTHS, TEMPLATES, STATISTICS, FACTORS))
    l_jobs = list(itertools.product(LINE_SHAPES, TEMPLATES, STATISTICS, (0.95, 1.05)))
    with ProcessPoolExecutor(args.workers) as pool:
        r_rows = list(pool.map(r_job, sorted(r_jobs, key=lambda j: j[2] != "peaks_3b")))
        l_rows = [row for rows in pool.map(line_job, sorted(l_jobs, key=lambda j: j[2] != "peaks_3b")) for row in rows]
    write(r_rows, "results/threshold_sensitivity_R.csv")
    write(l_rows, "results/threshold_sensitivity_lines.csv")

    # R and its change with the threshold
    A = {(r["width"], r["template"], r["statistic"], r["factor"]): r["A_star"] for r in r_rows}
    R = {(w, t, f): A[(w, t, "peaks_3b", f)] / A[(w, t, "threshold", f)] for w, t, f in itertools.product(WIDTHS, TEMPLATES, FACTORS)}
    md = ["# R and the NANOGrav 15 significances for a +-5% change of the collapse threshold", "",
          ("C_th scaled by 0.95 / 1.00 / 1.05 for both statistics. R = A*_peaks/A*_threshold (peaks of C_l), "
           "lognormal, f_PBH = 1, R-map conventions."), "",
          "| template | Delta | R (0.95) | R (1.00) | R (1.05) | A*_thr(0.95)/A*_thr(1) | A*_thr(1.05)/A*_thr(1) |",
          "|---|---|---|---|---|---|---|"]
    for t, w in itertools.product(TEMPLATES, WIDTHS):
        md.append(f"| {t} | {w:g} | {R[(w, t, 0.95)]:.4f} | {R[(w, t, 1.0)]:.4f} | {R[(w, t, 1.05)]:.4f} | "
                  f"{A[(w, t, 'threshold', 0.95)] / A[(w, t, 'threshold', 1.0)]:.3f} | "
                  f"{A[(w, t, 'threshold', 1.05)] / A[(w, t, 'threshold', 1.0)]:.3f} |")
    dr = max(abs(R[(w, t, f)] / R[(w, t, 1.0)] - 1) for w, t, f in R)
    ddr = max(abs(R[(w, t, f)] / R[(w, "gaussian", f)] / (R[(w, t, 1.0)] / R[(w, "gaussian", 1.0)]) - 1) for w, t, f in R)
    md += ["", f"Largest change of R: {dr:.2%}; of R/R_Gauss: {ddr:.2%}.", ""]

    # sigma_P for NANOGrav 15 (variant base), threshold -> peaks
    base = pt.load_lines()
    lines = {}
    for r in l_rows:
        lines.setdefault((r["shape"], r["template"], r["statistic"], r["factor"]), []).append((r["log10_kstar"], r["log10_A_star"]))
    for (s, t, st, v, f), line in base.items():
        if v == "base" and f == 1.0 and s in LINE_SHAPES and t in TEMPLATES and st in STATISTICS:
            lines[(s, t, st, 1.0)] = [tuple(x) for x in line]
    posts = pt.shape_posteriors()
    sig = {}
    for key, line in lines.items():
        _, mass, _ = pt.tension_metrics(posts[key[0]], np.array(sorted(line)))
        sig[key] = sigma_one_sided(mass)
    md += ["## sigma_P (NANOGrav 15, conventions of 2306.17149), threshold -> peaks of C_l", "",
           "| shape | template | C_th x 0.95 | C_th x 1.00 | C_th x 1.05 |", "|---|---|---|---|---|"]
    for s, t in itertools.product(LINE_SHAPES, TEMPLATES):
        cells = [f"{sig[(s, t, 'threshold', f)]:.2f} -> {sig[(s, t, 'peaks_3b', f)]:.2f}" for f in FACTORS]
        md.append(f"| {s} | {t} | " + " | ".join(cells) + " |")
    shift = [sig[(s, t, "peaks_3b", f)] - sig[(s, t, "threshold", f)] for s, t, f in itertools.product(LINE_SHAPES, TEMPLATES, FACTORS)
             if sig[(s, t, "peaks_3b", f)] <= pt.SIGMA_RELIABLE]
    dsig = [sig[(s, t, st, f)] - sig[(s, t, st, 1.0)] for s, t, st, f in itertools.product(LINE_SHAPES, TEMPLATES, STATISTICS, (0.95, 1.05))
            if max(sig[(s, t, st, f)], sig[(s, t, st, 1.0)]) <= pt.SIGMA_RELIABLE]
    md += ["", (f"Statistics shift threshold -> peaks over the cases with sigma_P <= {pt.SIGMA_RELIABLE:.0f}: "
                f"{min(shift):+.2f} to {max(shift):+.2f}. Change of sigma_P from the +-5% threshold: "
                f"{min(dsig):+.2f} to {max(dsig):+.2f}."), ""]
    Path("results/threshold_sensitivity.md").write_text("\n".join(md))
    print("\n".join(md))


if __name__ == "__main__":
    main()
