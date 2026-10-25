"""The two branches of the curvaton map behind the threshold sensitivity of R.

For r_dec < 1 the sudden-decay map zeta = F(zeta_g) is symmetric about zeta_g = -2r/3, where the
curvaton field changes sign (rho_sigma ~ sigma^2): troughs of zeta_g below -2r/3 are peaks of
zeta. f_PBH = 1 is split into branch 1 (zeta_g > -2r/3, F' > 0) and branch 2 (zeta_g < -2r/3,
F' < 0) at the amplitudes of threshold_sensitivity.py (curvaton r = 0.9, C_th x 0.95 / 1 / 1.05,
R-map conventions), and R is recomputed on each branch alone. Also checks that widening the
zeta_r window from 16 to 24 sigma_r leaves f_PBH unchanged.
Reads results/threshold_sensitivity_R.csv; writes results/curvaton_branches.csv and
results/curvaton_branches.md.

Usage: python scripts/curvaton_branches.py   (~20 min on 4 cores, after threshold_sensitivity.py)
"""

import csv
import itertools
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np

from ngpbh.abundance import AbundanceModel
from ngpbh.moments import MomentsConfig
from ngpbh.nongaussian import Curvaton, LocalNG
from ngpbh.spectra import LogNormal
from ngpbh.stats import CollapseModel, PeaksConfig, QuadratureConfig
from ngpbh.threshold import MuscoThreshold

R_DEC = 0.9
TEMPLATE = "curvaton r=0.9"
STATISTICS = {"threshold": ("threshold", PeaksConfig()), "peaks_3b": ("peaks", PeaksConfig(constraint="C_l"))}


class CurvatonBranch(LocalNG):
    """The curvaton map restricted to one side of zeta_g = -2r/3."""

    def __init__(self, branch):
        self.full = Curvaton(R_DEC)
        self.branch = branch

    def F(self, z):
        return self.full.F(z)

    def dF(self, z):
        return self.full.dF(z)

    def d2F(self, z):
        return self.full.d2F(z)

    def d3F(self, z):
        return self.full.d3F(z)

    def bounds(self):
        z_turn = -2.0 * R_DEC / 3.0
        return (z_turn, np.inf) if self.branch == 1 else (-np.inf, z_turn)


def solve(model, a_star):
    """A_* on one branch; f_PBH there is below 1 at the full A_*, so A_* bounds it from below."""
    for bracket in ((a_star, 4.0 * a_star), (1e-4, 2.0)):
        try:
            return model.solve_amplitude(1.0, bracket=bracket)
        except RuntimeError:
            continue
    return np.nan


def job(args):
    width, statistic, factor, a_star = args
    spec = LogNormal(1.0, 1e7, width)
    th = MuscoThreshold(transfer=True).evaluate(spec)
    stat, peaks = STATISTICS[statistic]
    model = AbundanceModel(spec, CollapseModel(C_th=factor * th.C_th, K=3.3, gamma=0.36), ng=Curvaton(R_DEC),
                           statistic=stat, moments_cfg=MomentsConfig(transfer=True, x_cut=4.49), peaks=peaks, n_mh=50)
    model.quad = QuadratureConfig(z_span=24.0)
    row = {"width": width, "statistic": statistic, "factor": factor, "A_star": a_star,
           "f_zspan24": model.fpbh(a_star)}
    model.quad = QuadratureConfig()
    for branch in (1, 2):
        model.ng = CurvatonBranch(branch)
        row[f"f_branch{branch}"] = model.fpbh(a_star)
        row[f"A_star_branch{branch}"] = solve(model, a_star)
    return row


def main():
    with open("results/threshold_sensitivity_R.csv") as fh:
        jobs = [(float(r["width"]), r["statistic"], float(r["factor"]), float(r["A_star"])) for r in csv.DictReader(fh)
                if r["template"] == TEMPLATE]
    with ProcessPoolExecutor(4) as pool:
        rows = list(pool.map(job, sorted(jobs, key=lambda j: j[1] != "peaks_3b")))
    rows.sort(key=lambda r: (r["width"], r["factor"], r["statistic"]))
    with open("results/curvaton_branches.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)

    idx = {(r["width"], r["statistic"], r["factor"]): r for r in rows}
    widths = sorted({r["width"] for r in rows})
    factors = sorted({r["factor"] for r in rows})
    md = ["# The two branches of the curvaton map (r = 0.9) and the threshold sensitivity of R", "",
          ("Branch 2: zeta_g < -2r/3 = -0.6, troughs of zeta_g that are peaks of zeta. Shares of f_PBH = 1 at "
           "the A_* of threshold_sensitivity.py; R = A*_peaks/A*_threshold (peaks of C_l) for the full map "
           "and for each branch alone. Lognormal, f_PBH = 1, R-map conventions."), "",
          "| Delta | C_th x | branch-2 share, threshold | branch-2 share, peaks | R | R, branch 1 | R, branch 2 |",
          "|---|---|---|---|---|---|---|"]
    for wd, f in itertools.product(widths, factors):
        t, p = idx[(wd, "threshold", f)], idx[(wd, "peaks_3b", f)]
        cells = [f"{t['f_branch2']:.3f}", f"{p['f_branch2']:.3f}", f"{p['A_star'] / t['A_star']:.4f}"]
        cells += [f"{p[f'A_star_branch{b}'] / t[f'A_star_branch{b}']:.4f}" for b in (1, 2)]
        md.append(f"| {wd:g} | {f:.2f} | " + " | ".join(cells) + " |")

    def spread(branch):
        out = []
        for wd in widths:
            r = [idx[(wd, "peaks_3b", f)][f"A_star_branch{branch}"] / idx[(wd, "threshold", f)][f"A_star_branch{branch}"]
                 for f in factors]
            out.append(max(abs(x / r[factors.index(1.0)] - 1.0) for x in r))
        return max(out)

    split = max(abs(r["f_branch1"] + r["f_branch2"] - 1.0) for r in rows)
    window = max(abs(r["f_zspan24"] - 1.0) for r in rows)
    md += ["", (f"Largest change of R with the threshold on one branch alone: {spread(1):.2%} (branch 1), "
                f"{spread(2):.2%} (branch 2). Checks: |f_1 + f_2 - 1| <= {split:.1e}; zeta_r window 24 "
                f"instead of 16 sigma_r: |f_PBH - 1| <= {window:.1e}."), ""]
    Path("results/curvaton_branches.md").write_text("\n".join(md))
    print("\n".join(md))


if __name__ == "__main__":
    main()
