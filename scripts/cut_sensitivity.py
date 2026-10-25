"""How much does R = A_*^peaks / A_*^threshold depend on the UV cut k r < x_cut of the moments?

Both statistics get the same cut, as in scan_amplitudes.py (peaks: maxima of C_l, Phase 3b).
With the transfer function the moments converge without a cut for these spectra, so "none" is
included. Writes results/r_cut_sensitivity.csv and prints R and R/R_Gauss per cut.

Usage: python scripts/cut_sensitivity.py [--workers 4]   (~8 min on 4 cores)
"""

import argparse
import csv
import itertools
from concurrent.futures import ProcessPoolExecutor

from scan_amplitudes import solve

CUTS = (3.5, 4.49, 6.0, 7.725, None)  # 4.49 and 7.725: first two zeros of the top-hat window
WIDTHS = (0.5, 1.0, 2.0)
TEMPLATES = ("gaussian", "fNL=-2", "fNL=+2")
STATISTICS = ("threshold", "peaks_3b")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--workers", type=int, default=4)
    ap.add_argument("--out", default="results/r_cut_sensitivity.csv")
    args = ap.parse_args()

    jobs = [(w, t, s, c, 1e7, 1.0) for c, w, t, s in itertools.product(CUTS, WIDTHS, TEMPLATES, STATISTICS)]
    jobs.sort(key=lambda j: j[2] != "peaks_3b")  # slowest first
    with ProcessPoolExecutor(args.workers) as pool:
        rows = list(pool.map(solve, jobs))
    with open(args.out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)

    A = {(r["x_cut"], r["width"], r["template"], r["statistic"]): r["A_star"] for r in rows}
    cuts = [c if c is not None else "none" for c in CUTS]
    R = {(c, w, t): A[(c, w, t, "peaks_3b")] / A[(c, w, t, "threshold")] for c in cuts for w in WIDTHS for t in TEMPLATES}
    print("R per cut:", "  ".join(f"{c!s:>6}" for c in cuts))
    for w, t in itertools.product(WIDTHS, TEMPLATES):
        print(f"Delta={w:3.1f} {t:9s}", "  ".join(f"{R[(c, w, t)]:.4f}" for c in cuts))
    dr = max(abs(R[(c, w, t)] / R[(4.49, w, t)] - 1) for c in cuts for w in WIDTHS for t in TEMPLATES)
    ddr = max(abs(R[(c, w, t)] / R[(c, w, "gaussian")] / (R[(4.49, w, t)] / R[(4.49, w, "gaussian")]) - 1)
              for c in cuts for w in WIDTHS for t in TEMPLATES)
    print(f"largest change from x_cut = 4.49: R {dr:.2%}, R/R_Gauss {ddr:.2%}")


if __name__ == "__main__":
    main()
