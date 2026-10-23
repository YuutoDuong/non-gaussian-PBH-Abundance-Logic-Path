"""Required amplitude A_* for f_PBH = 1 over spectral width and NG template (Phases 1b, 2, 4 preview).

All statistics use the same threshold (Musco et al. with the transfer function, evaluated per
width), the same top-hat window, transfer function, k < 4.49/r cut and K, gamma, so the ratio
R = A_*^peaks / A_*^threshold measures the statistics only.  "peaks" counts maxima of the
Gaussian C_g (Phase 3a, Gaussian-peak approximation); "peaks_3b" counts maxima of the
non-Gaussian C_l = F'(zeta_r) C_g (Phase 3b, docs/phase3b_derivation.tex).

Usage: python scripts/scan_amplitudes.py [--fpbh 1] [--kstar 1e7] [--workers 4]
"""

import argparse
import csv
import itertools
import time
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np

from ngpbh.abundance import AbundanceModel
from ngpbh.moments import MomentsConfig
from ngpbh.nongaussian import Logarithmic, Quadratic
from ngpbh.spectra import LogNormal
from ngpbh.stats import CollapseModel, PeaksConfig
from ngpbh.threshold import MuscoThreshold

STATISTICS = {  # label -> (statistic, PeaksConfig)
    "threshold": ("threshold", PeaksConfig()),
    "peaks": ("peaks", PeaksConfig(constraint="C_g")),
    "peaks_3b": ("peaks", PeaksConfig(constraint="C_l")),
}

WIDTHS = (0.1, 0.3, 0.5, 1.0, 2.0)
TEMPLATES = {
    "gaussian": None,
    "fNL=-2": Quadratic(-2.0),
    "fNL=-1": Quadratic(-1.0),
    "fNL=-0.5": Quadratic(-0.5),
    "fNL=+0.5": Quadratic(0.5),
    "fNL=+1": Quadratic(1.0),
    "fNL=+2": Quadratic(2.0),
    "USR mu=5/6": Logarithmic(5.0 / 6.0),
    "USR mu=5/2": Logarithmic(2.5),
    "USR mu=10": Logarithmic(10.0),
    "log mu=-5/2": Logarithmic(-2.5),
}


def solve(job):
    width, label, statistic, x_cut, kstar, fpbh = job
    spec = LogNormal(1.0, kstar, width)
    th = MuscoThreshold(transfer=True).evaluate(spec)
    stat, peaks = STATISTICS[statistic]
    model = AbundanceModel(
        spec,
        CollapseModel(C_th=th.C_th, K=3.3, gamma=0.36),
        ng=TEMPLATES[label],
        statistic=stat,
        moments_cfg=MomentsConfig(transfer=True, x_cut=x_cut),
        peaks=peaks,
        n_mh=50,
    )
    t0 = time.time()
    try:
        a = model.solve_amplitude(fpbh, bracket=(1e-4, 2.0))
    except RuntimeError as err:  # e.g. f_PBH never reaches the target
        a = np.nan
        print(f"  {width} {label} {statistic}: {err}", flush=True)
    return dict(
        width=width,
        template=label,
        f_nl=getattr(TEMPLATES[label], "f_nl", 0.0),
        statistic=statistic,
        x_cut=x_cut if x_cut is not None else "none",
        C_th=th.C_th,
        q=th.q,
        A_star=a,
        seconds=round(time.time() - t0, 1),
    )


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--fpbh", type=float, default=1.0)
    ap.add_argument("--kstar", type=float, default=1e7)
    ap.add_argument("--workers", type=int, default=4)
    ap.add_argument("--out", default="results/phase2_amplitudes.csv")
    args = ap.parse_args()

    jobs = [(w, t, s, 4.49, args.kstar, args.fpbh) for w, t, s in itertools.product(WIDTHS, TEMPLATES, STATISTICS)]
    # threshold statistics without the cut, as in 2211.01728 / 2306.17149, for reference
    jobs += [(w, t, "threshold", None, args.kstar, args.fpbh) for w, t in itertools.product(WIDTHS, TEMPLATES)]
    jobs.sort(key=lambda j: j[2] != "peaks_3b")  # slowest first, for load balance
    rows = []
    with ProcessPoolExecutor(args.workers) as pool:
        for row in pool.map(solve, jobs):
            rows.append(row)
            print(f"Delta={row['width']:3.1f} {row['template']:12s} {row['statistic']:9s} cut={row['x_cut']!s:5s} "
                  f"A*={row['A_star']:.5e} ({row['seconds']}s)", flush=True)
    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    with open(args.out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)


if __name__ == "__main__":
    main()
