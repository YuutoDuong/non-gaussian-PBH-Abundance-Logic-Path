"""Primordial-NG corrections to the SIGW spectrum for the Phase 5 spectral shapes.

For each shape (k* = 1, A = 1) and on a grid in x = k/k*: the Gaussian term T2 and the
disconnected NG terms (hybrid, reducible) by quadrature, and the connected ones (C, Z, planar,
nonplanar) by randomised quasi-Monte Carlo with standard errors (ngpbh.sigw_ng). With
Delta^2_g = A s(k/k*) and zeta = zeta_g + F zeta_g^2 + G zeta_g^3,

    Omega_GW,RD = A^2 T2 + A^3 [F^2 (hybrid + C + Z) + 12 G sbar T2] + F^4 A^4 (reducible + planar + nonplanar),

sbar = sigma^2/A. Writes results/sigw_ng_terms.csv.

Usage: python scripts/sigw_ng_terms.py [--workers 4] [--n-log2 16] [--n-rand 8]
"""

import argparse
import csv
import time
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np

from ngpbh.sigw_ng import LnSampler, disconnected_terms, mc_term

from pta_tension import SHAPES

# the NG15 bins for 10^5 <= k* <= 10^10 Mpc^-1 lie in 1.3e-4 <= x <= 180; the Gaussian template
# stops at x = 30 (beyond, the posterior mass is < 1e-7 for every shape)
X_GRID = np.concatenate([np.geomspace(1e-4, 0.1, 24, endpoint=False), np.geomspace(0.1, 10.0, 60, endpoint=False),
                         np.geomspace(10.0, 30.0, 6)])
CONNECTED = ("C", "Z", "planar", "nonplanar")


def shape_job(args):
    name, n_log2, n_rand = args
    sp = SHAPES[name]
    t0 = time.time()
    det = disconnected_terms(sp, X_GRID)
    sbar = LnSampler(sp).sbar
    rows = []
    for i, x in enumerate(X_GRID):
        row = dict(shape=name, x=x, sbar=sbar, T2=det["gaussian"][i], hybrid=det["hybrid"][i],
                   reducible=det["reducible"][i])
        for term in CONNECTED:
            val, err = mc_term(sp, term, x, n_log2=n_log2, n_rand=n_rand)
            row[term], row[term + "_err"] = val, err
        rows.append(row)
    return name, rows, time.time() - t0


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--workers", type=int, default=4)
    ap.add_argument("--n-log2", type=int, default=16)
    ap.add_argument("--n-rand", type=int, default=8)
    ap.add_argument("--out", default="results/sigw_ng_terms.csv")
    args = ap.parse_args()
    jobs = [(name, args.n_log2, args.n_rand) for name in SHAPES]
    rows = []
    with ProcessPoolExecutor(args.workers) as pool:
        for name, r, dt in pool.map(shape_job, jobs):
            rows += r
            peak = max(r, key=lambda x: x["T2"])
            t3 = peak["hybrid"] + peak["C"] + peak["Z"]
            t4 = peak["reducible"] + peak["planar"] + peak["nonplanar"]
            print(f"{name:10s} ({dt:.0f}s): at the T2 peak x = {peak['x']:.3g}, T3/T2 = {t3 / peak['T2']:.3f}, "
                  f"T4/T2 = {t4 / peak['T2']:.3f}", flush=True)
    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    with open(args.out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]))
        w.writeheader()
        for r in rows:
            w.writerow({k: (f"{v:.7e}" if isinstance(v, float) else v) for k, v in r.items()})


if __name__ == "__main__":
    main()
