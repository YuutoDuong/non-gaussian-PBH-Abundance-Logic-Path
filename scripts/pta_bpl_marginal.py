"""NG15 posterior in (log10 k*, log10 A) for the broken power law with alpha, beta, gamma marginalised.

This is the posterior 2306.17149 overlay on their Fig. 3 (from their Fig. 1). The shape
parameters are summed on a grid with trapezoid weights for uniform priors alpha in [1, 8],
beta in [0.25, 10], gamma in [0.25, 10] (their 1D marginals are prior-dominated and consistent
with these ranges; alpha < 1.1 is excluded by the data). Uniform priors log10 k* in [5, 10],
log10 A in [-3, 0], as in scripts/pta_tension.py.

Usage: python scripts/pta_bpl_marginal.py [--workers 4]
Writes results/phase5_posterior_bpl_marginal.npz (log10_kstar, log10_A, logpost).
"""

import argparse
import itertools
from concurrent.futures import ProcessPoolExecutor

import numpy as np
from scipy.special import logsumexp

from ngpbh.pta import FreeSpectrum, log_posterior_sigw
from ngpbh.sigw import SIGWTemplate
from ngpbh.spectra import BrokenPowerLaw

from pta_tension import LOG10_A, LOG10_KSTAR

ALPHA = np.array([1.0, 1.5, 2.5, 4.0, 6.0, 8.0])
BETA = np.array([0.25, 0.75, 1.5, 3.0, 5.0, 7.5, 10.0])
GAMMA = np.array([0.25, 0.75, 1.5, 3.0, 5.0, 7.5, 10.0])
GRID_KW = dict(h_factor=1.0, t_max_cap=1e8)  # Omega to 0.15%; x >= 1e-4 covers the 14 bins for k* <= 1e10


def trapezoid_weights(x):
    w = np.zeros_like(x)
    w[1:] += 0.5 * np.diff(x)
    w[:-1] += 0.5 * np.diff(x)
    return w / w.sum()


def job(abg):
    a, b, g = abg
    tpl = SIGWTemplate(BrokenPowerLaw(1.0, 1.0, a, b, g), x_min=1e-4, n_x=320, grid_kw=GRID_KW)
    return log_posterior_sigw(tpl, FreeSpectrum(n_bins=14), LOG10_KSTAR, LOG10_A)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--workers", type=int, default=4)
    args = ap.parse_args()
    grid = list(itertools.product(ALPHA, BETA, GAMMA))
    w = [np.log(x) for x in itertools.product(trapezoid_weights(ALPHA), trapezoid_weights(BETA), trapezoid_weights(GAMMA))]
    log_w = np.array([sum(x) for x in w])
    with ProcessPoolExecutor(args.workers) as pool:
        lps = np.array(list(pool.map(job, grid, chunksize=4)))
    logpost = logsumexp(lps + log_w[:, None, None], axis=0)
    evid = logsumexp(lps.reshape(len(grid), -1), axis=1) + log_w
    order = np.argsort(evid)[::-1][:5]
    print("highest-evidence shapes (alpha, beta, gamma):", [grid[i] for i in order])
    np.savez_compressed("results/phase5_posterior_bpl_marginal.npz", log10_kstar=LOG10_KSTAR, log10_A=LOG10_A,
                        logpost=logpost)


if __name__ == "__main__":
    main()
