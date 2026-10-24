"""Phase 5 validation: the NG15 free-spectrum likelihood and the SIGW posterior against published fits.

1. Power-law (SMBHB) fits to the 14-bin HD free spectrum vs NANOGrav 2306.16213.
2. SIGW lognormal fit, Delta marginalised, vs 2306.17149 Fig. 2 (NANOGrav15 numbers).
3. SIGW fits for the fixed shapes used in Phase 5 (2-sigma intervals, best fit).

Writes results/phase5_validation.md. Needs data/ng15 (python scripts/fetch_ng15.py).
"""

from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np

from ngpbh.pta import FreeSpectrum, log_posterior_sigw
from ngpbh.sigw import SIGWTemplate
from ngpbh.spectra import BrokenPowerLaw, LogNormal

LOG10_KSTAR = np.linspace(5.0, 10.0, 501)
LOG10_A = np.linspace(-3.0, 0.0, 301)  # A <= 1, as 2306.17149's posteriors
DELTAS = np.linspace(0.1, 4.0, 40)


def quantiles(x, p, qs=(0.0228, 0.1587, 0.5, 0.8413, 0.9772)):
    c = np.cumsum(p) / np.sum(p)
    return np.interp(qs, c, x)


def fmt(x, p, digits=2):
    lo2, lo1, med, hi1, hi2 = quantiles(x, p)
    return f"{med:.{digits}f} (+{hi2 - med:.{digits}f}/-{med - lo2:.{digits}f})"


def lognormal_logpost(delta):
    fs = FreeSpectrum(n_bins=14)
    return log_posterior_sigw(SIGWTemplate(LogNormal(1.0, 1.0, delta)), fs, LOG10_KSTAR, LOG10_A)


def main():
    fs = FreeSpectrum(n_bins=14)
    lines = ["# Phase 5 validation", "",
             "Likelihood: NANOGrav 15-yr HD free spectrum, 14 bins, KDE refit (ceffyl data, Zenodo 10495907).",
             "Intervals are median and 2-sigma (2.3%-97.7%) unless stated.", ""]

    # 1. power laws
    la = np.linspace(-16.0, -13.0, 3001)
    lnl = fs.loglike_rho(fs.log10_rho_powerlaw(la, 13.0 / 3.0))
    p = np.exp(lnl - lnl.max())
    q = quantiles(la, p, (0.05, 0.5, 0.95))
    a_med = 10 ** q[1]
    lines += ["## 1. Power law h_c = A (f/yr^-1)^((3-gamma)/2)", "",
              "| fit | this refit (median, 90%) | NANOGrav 2306.16213 (median, 90%) |", "|---|---|---|",
              f"| gamma = 13/3: A | {a_med / 1e-15:.2f} (+{(10 ** q[2] - a_med) / 1e-15:.2f}/-{(a_med - 10 ** q[0]) / 1e-15:.2f}) x 1e-15 "
              "| 2.4 (+0.7/-0.6) x 1e-15 |"]
    g = np.linspace(0.0, 7.0, 701)
    LA, G = np.meshgrid(la[::3], g, indexing="ij")
    lnl2 = fs.loglike_rho(fs.log10_rho_powerlaw(LA, G))
    p2 = np.exp(lnl2 - lnl2.max())
    qg = quantiles(g, p2.sum(0), (0.05, 0.5, 0.95))
    lines += [f"| gamma free (uniform 0-7): gamma | {qg[1]:.2f} (+{qg[2] - qg[1]:.2f}/-{qg[1] - qg[0]:.2f}) | 3.2 (+0.6/-0.6) |", ""]

    # 2. lognormal, Delta marginalised (uniform 0.1-4), log10 k* in [5, 10], log10 A in [-3, 0]
    with ProcessPoolExecutor(4) as pool:
        lps = np.array(list(pool.map(lognormal_logpost, DELTAS)))  # (n_delta, n_k, n_A)
    post = np.exp(lps - lps.max())
    lines += ["## 2. Lognormal SIGW, Delta marginalised (uniform 0.1-4)", "",
              "| parameter | this work | 2306.17149 Fig. 2, NANOGrav15 |", "|---|---|---|",
              f"| log10 k* | {fmt(LOG10_KSTAR, post.sum(axis=(0, 2)))} | 8.14 (+0.77/-0.76) |",
              f"| log10 A | {fmt(LOG10_A, post.sum(axis=(0, 1)))} | -0.56 (+0.53/-0.74) |",
              f"| Delta | {fmt(DELTAS, post.sum(axis=(1, 2)))} | 0.98 (+1.13/-0.74) |", ""]

    # 3. fixed shapes used in Phase 5
    shapes = {"LN Delta=0.5": LogNormal(1.0, 1.0, 0.5), "LN Delta=1": LogNormal(1.0, 1.0, 1.0),
              "LN Delta=2": LogNormal(1.0, 1.0, 2.0), "BPL (4,3,1)": BrokenPowerLaw(1.0, 1.0, 4.0, 3.0, 1.0)}
    lines += ["## 3. Fixed shapes (log10 k* in [5, 10], log10 A in [-3, 0])", "",
              "| shape | log10 k* | log10 A | best fit (log10 k*, log10 A) | max ln L |", "|---|---|---|---|---|"]
    for name, sp in shapes.items():
        lp = log_posterior_sigw(SIGWTemplate(sp), fs, LOG10_KSTAR, LOG10_A)
        pp = np.exp(lp - lp.max())
        i, j = np.unravel_index(lp.argmax(), lp.shape)
        lines.append(f"| {name} | {fmt(LOG10_KSTAR, pp.sum(1))} | {fmt(LOG10_A, pp.sum(0))} | "
                     f"({LOG10_KSTAR[i]:.2f}, {LOG10_A[j]:.2f}) | {lp.max():.2f} |")
    lines += ["", "2306.17149 Fig. 1 (BPL, alpha, beta, gamma marginalised), NANOGrav15: "
              "log10 k* = 7.92 (+0.89/-0.62), log10 A = -0.86 (+0.79/-0.65).", ""]
    Path("results").mkdir(exist_ok=True)
    Path("results/phase5_validation.md").write_text("\n".join(lines))
    print("\n".join(lines))


if __name__ == "__main__":
    main()
