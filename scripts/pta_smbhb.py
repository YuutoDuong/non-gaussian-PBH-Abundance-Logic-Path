"""NANOGrav 15 significances with an astrophysical SMBHB background added to the SIGW.

The SMBHB signal is a power law in the timing residuals, rho_B^2 = A_B^2 (f/f_yr)^(3-gamma_B) /
(12 pi^2 f^3 T), with NANOGrav's astrophysical prior on (log10 A_B, gamma_B): a 2D Gaussian with
mean (-15.6, 4.7) and covariance [[0.28, -0.0026], [-0.0026, 0.12]] (2306.16219 eq. A1, fitted to
the GWOnly-Ext SMBHB populations). The powers add, rho^2 = rho_SIGW^2 + rho_B^2, and the SMBHB
parameters are marginalised on a grid. The SIGW keeps the uniform priors of pta_tension.py in
log10 k* and log10 A.

sigma_P of the PBH-allowed region (variant base, f_PBH = 1, threshold -> peaks of C_l) is given
without and with the SMBHB, for log10 A in [-3, 0] (the reference) and in [-4, 0]. With the SMBHB,
small A, where the SMBHB alone explains the signal, lies inside the allowed region, so sigma_P
depends on how far down the prior on A reaches; the profile Delta chi^2 (best fit overall against
best fit in the allowed region) does not. Also reported: the Bayes factor of SIGW + SMBHB against
the SMBHB alone, to compare with NANOGrav's.
Writes results/pta_smbhb.csv, results/pta_smbhb.md and the posteriors (ln L, float32) in
results/pta_smbhb_posteriors.npz.

Usage: python scripts/pta_smbhb.py   (~15 min on 4 cores)
"""

import csv
import itertools
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np
import pta_tension as pt
from scipy.interpolate import CubicSpline
from scipy.special import logsumexp

from ngpbh import thermal
from ngpbh.pta import FreeSpectrum, sigma_one_sided
from ngpbh.sigw import SIGWTemplate

MU_B = np.array([-15.6, 4.7])
COV_B = np.array([[0.28, -0.0026], [-0.0026, 0.12]])
LOG10_KSTAR = np.linspace(5.0, 10.0, 251)
LOG10_A = np.linspace(-4.0, 0.0, 401)
LOG10_AB = np.arange(-17.75, -13.45 + 1e-9, 0.025)  # mu -+ 4 sigma
GAMMA_B = np.arange(3.3, 6.1 + 1e-9, 0.1)
STATS = ("threshold", "peaks_3b")
PRIORS = {"[-3,0]": -3.0, "[-4,0]": -4.0}


def sigw_log10_rho(shape, fs):
    """log10 rho_j of the SIGW at A = 1 on the k* grid, shape (n_kstar, n_bins)."""
    f = fs.freqs
    k = f / thermal.HZ_PER_INV_MPC
    base = thermal.H2_OMEGA_R * thermal.gw_dilution(f)
    template = SIGWTemplate(pt.SHAPES[shape].with_amplitude(1.0))
    with np.errstate(divide="ignore"):
        return np.array([fs.log10_rho(base * template(k / 10.0**lk)) for lk in LOG10_KSTAR])


def smbhb_nodes(fs):
    """Grid of (log10 A_B, gamma_B) with normalised prior weights and the SMBHB ln rho_j^2."""
    a, g = np.meshgrid(LOG10_AB, GAMMA_B, indexing="ij")
    d = np.stack([a.ravel() - MU_B[0], g.ravel() - MU_B[1]], axis=1)
    lw = -0.5 * np.einsum("ni,ij,nj->n", d, np.linalg.inv(COV_B), d)
    keep = lw > -8.0  # inside the 4-sigma ellipse
    lw = lw[keep] - logsumexp(lw[keep])
    lr2 = 2.0 * np.log(10.0) * fs.log10_rho_powerlaw(a.ravel()[keep], g.ravel()[keep])
    return lw, lr2


def loglike_grid(r0, fs, lr2_b=None, lw=None):
    """ln L(k*, A) without (lr2_b None) or with the SMBHB marginalised."""
    lr2_s = 2.0 * np.log(10.0) * (LOG10_A[None, :, None] + r0[:, None, :])  # ln rho_S^2
    if lr2_b is None:
        return fs.loglike_rho(lr2_s / (2.0 * np.log(10.0)))
    acc = np.full(lr2_s.shape[:2], -np.inf)
    for w, b in zip(lw, lr2_b):
        ll = fs.loglike_rho(np.logaddexp(lr2_s, b) / (2.0 * np.log(10.0)))
        acc = np.logaddexp(acc, w + ll)
    return acc


def metrics(logl, line, a_min):
    """sigma_P and the profile Delta chi^2 of A < A_*(k*) for the prior log10 A in [a_min, 0]."""
    keep = LOG10_A >= a_min - 1e-9
    lp, la = logl[:, keep], LOG10_A[keep]
    ok = np.isfinite(line[:, 1])
    a_star = CubicSpline(line[ok, 0], line[ok, 1])(LOG10_KSTAR)
    allowed = la[None, :] < a_star[:, None]
    p = np.exp(lp - lp.max())
    mass = p[allowed].sum() / p.sum()
    dchi2 = 2.0 * (lp.max() - lp[allowed].max())
    return sigma_one_sided(mass), dchi2


def shape_job(shape):
    fs = FreeSpectrum(n_bins=14, smooth=pt.SMOOTH)
    lw, lr2_b = smbhb_nodes(fs)
    r0 = sigw_log10_rho(shape, fs)
    return shape, loglike_grid(r0, fs), loglike_grid(r0, fs, lr2_b, lw)


def main():
    lines = pt.load_lines()
    rows, bf = [], {}
    with ProcessPoolExecutor(4) as pool:
        results = list(pool.map(shape_job, pt.SHAPES))
    np.savez_compressed("results/pta_smbhb_posteriors.npz", log10_kstar=LOG10_KSTAR, log10_A=LOG10_A,
                        **{f"{s}|{k}": ll.astype(np.float32) for s, l0, l1 in results for k, ll in (("nobhb", l0), ("bhb", l1))})
    for shape, ll0, ll1 in results:
        # Bayes factor SIGW + SMBHB vs SMBHB alone (uniform prior on the reference [-3, 0] x [5, 10]);
        # at log10 A = -4 the SIGW is negligible, so ll1 there is the SMBHB-only evidence
        keep = LOG10_A >= -3.0
        bf[shape] = float(np.exp(logsumexp(ll1[:, keep]) - np.log(ll1[:, keep].size) - ll1[:, 0].mean()))
        for template in pt.TEMPLATES:
            for st in STATS:
                line = lines.get((shape, template, st, "base", 1.0))
                if line is None:
                    continue
                row = {"shape": shape, "template": template, "statistic": st}
                row["sigma_p_nobhb"], row["dchi2_nobhb"] = metrics(ll0, line, -3.0)
                for name, a_min in PRIORS.items():
                    row[f"sigma_p_bhb{name}"], row[f"dchi2_bhb{name}"] = metrics(ll1, line, a_min)
                rows.append(row)
        print(f"{shape}: BF(SIGW+SMBHB vs SMBHB) = {bf[shape]:.1f}", flush=True)
    with open("results/pta_smbhb.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)

    idx = {(r["shape"], r["template"], r["statistic"]): r for r in rows}
    cut = pt.SIGMA_RELIABLE
    fmt = lambda x: f">{cut:.0f}" if x > cut else f"{x:.2f}"
    md = ["# NANOGrav 15 with an SMBHB background: sigma_P per statistic", "",
          ("SMBHB power law with NANOGrav's astrophysical prior (2306.16219 eq. A1), marginalised; SIGW priors "
           "uniform in log10 k* in [5, 10] and log10 A. Variant base (2306.17149 conventions), f_PBH = 1. "
           "Cells: threshold -> peaks of C_l. Delta chi^2: best fit overall minus best fit in the allowed region."), "",
          "| shape | template | no SMBHB | SMBHB, log10 A in [-3,0] | SMBHB, log10 A in [-4,0] | Delta chi^2 no SMBHB | Delta chi^2 SMBHB |",
          "|---|---|---|---|---|---|---|"]
    for shape, template in itertools.product(pt.SHAPES, pt.TEMPLATES):
        r = [idx.get((shape, template, st)) for st in STATS]
        if None in r:
            continue
        cells = [f"{fmt(r[0][k])} -> {fmt(r[1][k])}" for k in ("sigma_p_nobhb", "sigma_p_bhb[-3,0]", "sigma_p_bhb[-4,0]")]
        cells += [f"{r[0][k]:.1f} -> {r[1][k]:.1f}" for k in ("dchi2_nobhb", "dchi2_bhb[-3,0]")]
        md.append(f"| {shape} | {template} | " + " | ".join(cells) + " |")

    def shift(key):
        return [idx[(s, t, "peaks_3b")][key] - idx[(s, t, "threshold")][key] for s, t in itertools.product(pt.SHAPES, pt.TEMPLATES)
                if (s, t, "peaks_3b") in idx and idx[(s, t, "peaks_3b")]["sigma_p_nobhb"] <= cut]

    md += ["", "Statistics shift threshold -> peaks over the cases with sigma_P <= 4 without the SMBHB:", ""]
    for key, label in (("sigma_p_nobhb", "no SMBHB"), ("sigma_p_bhb[-3,0]", "SMBHB, [-3,0]"),
                       ("sigma_p_bhb[-4,0]", "SMBHB, [-4,0]")):
        s = shift(key)
        md.append(f"* {label}: {min(s):+.2f} to {max(s):+.2f} (median {np.median(s):+.2f})")
    md += ["", "Bayes factor SIGW + SMBHB against the SMBHB alone (SIGW prior [5,10] x [-3,0]):", ""]
    md += [f"* {s}: {v:.1f}" for s, v in bf.items()]
    Path("results/pta_smbhb.md").write_text("\n".join(md) + "\n")
    print("\n".join(md))


if __name__ == "__main__":
    main()
