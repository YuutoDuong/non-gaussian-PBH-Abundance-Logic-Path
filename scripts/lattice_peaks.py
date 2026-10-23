"""Phase 3 check on 3D Gaussian random fields: count maxima of C_l = F'(zeta_r) C_g.

zeta_g is drawn on an N^3 periodic grid with a lognormal spectrum; C_g and zeta_r are obtained
by FFT filtering with exactly the kernels of ngpbh.moments (top hat, shell, transfer function,
k <= 4.49 / r).  Local maxima of C_l (26 neighbours) above several thresholds are counted and
compared with the Kac-Rice predictions of Phase 3a (maxima of C_g) and Phase 3b (maxima of
C_l), evaluated with the lattice's own spectral moments (mode sums), so that only grid effects
separate prediction and measurement.  The Gaussian case calibrates those grid effects.

Usage: python scripts/lattice_peaks.py [--n 256] [--r 8] [--realisations 24]
"""

import argparse
import csv
import time
from pathlib import Path

import numpy as np
from scipy.ndimage import maximum_filter

from ngpbh import windows
from ngpbh.moments import CompactionMoments
from ngpbh.nongaussian import Gaussian, Logarithmic, Quadratic
from ngpbh.stats import CollapseModel, PeaksConfig, QuadratureConfig, log_beta

PHI = 2.0 / 3.0
X_CUT = 4.49


def k_grid(n):
    k1 = 2.0 * np.pi * np.fft.fftfreq(n)
    kz = 2.0 * np.pi * np.fft.rfftfreq(n)
    return np.sqrt(k1[:, None, None] ** 2 + k1[None, :, None] ** 2 + kz[None, None, :] ** 2)


def kernels(k, r):
    x = k * r
    keep = (x <= X_CUT) & (k > 0)
    t = windows.transfer_rd(x)
    kc = np.where(keep, 2.0 * PHI / 3.0 * x**2 * windows.tophat(x) * t, 0.0)
    kr = np.where(keep, windows.shell(x) * t, 0.0)
    return kc, kr


def lognormal_p3d(k, amp, k_star, width):
    """2 pi^2 P(k) / k^3 for P = amp / (sqrt(2 pi) width) exp(-ln^2(k/k*) / (2 width^2))."""
    with np.errstate(divide="ignore", invalid="ignore"):
        p = amp / (np.sqrt(2 * np.pi) * width) * np.exp(-np.log(k / k_star) ** 2 / (2 * width**2))
        out = 2.0 * np.pi**2 * p / k**3
    return np.where(k > 0, out, 0.0)


def lattice_moments(n, r, amp, k_star, width):
    """Exact ensemble moments of the lattice fields (rfft half-space weights), in units of r."""
    k = k_grid(n)
    kc, kr = kernels(k, r)
    p = lognormal_p3d(k, amp, k_star, width)
    w = np.full(k.shape, 2.0)
    w[..., 0] = 1.0
    if n % 2 == 0:
        w[..., -1] = 1.0
    wp = w * p / n**3

    def mom(order, a, b):
        return float(np.sum(wp * k ** (2 * order) * a * b)) * r ** (2 * order)

    return CompactionMoments(
        r, PHI,
        sigma_c2=mom(0, kc, kc), sigma_cr2=mom(0, kc, kr), sigma_r2=mom(0, kr, kr),
        s1_2=mom(1, kc, kc), s2_2=mom(2, kc, kc), scr1_2=mom(1, kc, kr),
        s1r_2=mom(1, kr, kr), s2cr_2=mom(2, kc, kr), s2r_2=mom(2, kr, kr),
    )


def realisation(n, r, amp, k_star, width, rng):
    k = k_grid(n)
    kc, kr = kernels(k, r)
    amp_k = np.sqrt(lognormal_p3d(k, amp, k_star, width))
    noise_k = np.fft.rfftn(rng.standard_normal((n, n, n)))
    zeta_k = noise_k * amp_k
    c_g = np.fft.irfftn(zeta_k * kc, s=(n, n, n))
    zeta_r = np.fft.irfftn(zeta_k * kr, s=(n, n, n))
    return c_g, zeta_r


def count_maxima(field, thresholds):
    is_max = maximum_filter(field, size=3, mode="wrap") == field
    vals = field[is_max]
    return np.array([np.sum(vals > t) for t in thresholds])


def predictions(mom, ng, thresholds, constraint):
    """Expected number of maxima of C_l above each threshold, per grid cell."""
    out = []
    for t in thresholds:
        col = CollapseModel(C_th=0.0, K=1.0, gamma=0.0, estimator="linear", C_l_th=t)
        lb = log_beta(mom, col, ng, "peaks", PeaksConfig(constraint=constraint), QuadratureConfig(n_cl=200))
        out.append(np.exp(lb) / (4.0 * np.pi / 3.0 * mom.r**3))
    return np.array(out)


CASES = {
    "gaussian": Gaussian(),
    "fNL=-2": Quadratic(-2.0),
    "fNL=+2": Quadratic(2.0),
    "mu=5/6": Logarithmic(5.0 / 6.0),
}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--n", type=int, default=256)
    ap.add_argument("--r", type=float, default=8.0, help="smoothing radius in grid cells")
    ap.add_argument("--width", type=float, default=1.0)
    ap.add_argument("--sigma-c", type=float, default=0.15)
    ap.add_argument("--realisations", type=int, default=24)
    ap.add_argument("--seed", type=int, default=20261023)
    ap.add_argument("--out", default="results/phase3b_lattice.csv")
    args = ap.parse_args()

    n, r, width = args.n, args.r, args.width
    k_star = 2.74 / r
    unit = lattice_moments(n, r, 1.0, k_star, width)
    amp = args.sigma_c**2 / unit.sigma_c2
    mom = lattice_moments(n, r, amp, k_star, width)
    nus = np.array([1.5, 2.0, 2.5, 3.0])
    thresholds = nus * mom.sigma_c
    print(f"N={n} r={r} Delta={width} A={amp:.4g} sigma_c={mom.sigma_c:.4f} sigma_r={mom.sigma_r:.4f} "
          f"rho={mom.gamma_cr:.3f} L/r={n / r:.1f}", flush=True)

    rng = np.random.default_rng(args.seed)
    counts = {name: [] for name in CASES}
    t0 = time.time()
    for i in range(args.realisations):
        c_g, zeta_r = realisation(n, r, amp, k_star, width, rng)
        for name, ng in CASES.items():
            counts[name].append(count_maxima(ng.dF(zeta_r) * c_g, thresholds))
        print(f"  realisation {i + 1}/{args.realisations} ({time.time() - t0:.0f}s)", flush=True)

    rows = []
    cells = float(n) ** 3
    for name, ng in CASES.items():
        c = np.array(counts[name], dtype=float)
        mean, err = c.mean(axis=0) / cells, c.std(axis=0, ddof=1) / np.sqrt(len(c)) / cells
        ng_arg = None if name == "gaussian" else ng
        pa = predictions(mom, ng_arg, thresholds, "C_g")
        pb = predictions(mom, ng_arg, thresholds, "C_l")
        for j, nu in enumerate(nus):
            row = dict(case=name, nu=nu, threshold=thresholds[j], per_realisation=c[:, j].mean(),
                       n_lattice=mean[j], err=err[j], n_3a=pa[j], n_3b=pb[j],
                       lattice_over_3a=mean[j] / pa[j], lattice_over_3b=mean[j] / pb[j], err_rel=err[j] / mean[j])
            rows.append(row)
            print(f"{name:9s} nu={nu:.1f}: lattice {mean[j]:.4e} +- {err[j] / mean[j]:.1%}  "
                  f"3a {pa[j]:.4e} (ratio {mean[j] / pa[j]:.3f})  3b {pb[j]:.4e} (ratio {mean[j] / pb[j]:.3f})",
                  flush=True)
    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    with open(args.out, "w", newline="") as fh:
        wr = csv.DictWriter(fh, fieldnames=list(rows[0]))
        wr.writeheader()
        wr.writerows(rows)


if __name__ == "__main__":
    main()
