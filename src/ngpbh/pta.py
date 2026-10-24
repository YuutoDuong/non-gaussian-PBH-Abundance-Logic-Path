"""NANOGrav 15-year likelihood from the free-spectrum refit (ceffyl; Lamb, Taylor & van Haasteren
2303.15442), and posterior summaries on a parameter grid.

The data are KDE log-densities of the Hellings-Downs-correlated free-spectrum posteriors, one per
frequency bin f_j = j/T (PTArcade data, Zenodo record 10495907, `ng15_30f_fs{hd}_ceffyl`;
`scripts/fetch_ng15.py` downloads and checks it). A spectrum h^2 Omega_GW(f) enters through the
timing-residual power per bin,

    rho_j^2 = h_c^2(f_j) / (12 pi^2 f_j^3 T) = H_100^2 h^2 Omega_GW(f_j) / (8 pi^4 f_j^5 T),

and ln L = sum_j ln p_j(log10 rho_j) over the lowest n_bins bins (14 as in the NANOGrav analyses
and 2306.17149). As in ceffyl, rho below the KDE grid takes the density at its lower edge (the
data only bound rho from above there) and rho above it has zero likelihood. ceffyl takes the
nearest grid point; here the log density is interpolated linearly, which removes the 0.0014-dex
staircase and otherwise agrees.

The KDEs use a compact (Epanechnikov) kernel with bandwidths down to 0.004 dex, so outside the
sampled range the density is floored at ln p = -36 and in the sparse tails it is a comb of
single-sample spikes between floors (bins 1 and 2, the most constraining ones, have almost no
low-rho tail). That is harmless in the bulk but makes tail probabilities (beyond ~2.5 sigma)
lumpy and cliff-like. `smooth` (in dex) replaces each density by max(KDE, KDE * Gaussian): the
bulk is unchanged, holes are filled and the tails decay smoothly beyond the last samples.
"""

import hashlib
from pathlib import Path

import numpy as np
from scipy.special import erfinv

from .thermal import H100_HZ

DATA_URL = "https://zenodo.org/api/records/10495907/files/ng15_30f_fs%7Bhd%7D_ceffyl.zip/content"
DATA_MD5 = "4e67348394468caa4ac21207df33ccba"
DATA_SUBDIR = "ng15_30f_fs{hd}_ceffyl"
DEFAULT_DIR = Path(__file__).resolve().parents[2] / "data" / "ng15"
F_YR = 1.0 / 3.15576e7


def md5sum(path):
    h = hashlib.md5()
    with open(path, "rb") as fh:
        for block in iter(lambda: fh.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


class FreeSpectrum:
    """KDE representation of a PTA free spectrum and the likelihood of a model h^2 Omega_GW(f)."""

    def __init__(self, directory=None, n_bins=14, smooth=None):
        d = Path(directory) if directory is not None else DEFAULT_DIR / DATA_SUBDIR
        if not d.is_dir():
            raise FileNotFoundError(f"{d} not found; run scripts/fetch_ng15.py first")
        freqs = np.load(d / "freqs.npy", allow_pickle=False)
        self.T = 1.0 / freqs[0]
        self.freqs = freqs[:n_bins]
        self.grid = np.load(d / "log10rhogrid.npy", allow_pickle=False)
        dens = np.load(d / "density.npy", allow_pickle=False)[0, :n_bins]
        self.log_density = np.nan_to_num(dens, nan=-36.0)  # ceffyl's cleaning
        self.n_bins = n_bins
        self._dx = self.grid[1] - self.grid[0]
        self._log_dx = np.log(self._dx)
        self.smooth = smooth
        if smooth:
            half = int(np.ceil(5 * smooth / self._dx))
            kern = np.exp(-0.5 * (np.arange(-half, half + 1) * self._dx / smooth) ** 2)
            kern /= kern.sum()
            blurred = np.array([np.convolve(np.exp(d), kern, mode="same") for d in self.log_density])
            with np.errstate(divide="ignore"):
                self.log_density = np.maximum(self.log_density, np.log(blurred))

    def log10_rho(self, h2omega):
        """log10 rho_j [s] for h^2 Omega_GW given at the bin frequencies (last axis)."""
        f = self.freqs
        return 0.5 * np.log10(H100_HZ**2 * np.asarray(h2omega) / (8.0 * np.pi**4 * f**5 * self.T))

    def log10_rho_powerlaw(self, log10_A, gamma):
        """log10 rho_j for h_c = A (f/f_yr)^((3-gamma)/2) (the SMBHB benchmark, gamma = 13/3)."""
        f = self.freqs
        a = np.asarray(log10_A)[..., None]
        g = np.asarray(gamma)[..., None]
        log10_hc2 = 2 * a + (3.0 - g) * np.log10(f / F_YR)
        return 0.5 * (log10_hc2 - np.log10(12.0 * np.pi**2 * f**3 * self.T))

    def loglike_rho(self, log10_rho):
        """sum_j ln p_j(log10 rho_j); log10_rho has the bins on its last axis."""
        x = (np.asarray(log10_rho, dtype=float) - self.grid[0]) / self._dx
        n = self.grid.size
        above = x > n - 1
        x = np.clip(x, 0.0, n - 1.0 - 1e-9)
        i = x.astype(int)
        frac = x - i
        rows = np.arange(self.n_bins)
        lp = (1.0 - frac) * self.log_density[rows, i] + frac * self.log_density[rows, i + 1]
        out = np.sum(lp + self._log_dx, axis=-1)
        return np.where(np.any(above, axis=-1), -np.inf, out)

    def loglike(self, h2omega):
        return self.loglike_rho(self.log10_rho(h2omega))


def log_posterior_sigw(template, freespec, log10_kstar, log10_A):
    """ln L on a (k*, A) grid for one SIGW shape; uniform priors in log10 k* and log10 A.

    template is a `sigw.SIGWTemplate`; h^2 Omega_GW = h^2 Omega_r c_g(f) A^2 Omega~(k/k*), so for
    each k* the A-dependence is a shift of log10 rho by log10 A.
    """
    from . import thermal

    f = freespec.freqs
    k = f / thermal.HZ_PER_INV_MPC
    base = thermal.H2_OMEGA_R * thermal.gw_dilution(f)
    out = np.empty((len(log10_kstar), len(log10_A)))
    for a, lk in enumerate(log10_kstar):
        om = base * template(k / 10.0**lk)  # h^2 Omega_GW at A = 1
        with np.errstate(divide="ignore"):
            r0 = freespec.log10_rho(om)
        out[a] = freespec.loglike_rho(np.asarray(log10_A)[:, None] + r0[None, :])
    return out


def credible_level(log_post, mask):
    """Smallest highest-posterior-density credible mass whose region reaches the masked cells.

    Cells must have equal prior volume. A value of 0.954 means the masked region lies outside the
    95.4% HPD region and touches it (2306.17149's "excluded at 2 sigma").
    """
    lp = np.asarray(log_post, dtype=float).ravel()
    p = np.exp(lp - lp.max())
    p /= p.sum()
    order = np.argsort(p)[::-1]
    level = np.empty_like(p)
    level[order] = np.cumsum(p[order])
    m = np.asarray(mask, dtype=bool).ravel()
    return float(level[m].min()) if m.any() else 1.0


def posterior_mass(log_post, mask):
    lp = np.asarray(log_post, dtype=float)
    p = np.exp(lp - lp.max())
    return float(p[np.asarray(mask, dtype=bool)].sum() / p.sum())


def sigma_one_sided(p):
    """Gaussian-equivalent significance of a tail probability p (one-sided): 0.0228 -> 2."""
    from scipy.special import ndtri

    return float(-ndtri(np.clip(p, 1e-300, 1.0)))


def sigma_1d(level):
    """Gaussian-equivalent significance of a credible level, two-sided 1D (0.683 -> 1, 0.954 -> 2)."""
    return float(np.sqrt(2.0) * erfinv(np.clip(level, 0.0, 1.0 - 1e-16)))


def sigma_2d(level):
    """Significance for 2D Gaussian contours (0.393 -> 1, 0.865 -> 2, 0.989 -> 3)."""
    return float(np.sqrt(-2.0 * np.log1p(-np.clip(level, 0.0, 1.0 - 1e-16))))
