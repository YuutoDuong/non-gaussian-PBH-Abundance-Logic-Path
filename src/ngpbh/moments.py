"""Spectral integrals for the compaction-function statistics at smoothing scale r.

With C_g = -2 Phi r zeta_g'(r) and zeta_r = zeta_g on the shell of radius r, the Fourier
kernels are C_g,k = (2 Phi / 3) (kr)^2 W(kr) zeta_k and zeta_r,k = W_s(kr) zeta_k, giving
(x = k r, P^T = T^2(k, r) P):

    sigma_c^2   = <C_g^2>             = (4 Phi^2/9) int dlnk x^4 W^2 P^T
    sigma_cr^2  = <C_g zeta_r>        = (2 Phi/3)   int dlnk x^2 W W_s P^T
    sigma_r^2   = <zeta_r^2>          =             int dlnk W_s^2 P^T

(2211.01728 eqs. 50-52, 2402.11033 eqs. 64).  For peaks of the field C_g(x; r) we also need
r^2-scaled spectral moments of C_g and the correlation of zeta_r with -lap C_g:

    s1^2 = (4 Phi^2/9) int dlnk x^6 W^2 P^T     (= r^2 <|grad C_g|^2>)
    s2^2 = (4 Phi^2/9) int dlnk x^8 W^2 P^T     (= r^4 <(lap C_g)^2>)
    scr1^2 = (2 Phi/3) int dlnk x^4 W W_s P^T   (= r^2 <(-lap C_g) zeta_r>)

Peaks of the non-Gaussian C_l = F'(zeta_r) C_g (docs/phase3b_derivation.tex) also need the
gradient and Laplacian moments of zeta_r:

    s1r^2  =           int dlnk x^2 W_s^2 P^T   (= r^2 <|grad zeta_r|^2>)
    s2cr^2 = (2 Phi/3) int dlnk x^6 W W_s P^T   (= r^4 <lap C_g lap zeta_r>)
    s2r^2  =           int dlnk x^4 W_s^2 P^T   (= r^4 <(lap zeta_r)^2>)
"""

from dataclasses import dataclass, replace
from typing import Optional
import warnings

import numpy as np
from scipy.integrate import simpson

from . import windows

MOMENT_NAMES = ("sigma_c2", "sigma_cr2", "sigma_r2", "s1_2", "s2_2", "scr1_2", "s1r_2", "s2cr_2", "s2r_2")


@dataclass(frozen=True)
class MomentsConfig:
    window: str = "tophat"
    transfer: bool = True
    x_cut: Optional[float] = None  # hard cut k <= x_cut / r; pbhstat uses 4.49
    phi: float = 2.0 / 3.0
    points_per_period: int = 24
    shape_points: int = 40
    rel_support: float = 1e-16
    x_max: float = 2.0e3


@dataclass(frozen=True)
class CompactionMoments:
    r: float
    phi: float
    sigma_c2: float
    sigma_cr2: float
    sigma_r2: float
    s1_2: float
    s2_2: float
    scr1_2: float
    s1r_2: float = np.nan
    s2cr_2: float = np.nan
    s2r_2: float = np.nan

    def scaled(self, factor):
        """Moments of A * P from moments of P (all are linear in the amplitude)."""
        return replace(self, **{n: getattr(self, n) * factor for n in MOMENT_NAMES})

    @property
    def sigma_c(self):
        return np.sqrt(self.sigma_c2)

    @property
    def sigma_r(self):
        return np.sqrt(self.sigma_r2)

    @property
    def gamma_cr(self):
        """Correlation coefficient of C_g and zeta_r."""
        return self.sigma_cr2 / np.sqrt(self.sigma_c2 * self.sigma_r2)

    @property
    def gamma_peak(self):
        """BBKS gamma = s1^2 / (s0 s2) of the field C_g(x; r)."""
        return self.s1_2 / np.sqrt(self.sigma_c2 * self.s2_2)

    @property
    def rstar_over_r(self):
        """BBKS R_* = sqrt(3) s1 / s2, in units of r."""
        return np.sqrt(3.0 * self.s1_2 / self.s2_2)

    @property
    def gamma_xz(self):
        """Correlation of x = -lap C_g / s2 with z = zeta_r / sigma_r."""
        return self.scr1_2 / np.sqrt(self.s2_2 * self.sigma_r2)


def _kernels(x, cfg):
    """The integrands of MOMENT_NAMES divided by P(k), on an array x = k r."""
    W = windows.WINDOWS[cfg.window](x)
    Ws = windows.shell(x)
    if cfg.transfer:
        T2 = windows.transfer_rd(x) ** 2
    else:
        T2 = np.ones_like(x)
    if cfg.x_cut is not None:
        keep = x <= cfg.x_cut
        W = W * keep
        Ws = Ws * keep
    c = 4.0 * cfg.phi**2 / 9.0
    cr = 2.0 * cfg.phi / 3.0
    x2 = x * x
    x4 = x2 * x2
    return np.stack(
        [
            c * x4 * W * W * T2,
            cr * x2 * W * Ws * T2,
            Ws * Ws * T2,
            c * x4 * x2 * W * W * T2,
            c * x4 * x4 * W * W * T2,
            cr * x4 * W * Ws * T2,
            x2 * Ws * Ws * T2,
            cr * x4 * x2 * W * Ws * T2,
            x4 * Ws * Ws * T2,
        ]
    )


def _lnk_grids(spectrum, r, cfg):
    """Two grids covering the support: uniform in ln k where the spectral shape limits the
    resolution, uniform in k where window oscillations do (their phase grows like k r)."""
    k_lo, k_hi = spectrum.support(cfg.rel_support)
    x_hi = cfg.x_cut if cfg.x_cut is not None else cfg.x_max
    k_hi = min(k_hi, x_hi / r)
    if k_hi <= k_lo:
        return None, None
    h_shape = min(0.02, spectrum.lnk_scale / cfg.shape_points)
    dk = np.pi / (cfg.points_per_period * r)
    k_t = dk / h_shape
    grid_ln = grid_k = None
    if k_t > k_lo:
        top = min(k_t, k_hi)
        n = max(int(np.ceil(np.log(top / k_lo) / h_shape)), 8)
        n += n % 2  # even number of intervals for Simpson
        grid_ln = np.linspace(np.log(k_lo), np.log(top), n + 1)
    if k_t < k_hi:
        bottom = max(k_t, k_lo)
        n = max(int(np.ceil((k_hi - bottom) / dk)), 8)
        n += n % 2
        grid_k = np.linspace(bottom, k_hi, n + 1)
    return grid_ln, grid_k


def compaction_moments(spectrum, r, cfg=MomentsConfig()):
    """All CompactionMoments of `spectrum` at smoothing scale r [Mpc]."""
    if spectrum.is_delta:
        vals = spectrum.A * _kernels(np.array([spectrum.k_star * r]), cfg)[:, 0]
        return CompactionMoments(r, cfg.phi, *vals)
    if not cfg.transfer and cfg.x_cut is None and spectrum.support(cfg.rel_support)[1] * r > cfg.x_max:
        warnings.warn(
            "no transfer function and no k cut: the top-hat moments are dominated by "
            "sub-horizon modes and are truncated at x_max",
            RuntimeWarning,
        )
    grid_ln, grid_k = _lnk_grids(spectrum, r, cfg)
    total = np.zeros(len(MOMENT_NAMES))
    if grid_ln is not None:
        k = np.exp(grid_ln)
        f = _kernels(k * r, cfg) * spectrum(k)
        total += simpson(f, x=grid_ln, axis=-1)
    if grid_k is not None:
        f = _kernels(grid_k * r, cfg) * spectrum(grid_k) / grid_k
        total += simpson(f, x=grid_k, axis=-1)
    return CompactionMoments(r, cfg.phi, *total)


def moments_on_grid(spectrum, r_values, cfg=MomentsConfig()):
    return [compaction_moments(spectrum, float(r), cfg) for r in np.atleast_1d(r_values)]
