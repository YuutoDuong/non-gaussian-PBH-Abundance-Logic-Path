"""Mass fraction beta at a fixed smoothing scale r: threshold statistics and peaks theory.

Both statistics integrate the same bivariate Gaussian P_G(C_g, zeta_r) over the same type-I
domain {C(C_l) > C_th, C_l < 2 Phi}, with C_l = F'(zeta_r) C_g and C = C_l - C_l^2/(4 Phi):

    beta = int dC_g dzeta_r  K (C - C_th)^gamma  P_G(C_g, zeta_r)  w(C_g, zeta_r).

Threshold statistics (2211.01728 eq. 56, 2402.11033 eq. 62) have w = 1.  Peaks theory counts
maxima of the Gaussian field C_g(x; r) inside the volume V = (4 pi/3) r^3 (2406.20089 eq. 3.12):
with nu = C_g / sigma_c, z = zeta_r / sigma_r and x = -lap C_g / s2, the BBKS density with one
extra correlated scalar is

    n_pk(nu, z) = sqrt(2 pi) / ((2 pi)^2 R_*^3)  P(nu, z)  G(m, s),
    G(m, s) = int_0^inf dx f(x) N(x; m, s^2),

where m and s^2 are the mean and variance of x given (nu, z).  The gradient and the traceless
Hessian are independent of all three scalars, so BBKS's f(x) is unchanged.  With F' < 0 a
maximum of C_l sits at a minimum of C_g, so G(-m, s) is used there.  Counting peaks of C_g
rather than of the non-Gaussian C_l is the Gaussian-peak approximation (Young 2201.13345); it is
exact for F' = const and is the assumption Phase 3 tests.

Everything is computed as log beta because beta is exponentially small.
"""

from dataclasses import dataclass
from typing import Optional

import numpy as np
from scipy.special import erf, logsumexp

from ._quad import leggauss as _leggauss
from .nongaussian import LocalNG
from .threshold import cl_threshold

LOG_2PI = np.log(2.0 * np.pi)

_F_SERIES = (0.075603752829436973499, -0.047252345518398108437, 0.018793546512999247674)


def bbks_f(x):
    """BBKS (A15) peak-curvature weight; the x^8 series avoids cancellation at small x."""
    x = np.asarray(x, dtype=float)
    out = np.zeros_like(x)
    small = (x > 0) & (x < 0.2)
    xs2 = x[small] ** 2
    out[small] = xs2**4 * (_F_SERIES[0] + xs2 * (_F_SERIES[1] + xs2 * _F_SERIES[2]))
    big = x >= 0.2
    xb = x[big]
    a = np.sqrt(2.5)
    out[big] = 0.5 * (xb**3 - 3.0 * xb) * (erf(a * xb) + erf(0.5 * a * xb)) + np.sqrt(
        2.0 / (5.0 * np.pi)
    ) * ((7.75 * xb**2 + 1.6) * np.exp(-0.625 * xb**2) + (0.5 * xb**2 - 1.6) * np.exp(-2.5 * xb**2))
    return out


def log_G(m, s, n_x=200):
    """log of G(m, s) = int_0^inf f(x) N(x; m, s^2) dx, vectorised over m."""
    m = np.atleast_1d(np.asarray(m, dtype=float))
    out = np.full(m.shape, -np.inf)
    lo = np.maximum(m - 12.0 * s, 0.0)
    hi = np.maximum(m + 12.0 * s, lo + 1e-12)
    ok = m + 12.0 * s > 0
    if not np.any(ok):
        return out
    nodes, weights = _leggauss(n_x)
    lo_, hi_ = lo[ok, None], hi[ok, None]
    x = 0.5 * (hi_ - lo_) * nodes + 0.5 * (hi_ + lo_)
    fx = bbks_f(x)
    log_gauss = -0.5 * ((x - m[ok, None]) / s) ** 2 - np.log(np.sqrt(2.0 * np.pi) * s)
    with np.errstate(divide="ignore"):
        terms = np.log(fx) + log_gauss + np.log(0.5 * (hi_ - lo_) * weights)
    out[ok] = logsumexp(terms, axis=1)
    return out


class _LogGTable:
    """log G(m, s) on a dense m grid for fixed s; the peak sum queries it ~1e6 times.

    Above the table f(x) = x^3 - 3x up to terms of order exp(-5x^2/8), so
    G = m^3 - 3m + 3 m s^2 there."""

    def __init__(self, s, m_lo=None, m_hi=30.0, dm=0.01):
        self.s = s
        m_lo = -12.0 * s if m_lo is None else m_lo
        self.m = np.arange(m_lo, m_hi + dm, dm)
        self.vals = log_G(self.m, s)

    def __call__(self, m):
        m = np.asarray(m, dtype=float)
        out = np.interp(m, self.m, self.vals, left=-np.inf, right=np.nan)
        high = m > self.m[-1]
        if np.any(high):
            mh = m[high]
            out[high] = np.log(mh**3 - 3.0 * mh + 3.0 * mh * self.s**2)
        return out


_TABLES = {}


def _table_for(s):
    """Shared tables keyed by s rounded to 1e-4 (d log G / ds is O(1/m^2))."""
    key = round(float(s), 4)
    if key not in _TABLES:
        _TABLES[key] = _LogGTable(key)
    return _TABLES[key]


@dataclass(frozen=True)
class CollapseModel:
    """Critical collapse M/M_H = K (X - X_th)^gamma with estimator X.

    estimator="nonlinear": X = C = C_l - C_l^2/(4 Phi), threshold C_th (2211.01728).
    estimator="linear":    X = C_l, threshold C_l_th (pbhstat's PS and peaks models).
    In both cases only type-I perturbations, C_l < 2 Phi, are integrated.
    """

    C_th: float
    K: float = 3.3
    gamma: float = 0.36
    estimator: str = "nonlinear"
    C_l_th: Optional[float] = None
    phi: float = 2.0 / 3.0

    @property
    def cl_lower(self):
        if self.C_l_th is not None:
            return self.C_l_th
        return float(cl_threshold(self.C_th, self.phi))

    @property
    def cl_upper(self):
        return 2.0 * self.phi

    def log_mass_ratio(self, c_l):
        if self.estimator == "nonlinear":
            excess = (c_l - c_l**2 / (4.0 * self.phi)) - self.C_th
        elif self.estimator == "linear":
            excess = c_l - self.cl_lower
        else:
            raise ValueError("estimator must be 'nonlinear' or 'linear'")
        with np.errstate(divide="ignore", invalid="ignore"):
            return np.log(self.K) + self.gamma * np.log(np.maximum(excess, 0.0))


@dataclass(frozen=True)
class PeaksConfig:
    density: str = "bbks"  # or "high_peak": G -> (gamma nu)^3, as in pbhstat and 2406.20089
    volume: str = "tophat"  # V = 4 pi r^3 / 3; "gaussian": (2 pi)^(3/2) r^3 (Young 2201.13345)

    @property
    def log_volume(self):
        if self.volume == "tophat":
            return np.log(4.0 * np.pi / 3.0)
        if self.volume == "gaussian":
            return 1.5 * LOG_2PI
        raise ValueError("volume must be 'tophat' or 'gaussian'")


@dataclass(frozen=True)
class QuadratureConfig:
    n_cl: int = 160  # Gauss-Legendre nodes on t, C_l = lo + (hi - lo) t^cl_power
    cl_power: float = 3.0
    z_span: float = 16.0  # zeta_r range in units of sigma_r
    dz: float = 0.02  # trapezoid step in z, refined automatically for narrow spectra
    ps_factor: float = 1.0  # extra factor on threshold statistics (pbhstat uses 2)

    def refined(self, factor=2):
        return QuadratureConfig(self.n_cl * factor, self.cl_power, self.z_span, self.dz / factor, self.ps_factor)


def _cl_nodes(collapse, quad):
    t, w = _leggauss(quad.n_cl)
    t = 0.5 * (t + 1.0)
    w = 0.5 * w
    lo, hi = collapse.cl_lower, collapse.cl_upper
    p = quad.cl_power
    c_l = lo + (hi - lo) * t**p
    log_w = np.log((hi - lo) * p * t ** (p - 1.0) * w)
    return c_l, log_w


@dataclass
class PeakGeometry:
    """Conditional statistics of x = -lap C_g / s2 given (nu, z), from CompactionMoments."""

    gamma: float
    gamma_cr: float
    gamma_xz: float
    log_prefactor: float  # log[V sqrt(2 pi) / ((2 pi)^2 R_*^3)], lengths in units of r
    a: float  # E[x | nu, z] = a nu + b z
    b: float
    s: float
    table: Optional[_LogGTable] = None
    _gaussian: Optional["PeakGeometry"] = None

    @classmethod
    def from_moments(cls, mom, peaks, eps=1e-12):
        g, rho, gxz = mom.gamma_peak, mom.gamma_cr, mom.gamma_xz
        det = max(1.0 - rho**2, eps)
        a = (g - rho * gxz) / det
        b = (gxz - rho * g) / det
        s2 = max(1.0 - (g * a + gxz * b), eps)
        log_pref = peaks.log_volume + 0.5 * LOG_2PI - 2.0 * LOG_2PI - 3.0 * np.log(mom.rstar_over_r)
        geo = cls(g, rho, gxz, log_pref, a, b, np.sqrt(s2))
        if peaks.density == "bbks":
            geo.table = _table_for(geo.s)
        return geo

    def gaussian_only(self):
        """Geometry for counting peaks without the zeta_r variable (Gaussian zeta)."""
        if self._gaussian is None:
            s = np.sqrt(max(1.0 - self.gamma**2, 1e-12))
            geo = PeakGeometry(self.gamma, 0.0, 0.0, self.log_prefactor, self.gamma, 0.0, s)
            if self.table is not None:
                geo.table = _table_for(s)
            self._gaussian = geo
        return self._gaussian

    def log_weight(self, nu, z, sign, peaks):
        if peaks.density == "high_peak":
            gn = self.gamma * nu * sign
            with np.errstate(divide="ignore", invalid="ignore"):
                return self.log_prefactor + 3.0 * np.log(np.where(gn > 0, gn, 0.0))
        m = sign * (self.a * nu + self.b * z)
        return self.log_prefactor + self.table(m)


def log_beta(mom, collapse, ng: Optional[LocalNG] = None, statistic="threshold",
             peaks=PeaksConfig(), quad=QuadratureConfig(), geometry=None):
    """log beta at the smoothing scale of `mom` (CompactionMoments at the actual amplitude).

    `geometry` (a PeakGeometry built from the same spectrum shape) can be passed to reuse
    the tabulated G(m, s); it depends only on correlation coefficients, not on the amplitude.
    """
    if statistic not in ("threshold", "peaks"):
        raise ValueError("statistic must be 'threshold' or 'peaks'")
    sc, sr = mom.sigma_c, mom.sigma_r
    c_l, log_wc = _cl_nodes(collapse, quad)
    log_mass = collapse.log_mass_ratio(c_l)
    if statistic == "peaks" and geometry is None:
        geometry = PeakGeometry.from_moments(mom, peaks)

    if ng is None or ng.is_gaussian:
        nu = c_l / sc
        terms = -0.5 * nu**2 - 0.5 * LOG_2PI - np.log(sc) + log_mass + log_wc
        if statistic == "peaks":
            terms = terms + geometry.gaussian_only().log_weight(nu, 0.0, 1.0, peaks)
        else:
            terms = terms + np.log(quad.ps_factor)
        return float(logsumexp(terms))

    rho = mom.gamma_cr
    one_m = 1.0 - rho**2
    if one_m < 1e-10:
        raise ValueError("C_g and zeta_r are perfectly correlated; use a finite-width spectrum")
    # the domain of F always contains zeta = 0; scaling a bound by (1 - 1e-9) moves it inside
    lo_b, hi_b = ng.bounds()
    z_lo = max(-quad.z_span, lo_b / sr * (1 - 1e-9))
    z_hi = min(quad.z_span, hi_b / sr * (1 - 1e-9))

    def terms_on(z, log_wz):
        dF = ng.dF(sr * z)
        keep = dF != 0.0
        z, dF, log_wz = z[keep], dF[keep], log_wz[keep]
        nu = c_l[:, None] / (sc * dF[None, :])
        zz = z[None, :]
        log_p2 = -(nu**2 - 2.0 * rho * nu * zz + zz**2) / (2.0 * one_m) - LOG_2PI - 0.5 * np.log(one_m)
        t = log_p2 - np.log(sc * np.abs(dF))[None, :] + (log_mass + log_wc)[:, None] + log_wz[None, :]
        if statistic == "peaks":
            return t + geometry.log_weight(nu, zz, np.sign(dF)[None, :], peaks), z
        return t + np.log(quad.ps_factor), z

    def trapezoid(a, b, step):
        n = max(int(np.ceil((b - a) / step)), 2) + 1
        z = np.linspace(a, b, n)
        log_wz = np.full(n, np.log(z[1] - z[0]))
        log_wz[[0, -1]] -= np.log(2.0)
        return z, log_wz

    width = np.sqrt(one_m)
    # coarse pass: locate the z window that carries the integral (the ridge has width >= width)
    zc, wc = trapezoid(z_lo, z_hi, min(0.25, 2.0 * width))
    tc, zc_kept = terms_on(zc, wc)
    relevant = np.any(tc > tc.max() - 60.0, axis=0)
    step_c = zc[1] - zc[0]
    a = max(z_lo, zc_kept[relevant].min() - 2.0 * step_c)
    b = min(z_hi, zc_kept[relevant].max() + 2.0 * step_c)
    z, log_wz = trapezoid(a, b, min(quad.dz, width / 12.0))
    terms, _ = terms_on(z, log_wz)
    return float(logsumexp(terms))
