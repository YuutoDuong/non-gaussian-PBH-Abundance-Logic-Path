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
exact for F' = const.  PeaksConfig(constraint="C_g") selects it (Phase 3a, the default).

PeaksConfig(constraint="C_l") counts maxima of C_l itself (Phase 3b, docs/phase3b_derivation.tex).
Since grad C_l = F'(z) grad(C_g + lam zeta_r) with lam = C_l F''(z) / F'(z)^2, these are BBKS
maxima of the Gaussian field Y = C_g + lam zeta_r, with lam frozen at each (C_l, zeta_r) node:
the same formula with R_*, m and s of Y.  The rank-one Hessian term kappa grad zeta_r grad zeta_r^T
is dropped; it vanishes for the logarithmic template and is suppressed by one power of sigma
otherwise (see the derivation and tests/test_ngpeaks.py).

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


def bbks_shape_moments(x, n=96, y_max=3.0):
    """f(x) and f2(x) = E[e2(-zeta) Theta_peak | x] by quadrature over the BBKS shape variables.

    With eigenvalues l1 >= l2 >= l3 of -zeta, x = sum l, y = (l1 - l3)/2, z = (l1 - 2 l2 + l3)/2,
    the traceless part has density ~ y (y^2 - z^2) exp[-5/2 (3y^2 + z^2)] on y >= |z|, the peak
    condition is l3 = (x - 3y + z)/3 > 0, and e2 = l1 l2 + l1 l3 + l2 l3 = x^2/3 - (3y^2 + z^2)/3.
    f = 27 E[l1 l2 l3 Theta | x] reproduces bbks_f to 1e-10, which fixes the normalisation;
    f2 -> x^2/3 - 1/3 for large x.  f2 sets the rank-one correction of Phase 3b peaks."""
    t, w = _leggauss(n)
    z_max = y_max  # y >= |z|
    x = np.atleast_1d(np.asarray(x, dtype=float))[:, None, None]
    z_lo, z_hi = np.maximum(-x / 4.0, -z_max), np.minimum(x / 2.0, z_max)
    z = 0.5 * (z_hi - z_lo) * t[None, :, None] + 0.5 * (z_hi + z_lo)
    wz = 0.5 * (z_hi - z_lo) * w[None, :, None]
    y_lo = np.abs(z)
    y_hi = np.maximum(np.minimum((x + z) / 3.0, y_max), y_lo)
    y = 0.5 * (y_hi - y_lo) * t[None, None, :] + 0.5 * (y_hi + y_lo)
    wy = 0.5 * (y_hi - y_lo) * w[None, None, :]
    weight = y * (y**2 - z**2) * np.exp(-2.5 * (3.0 * y**2 + z**2)) * wz * wy
    det = (x + 3.0 * y + z) * (x - 2.0 * z) * (x - 3.0 * y + z) / 27.0
    e2 = x**2 / 3.0 - (3.0 * y**2 + z**2) / 3.0
    # normalisation: the same density over its full domain
    zz = z_max * t[:, None]
    yy = 0.5 * (y_max - np.abs(zz)) * t[None, :] + 0.5 * (y_max + np.abs(zz))
    norm = np.sum(yy * (yy**2 - zz**2) * np.exp(-2.5 * (3.0 * yy**2 + zz**2))
                  * z_max * w[:, None] * 0.5 * (y_max - np.abs(zz)) * w[None, :])
    return 27.0 * np.sum(det * weight, axis=(1, 2)) / norm, np.sum(e2 * weight, axis=(1, 2)) / norm


class _F2:
    """f2(x) on a fine grid (linear interpolation, relative error < 1e-5), x^2/3 - 1/3 beyond."""

    def __init__(self, x_hi=14.0, dx=0.005):
        self.x = np.arange(0.0, x_hi + dx / 2, dx)
        self.vals = bbks_shape_moments(self.x, n=48)[1]  # n = 48 agrees with n = 96 to 1e-10
        self.vals[0] = 0.0

    def __call__(self, x):
        x = np.asarray(x, dtype=float)
        out = np.interp(x, self.x, self.vals, left=0.0)
        high = x > self.x[-1]
        out[high] = x[high] ** 2 / 3.0 - 1.0 / 3.0
        return out


_F2_CACHE = []


def _f2():
    if not _F2_CACHE:
        _F2_CACHE.append(_F2())
    return _F2_CACHE[0]


def _log_gauss_average(m, s, log_f, n_x):
    """log int_0^inf exp(log_f(x)) N(x; m, s^2) dx, vectorised over m (s scalar or like m)."""
    m = np.atleast_1d(np.asarray(m, dtype=float))
    s = np.broadcast_to(np.asarray(s, dtype=float), m.shape)
    out = np.full(m.shape, -np.inf)
    lo = np.maximum(m - 12.0 * s, 0.0)
    hi = np.maximum(m + 12.0 * s, lo + 1e-12)
    ok = m + 12.0 * s > 0
    if not np.any(ok):
        return out
    nodes, weights = _leggauss(n_x)
    lo_, hi_, s_ = lo[ok, None], hi[ok, None], s[ok, None]
    x = 0.5 * (hi_ - lo_) * nodes + 0.5 * (hi_ + lo_)
    log_gauss = -0.5 * ((x - m[ok, None]) / s_) ** 2 - np.log(np.sqrt(2.0 * np.pi) * s_)
    with np.errstate(divide="ignore"):
        terms = log_f(x) + log_gauss + np.log(0.5 * (hi_ - lo_) * weights)
    out[ok] = logsumexp(terms, axis=1)
    return out


def log_G(m, s, n_x=200):
    """log of G(m, s) = int_0^inf f(x) N(x; m, s^2) dx, vectorised over m (s scalar or like m)."""
    with np.errstate(divide="ignore"):
        return _log_gauss_average(m, s, lambda x: np.log(bbks_f(x)), n_x)


def log_G2(m, s, n_x=200):
    """log of G2(m, s) = int_0^inf f2(x) N(x; m, s^2) dx (rank-one correction, Phase 3b)."""
    f2 = _f2()
    with np.errstate(divide="ignore"):
        return _log_gauss_average(m, s, lambda x: np.log(f2(x)), n_x)


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


class _LogGTable2D:
    """log G(m, s) and R2(m, s) = 27 G2 / G on a uniform (ln s, m) grid, for Phase 3b where s
    changes from node to node.

    Bilinear interpolation: it cannot overshoot next to the region where G underflows (stored
    as -700), and its error (< 5e-4 in log G for m > 1, < 1.1e-4 for m > 3, where peaks above
    threshold sit) moves A_* by < 1e-5.  48 quadrature nodes are accurate to 1e-6.  Above m_hi
    the asymptotic forms are used (G = m^3 - 3m + 3 m s^2, G2 = (m^2 + s^2 - 1)/3).  Below
    m_lo, G < exp(-70) for every s <= 1 and is set to zero, as the 1D table does below -12 s;
    below s_lo, direct quadrature."""

    def __init__(self, s_lo=5e-3, m_lo=-12.0, m_hi=40.0, n_s=161, dm=0.02, n_x=48):
        self.ln_s0, self.ln_s1, self.n_s = np.log(s_lo), 0.0, n_s
        self.m_lo, self.m_hi = m_lo, m_hi
        self.n_m = int(round((m_hi - m_lo) / dm)) + 1
        self.dm = (m_hi - m_lo) / (self.n_m - 1)
        self.d_ln_s = (self.ln_s1 - self.ln_s0) / (n_s - 1)
        m = np.linspace(m_lo, m_hi, self.n_m)
        s = np.exp(np.linspace(self.ln_s0, self.ln_s1, n_s))
        log_g = np.array([log_G(m, si, n_x) for si in s])
        log_g2 = np.array([log_G2(m, si, n_x) for si in s])
        self.vals = np.maximum(log_g, -700.0)
        with np.errstate(invalid="ignore"):
            r2 = 27.0 * np.exp(log_g2 - log_g)
        # where G underflows R2 is irrelevant; keep it finite and bounded by its last valid value
        self.r2 = np.where(np.isfinite(r2), r2, 0.0)

    @staticmethod
    def _asymptotic(m, s):
        g = m**3 - 3.0 * m + 3.0 * m * s**2
        return np.log(g), 9.0 * (m**2 + s**2 - 1.0) / g

    def lookup(self, m, s):
        """(log G, R2) at arrays m, s."""
        m, s = np.broadcast_arrays(np.asarray(m, dtype=float), np.asarray(s, dtype=float))
        shape = m.shape
        m, s = m.ravel(), np.minimum(s.ravel(), 1.0)
        log_g, r2 = np.full(m.shape, -np.inf), np.zeros(m.shape)
        ln_s = np.log(np.maximum(s, 1e-300))
        inside = (ln_s >= self.ln_s0) & (m >= self.m_lo) & (m <= self.m_hi)
        u = (ln_s[inside] - self.ln_s0) / self.d_ln_s
        v = (m[inside] - self.m_lo) / self.dm
        i = np.minimum(u.astype(np.intp), self.n_s - 2)
        j = np.minimum(v.astype(np.intp), self.n_m - 2)
        fu, fv = u - i, v - j
        k = i * self.n_m + j
        for out, t in ((log_g, self.vals.ravel()), (r2, self.r2.ravel())):
            t00, t01 = np.take(t, k), np.take(t, k + 1)
            t10, t11 = np.take(t, k + self.n_m), np.take(t, k + self.n_m + 1)
            out[inside] = (1 - fu) * (t00 + fv * (t01 - t00)) + fu * (t10 + fv * (t11 - t10))
        high = m > self.m_hi
        log_g[high], r2[high] = self._asymptotic(m[high], s[high])
        rest = (ln_s < self.ln_s0) & (m >= self.m_lo) & (m <= self.m_hi)
        if np.any(rest):
            log_g[rest] = log_G(m[rest], s[rest])
            with np.errstate(invalid="ignore"):
                r2_rest = 27.0 * np.exp(log_G2(m[rest], s[rest]) - log_g[rest])
            r2[rest] = np.where(np.isfinite(r2_rest), r2_rest, 0.0)
        return log_g.reshape(shape), r2.reshape(shape)

    def __call__(self, m, s):
        return self.lookup(m, s)[0]


_TABLE_2D = []


def _table_2d():
    if not _TABLE_2D:
        _TABLE_2D.append(_LogGTable2D())
    return _TABLE_2D[0]


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
    constraint: str = "C_g"  # field whose maxima are counted: "C_g" (Phase 3a) or "C_l" (Phase 3b)
    rank_one: bool = True  # Phase 3b: include the kappa grad zeta_r grad zeta_r^T Hessian term

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


@dataclass
class NGPeakGeometry:
    """Phase 3b: maxima of C_l = F'(zeta_r) C_g are BBKS maxima of Y = C_g + lam zeta_r,
    corrected for the rank-one Hessian term.

    Stores the unit-amplitude CompactionMoments; log_weight_ng evaluates, at each (C_l, zeta_r)
    node, R_*, the conditional mean m = a nu + b z and variance s^2 of x = -lap Y / sigma_2(Y)
    (docs/phase3b_derivation.tex, eqs. 12-17), and the factor 1 - tau R2 with
    tau = sign(F') kappa v / sigma_2(Y) and R2 = 27 G2 / G (Sec. 7).  Y is normalised by
    (1 + lam^2)^(1/2), i.e. lam = tan(theta), so that F' -> 0 stays finite.  lam = 0
    reproduces PeakGeometry."""

    rho: float
    sc2: float
    sr2: float
    s1: float
    scr1: float
    s1r: float
    s2: float
    s2cr: float
    s2r: float
    log_volume: float
    base: PeakGeometry

    @classmethod
    def from_moments(cls, mom, peaks):
        if peaks.density != "bbks":
            raise ValueError("Phase 3b peaks need density='bbks'")
        if not np.isfinite(mom.s2r_2):
            raise ValueError("moments lack the zeta_r gradient moments needed for Phase 3b")
        return cls(mom.gamma_cr, mom.sigma_c2, mom.sigma_r2, mom.s1_2, mom.scr1_2, mom.s1r_2,
                   mom.s2_2, mom.s2cr_2, mom.s2r_2, peaks.log_volume,
                   PeakGeometry.from_moments(mom, peaks))

    def gaussian_only(self):
        return self.base.gaussian_only()

    def log_weight_ng(self, nu, z, c_l, dF, d2F, d3F, sqrt_amp, rank_one=True, eps=1e-12):
        """Peak weight at nodes (nu, z) = (C_g / sigma_c, zeta_r / sigma_r) with C_l = F' C_g
        and F', F'', F''' at zeta_r (arrays broadcastable against nu); sqrt_amp = A^(1/2)."""
        # lam = C_l F'' / F'^2 = tan(theta), theta in [-pi/2, pi/2]
        num, den = c_l * d2F, dF**2
        hyp = np.hypot(num, den)
        ct, st = den / hyp, num / hyp
        s1l = ct**2 * self.s1 + 2.0 * ct * st * self.scr1 + st**2 * self.s1r
        s2l = ct**2 * self.s2 + 2.0 * ct * st * self.s2cr + st**2 * self.s2r
        g = (ct * self.s1 + st * self.scr1) / np.sqrt(s2l * self.sc2)
        gz = (ct * self.scr1 + st * self.s1r) / np.sqrt(s2l * self.sr2)
        det = 1.0 - self.rho**2
        a = (g - self.rho * gz) / det
        b = (gz - self.rho * g) / det
        s = np.sqrt(np.maximum(1.0 - a * g - b * gz, eps))
        sign = np.sign(dF)
        log_pref = self.log_volume - 1.5 * LOG_2PI - 1.5 * np.log(3.0 * s1l / s2l)
        log_g, r2 = _table_2d().lookup(sign * (a * nu + b * z), s)
        if not rank_one:
            return log_pref + log_g
        # tau = sign kappa v / sigma_2(Y), kappa cos^3(theta) = C_l F''' ct^3 / F'^2 - 2 st^2 ct F' / C_l
        with np.errstate(over="ignore", invalid="ignore", divide="ignore"):
            kappa_c3 = np.nan_to_num(c_l * d3F * ct**3 / dF**2) - 2.0 * st**2 * ct * dF / c_l
        grad = self.s1 * self.s1r - self.scr1**2
        tau = sign * sqrt_amp * grad / (3.0 * s1l * np.sqrt(s2l)) * kappa_c3
        u = tau * r2
        # exact to first order in tau; beyond u = 1/2 a smooth positive continuation (see tests)
        log_q = np.where(u <= 0.5, np.log1p(-np.minimum(u, 0.5)), np.log(0.5) - 2.0 * (u - 0.5))
        return log_pref + log_g + log_q


def make_peak_geometry(mom, peaks):
    """The geometry log_beta needs for `peaks`, from moments at any amplitude."""
    if peaks.constraint == "C_l":
        return NGPeakGeometry.from_moments(mom, peaks)
    if peaks.constraint == "C_g":
        return PeakGeometry.from_moments(mom, peaks)
    raise ValueError("constraint must be 'C_g' or 'C_l'")


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
        geometry = make_peak_geometry(mom, peaks)

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
    sqrt_amp = np.sqrt(mom.sigma_c2 / geometry.sc2) if isinstance(geometry, NGPeakGeometry) else 1.0
    # the domain of F always contains zeta = 0; scaling a bound by (1 - 1e-9) moves it inside
    lo_b, hi_b = ng.bounds()
    z_lo = max(-quad.z_span, lo_b / sr * (1 - 1e-9))
    z_hi = min(quad.z_span, hi_b / sr * (1 - 1e-9))

    def terms_on(z, log_wz, floor=-np.inf):
        """Integrand on the (C_l, z) grid; Phase 3b peak weights only where t > floor."""
        dF = ng.dF(sr * z)
        keep = dF != 0.0
        z, dF, log_wz = z[keep], dF[keep], log_wz[keep]
        nu = c_l[:, None] / (sc * dF[None, :])
        zz = z[None, :]
        log_p2 = -(nu**2 - 2.0 * rho * nu * zz + zz**2) / (2.0 * one_m) - LOG_2PI - 0.5 * np.log(one_m)
        t = log_p2 - np.log(sc * np.abs(dF))[None, :] + (log_mass + log_wc)[:, None] + log_wz[None, :]
        sign = np.sign(dF)[None, :]
        if statistic == "peaks" and isinstance(geometry, NGPeakGeometry):
            row, col = np.nonzero(t > floor)
            w = np.full(t.shape, -np.inf)
            w[row, col] = geometry.log_weight_ng(nu[row, col], z[col], c_l[row], dF[col],
                                                 ng.d2F(sr * z)[col], ng.d3F(sr * z)[col],
                                                 sqrt_amp, peaks.rank_one)
            return t + w, z, w
        if statistic == "peaks":
            return t + geometry.log_weight(nu, zz, sign, peaks), z, None
        return t + np.log(quad.ps_factor), z, None

    def trapezoid(a, b, step):
        n = max(int(np.ceil((b - a) / step)), 2) + 1
        z = np.linspace(a, b, n)
        log_wz = np.full(n, np.log(z[1] - z[0]))
        log_wz[[0, -1]] -= np.log(2.0)
        return z, log_wz

    width = np.sqrt(one_m)
    # coarse pass: locate the z window that carries the integral (the ridge has width >= width)
    zc, wc = trapezoid(z_lo, z_hi, min(0.25, 2.0 * width))
    tc, zc_kept, w_c = terms_on(zc, wc)
    if not np.isfinite(tc.max()):
        # no node in the window can host a peak (Phase 3b at absurdly small amplitudes, where
        # the ridge z = rho nu lies far outside z_span): beta is zero to working precision
        return -np.inf
    relevant = np.any(tc > tc.max() - 60.0, axis=0)
    step_c = zc[1] - zc[0]
    a = max(z_lo, zc_kept[relevant].min() - 2.0 * step_c)
    b = min(z_hi, zc_kept[relevant].max() + 2.0 * step_c)
    z, log_wz = trapezoid(a, b, min(quad.dz, width / 12.0))
    floor = -np.inf
    if w_c is not None:
        # a fine node matters only if t + w > max - 50; w varies smoothly, so bound it by the
        # largest coarse weight plus a margin (the fine step is >= 12x smaller than the coarse)
        floor = tc.max() - 50.0 - (w_c[np.isfinite(w_c)].max() + 10.0) - np.log(step_c / (z[1] - z[0]))
    terms, _, _ = terms_on(z, log_wz, floor)
    return float(logsumexp(terms))
