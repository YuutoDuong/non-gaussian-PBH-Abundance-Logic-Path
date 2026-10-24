"""Scalar-induced gravitational waves from a peaked curvature spectrum (Gaussian zeta, radiation era).

    Omega_GW,RD(k) = P_h(k)/24 = (1/6) int_1^inf dt int_0^1 ds K(t,s) P(k(t+s)/2) P(k(t-s)/2)
    K(t,s) = [(t^2-1)(1-s^2)/(t^2-s^2)]^2 I^2_{t,s}                     (2306.17149 eqs. 6-7)
    h^2 Omega_GW(f) = h^2 Omega_r (g*/g*0) (g*s/g*s0)^(-4/3) Omega_GW,RD(k),  f = k/(2 pi)  (eq. 3)

The same kernel as Kohri & Terada (1804.08577) in t = u + v, s = |u - v|. Checks in the tests:
the Dirac-delta closed form, Omega_GW,RD = 0.8222 A^2 for a flat spectrum, and the k^3 infrared
tail. Primordial NG corrections to the SIGW (relative size ~ A (3 f_NL/5)^2) are not included,
as in 2306.17149.

Quadrature: t is split at the resonance t = sqrt(3), where the kernel has a log singularity and
a step, and mapped to w = ln|t - sqrt(3)| on each side; Gauss-Legendre panels in w and in s.
The kernel is computed once per grid, so a spectrum costs two evaluations of P per node and k.
"""

from dataclasses import dataclass

import numpy as np

from . import thermal

SQRT3 = np.sqrt(3.0)
_N = np.arange(3, 15)
_G_SERIES = (2.0 - _N) / (2.0 * _N * (_N - 1.0))  # g(e)/e^3 = sum c_n e^(n-3), c_n = (2-n)/(2n(n-1))


def kernel(t, s, t2m3=None):
    """K(t, s) of 2306.17149 eqs. (6)-(7).

    The logarithmic bracket of eq. (7) equals y g(e) with y = 3 - s^2, e = (t^2 - s^2)/y and
    g(e) = e + ln|1-e| - e ln|1-e|/2, which vanishes like -e^3/12 at the corner t, s -> 1 while
    the prefactor diverges like e^-6; the ratio g/e^3 is taken from its series for small e.
    t2m3 = t^2 - 3 can be passed exactly when t is within rounding of the resonance sqrt(3).
    """
    t, s = np.broadcast_arrays(np.asarray(t, dtype=float), np.asarray(s, dtype=float))
    t2, s2 = t * t, s * s
    t2m3 = t2 - 3.0 if t2m3 is None else np.broadcast_to(t2m3, t.shape)
    y = 3.0 - s2
    xp = t2m3 + y  # t^2 - s^2
    xm = t2m3 - y  # t^2 + s^2 - 6
    e = xp / y
    small = e < 0.05
    with np.errstate(divide="ignore", invalid="ignore"):
        log_term = np.log(np.abs(t2m3) / y)  # ln|1 - e|
        g_e3 = np.where(small, 0.0, (e + log_term - 0.5 * e * log_term) / e**3)
    if np.any(small):
        g_e3 = np.where(small, np.polynomial.polynomial.polyval(np.where(small, e, 0.0), _G_SERIES), g_e3)
    with np.errstate(divide="ignore", invalid="ignore"):
        res = np.where(t2m3 > 0.0, 0.25 * np.pi**2 * xm**2 / xp**6, 0.0)
    i2 = 288.0 * xm**2 * (res + g_e3**2 / y**4)
    proj = ((t2 - 1.0) * (1.0 - s2) / xp) ** 2
    return proj * i2


def omega_rd_delta(x):
    """Omega_GW,RD / A^2 for P = A delta(ln k/k*), x = k/k* (Kohri & Terada 1804.08577)."""
    x = np.asarray(x, dtype=float)
    out = np.zeros_like(x)
    ok = (x > 0) & (x < 2.0)
    out[ok] = kernel(2.0 / x[ok], 0.0) / (6.0 * x[ok] ** 2)
    return out


def _gl_panels(edges, n):
    x, w = np.polynomial.legendre.leggauss(n)
    a, b = edges[:-1, None], edges[1:, None]
    return (0.5 * (b - a) * x + 0.5 * (b + a)).ravel(), (0.5 * (b - a) * w).ravel()


def _w_edges(w_hi, h):
    coarse = np.array([-40.0, -30.0, -22.0, -16.0, -12.0, -9.0, -7.0, -5.5])
    n = max(int(np.ceil((w_hi - coarse[-1]) / h)), 1)
    return np.concatenate([coarse, np.linspace(coarse[-1], w_hi, n + 1)[1:]])


@dataclass
class SIGWGrid:
    """Quadrature nodes in (t, s) with the kernel and weights folded together."""

    t_max: float = 1e8
    h: float = 0.25  # Gauss-Legendre panel width in w = ln|t - sqrt(3)|
    n_gl: int = 12
    n_s: int = 96

    def __post_init__(self):
        w1, a1 = _gl_panels(_w_edges(np.log(SQRT3 - 1.0), self.h), self.n_gl)
        w2, a2 = _gl_panels(_w_edges(np.log(self.t_max - SQRT3), self.h), self.n_gl)
        d = np.concatenate([-np.exp(w1), np.exp(w2)])  # t - sqrt(3), exact
        t = SQRT3 + d
        wt = np.concatenate([a1, a2]) * np.abs(d)
        s, ws = _gl_panels(np.array([0.0, 0.5, 1.0]), self.n_s // 2)
        self.t, self.s = t[:, None], s[None, :]
        k = kernel(self.t, self.s, t2m3=(d * (2.0 * SQRT3 + d))[:, None])
        self.weight = wt[:, None] * ws[None, :] * k / 6.0
        self.u = 0.5 * (self.t + self.s)  # momenta of the two scalar modes, in units of k
        self.v = 0.5 * (self.t - self.s)


def grid_for(spectrum, k_min, rel=1e-10, h_factor=0.5, t_max_cap=np.inf, **kw):
    """A grid that resolves the spectrum's narrowest feature and reaches its UV support at k_min.

    For slowly decaying UV tails (BPL with small beta) the support is astronomically far; the
    integrand falls like t^-(4 + 2 beta) there, so t_max can be capped (1e8 is ample for x >= 1e-4).
    """
    k_hi = spectrum.support(rel)[1]
    h = min(0.5, h_factor * spectrum.lnk_scale)
    return SIGWGrid(t_max=min(max(2.0 * k_hi / k_min, 10.0), t_max_cap), h=h, **kw)


def omega_rd(spectrum, k, grid=None, chunk=8):
    """Omega_GW,RD(k) for the curvature spectrum (any PowerSpectrum with a pointwise shape)."""
    k = np.atleast_1d(np.asarray(k, dtype=float))
    if grid is None:
        grid = grid_for(spectrum, k.min())
    out = np.empty_like(k)
    for i in range(0, k.size, chunk):
        kk = k[i : i + chunk, None, None]
        p = spectrum(kk * grid.u) * spectrum(kk * grid.v)
        out[i : i + chunk] = np.sum(grid.weight * p, axis=(1, 2))
    return out


def omega_rd_ir(spectrum, k, n_q=4000):
    """Leading infrared tail k << k*: K -> 72 (1-s^2)^2 t^-4 [pi^2 + ln^2(e^2 (3-s^2)/t^2)].

    Drops only the s-dependence of P (2306.17149 eq. S1 also drops it in the logarithm).
    """
    k = np.atleast_1d(np.asarray(k, dtype=float))
    lo, hi = spectrum.support(1e-12)
    lnq = np.linspace(np.log(lo), np.log(hi), n_q)
    q = np.exp(lnq)
    s, ws = np.polynomial.legendre.leggauss(64)
    s, ws = 0.5 * (s + 1.0), 0.5 * ws
    out = np.empty_like(k)
    for i, kk in enumerate(k):
        t = 2.0 * q / kk
        logs = np.log(np.e**2 * (3.0 - s[None, :] ** 2) / t[:, None] ** 2)
        ang = np.sum(ws * (1.0 - s**2) ** 2 * (np.pi**2 + logs**2), axis=1)
        # dt = t dln q;  Omega = (1/6) int dt 72 t^-4 P(q)^2 ang
        out[i] = 12.0 * np.trapezoid(t**-3 * spectrum(q) ** 2 * ang, lnq)
    return out


def h2_omega_today(f, omega_rd_values):
    """Present-day h^2 Omega_GW at frequencies f [Hz] from Omega_GW,RD at k = 2 pi f."""
    return thermal.H2_OMEGA_R * thermal.gw_dilution(f) * omega_rd_values


@dataclass
class SIGWTemplate:
    """Omega_GW,RD(k)/A^2 tabulated in x = k/k* for one spectral shape; scales as A^2."""

    spectrum: object  # PowerSpectrum; amplitude and k* are taken as 1 internally
    x_min: float = 1e-4  # with x_max: the 14 NG15 bins for 10^5 <= k* <= 10^10 Mpc^-1
    x_max: float = 30.0
    n_x: int = 600
    grid_kw: dict = None  # passed to grid_for

    def __post_init__(self):
        unit = self.spectrum.with_amplitude(1.0)
        if hasattr(unit, "k_star"):
            from dataclasses import replace

            unit = replace(unit, k_star=1.0)
        self._unit = unit
        self.ln_x = np.linspace(np.log(self.x_min), np.log(self.x_max), self.n_x)
        om = omega_rd(unit, np.exp(self.ln_x), grid=grid_for(unit, self.x_min, **(self.grid_kw or {})))
        self.ln_omega = np.log(np.maximum(om, 1e-300))

    def __call__(self, x):
        """Omega_GW,RD / A^2 at x = k/k* (log-log interpolation)."""
        return np.exp(np.interp(np.log(x), self.ln_x, self.ln_omega, left=-np.inf, right=-np.inf))

    def h2_omega(self, f, A, k_star):
        """h^2 Omega_GW today at frequencies f [Hz] for amplitude A and peak k* [Mpc^-1]."""
        x = np.asarray(f, dtype=float) / thermal.HZ_PER_INV_MPC / k_star
        return thermal.H2_OMEGA_R * thermal.gw_dilution(f) * A**2 * self(x)
