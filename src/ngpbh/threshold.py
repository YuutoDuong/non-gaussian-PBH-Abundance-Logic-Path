"""Collapse thresholds on the peak value of the compaction function.

delta_c(q) is the Escriva-Germani-Sheth / Musco et al. relation between the threshold and
the curvature q = -r_m^2 C''(r_m) / (4 C(r_m)) of the compaction profile
(1907.13311, 2011.03014):

    delta_c(q) = (4/15) e^{-1/q} q^{1 - 5/(2q)} / gamma_lower(5/(2q), 1/q),

with delta_c -> 2/5 for q -> 0 (broadest profiles) and 2/3 for q -> infinity.

MuscoThreshold evaluates q on the mean profile of the spectrum, following the summary
of 2011.03014 given in the supplement of 2306.17149 (eqs. S4-S6).
"""

from dataclasses import dataclass

import numpy as np
from scipy.integrate import quad, simpson
from scipy.optimize import brentq
from scipy.special import gammainc, gammaln

from . import windows
from ._quad import leggauss as _leggauss
from .moments import MomentsConfig, _lnk_grids


def delta_c_of_q(q):
    q = np.asarray(q, dtype=float)
    s = 2.5 / q
    log_num = np.log(4.0 / 15.0) - 1.0 / q + (1.0 - s) * np.log(q)
    log_den = np.log(gammainc(s, 1.0 / q)) + gammaln(s)
    return np.exp(log_num - log_den)


def cl_threshold(C_th, phi=2.0 / 3.0):
    """Type-I root of C = C_l - C_l^2 / (4 Phi) = C_th."""
    return 2.0 * phi * (1.0 - np.sqrt(1.0 - C_th / phi))


@dataclass(frozen=True)
class ThresholdResult:
    C_th: float
    q: float = np.nan
    rm_kpeak: float = np.nan


@dataclass(frozen=True)
class ConstantThreshold:
    C_th: float

    def evaluate(self, spectrum):
        return ThresholdResult(self.C_th)


def _osc_integral(spectrum, r, transfer, a_sin, a_cos, x_split=2.0, x_max=2.0e3):
    """int dlnk P^T(k) [a_sin(k) sin(kr) + a_cos(k) cos(kr)] for smooth a_sin, a_cos.

    The mean profile carries no window, so without the transfer function the integrand only
    decays with the spectrum; the kr > x_split part then uses QUADPACK's Fourier quadrature.
    With the transfer function the integrand falls like (kr)^-3 and is truncated at x_max."""

    def integrand(k):
        x = k * r
        val = a_sin(k) * np.sin(x) + a_cos(k) * np.cos(x)
        if transfer:
            val = val * windows.transfer_rd(x) ** 2
        return val

    if spectrum.is_delta:
        k = np.array([spectrum.k_star])
        return float(spectrum.A * integrand(k)[0])
    if transfer:
        cfg = MomentsConfig(transfer=True, x_cut=None, x_max=x_max)
        grid_ln, grid_k = _lnk_grids(spectrum, r, cfg)
        total = 0.0
        if grid_ln is not None:
            k = np.exp(grid_ln)
            total += simpson(integrand(k) * spectrum(k), x=grid_ln)
        if grid_k is not None:
            total += simpson(integrand(grid_k) * spectrum(grid_k) / grid_k, x=grid_k)
        return total
    k_lo, k_hi = spectrum.support(1e-14)
    k_split = min(max(x_split / r, k_lo), k_hi)
    total = 0.0
    if k_split > k_lo:
        nodes, weights = _leggauss(256)
        a, b = np.log(k_lo), np.log(k_split)
        lnk = 0.5 * (b - a) * nodes + 0.5 * (b + a)
        k = np.exp(lnk)
        total += 0.5 * (b - a) * np.sum(weights * integrand(k) * spectrum(k))
    if k_hi > k_split:
        opts = dict(weight=None, limit=500, epsabs=0.0, epsrel=1e-10)
        for amp, w in ((a_sin, "sin"), (a_cos, "cos")):
            fn = lambda k, amp=amp: float(amp(np.array([k]))[0] * spectrum(np.array([k]))[0] / k)
            opts["weight"] = w
            total += quad(fn, k_split, k_hi, wvar=r, **opts)[0]
    return total


def rm_equation(spectrum, r, transfer=False):
    """Eq. S4: int dlnk [(k^2 r^2 - 1) sin(kr)/(kr) + cos(kr)] P^T = 0 at r = r_m.

    Equals -r d(r psi')/dr for the mean profile psi, so it changes sign from + to - at r_m."""
    return _osc_integral(
        spectrum,
        r,
        transfer,
        a_sin=lambda k: ((k * r) ** 2 - 1.0) / (k * r),
        a_cos=lambda k: np.ones_like(k),
    )


def musco_rm(spectrum, transfer=False, x_range=(0.5, 60.0), n_scan=80):
    """Innermost maximum r_m of the mean compaction profile."""
    kp = spectrum.k_peak
    xs = np.geomspace(*x_range, n_scan)
    vals = np.array([rm_equation(spectrum, x / kp, transfer) for x in xs])
    idx = np.where((vals[:-1] > 0) & (vals[1:] <= 0))[0]
    if len(idx) == 0:
        raise RuntimeError("no maximum of the mean compaction profile found in x_range")
    i = idx[0]
    return brentq(lambda x: rm_equation(spectrum, x / kp, transfer), xs[i], xs[i + 1], xtol=1e-12) / kp


def musco_q(spectrum, r_m, transfer=False):
    """Eq. S5: F(q)[1 + F(q)] q = -(1/2)[1 + r_m int dk k cos(k r_m) P^T / int dk sin(k r_m) P^T],
    with F(q) = sqrt(1 - (3/2) delta_c(q))."""
    zero = lambda k: np.zeros_like(k)
    num = _osc_integral(spectrum, r_m, transfer, a_sin=zero, a_cos=lambda k: k**2)
    den = _osc_integral(spectrum, r_m, transfer, a_sin=lambda k: k, a_cos=zero)
    rhs = -0.5 * (1.0 + r_m * num / den)

    def lhs(q):
        F = np.sqrt(max(1.0 - 1.5 * delta_c_of_q(q), 0.0))
        return F * (1.0 + F) * q - rhs

    return brentq(lhs, 1e-2, 1e4, xtol=1e-12)


@dataclass(frozen=True)
class MuscoThreshold:
    """Shape-dependent threshold delta_c(q) with q from the mean profile (2011.03014)."""

    transfer: bool = False

    def evaluate(self, spectrum):
        r_m = musco_rm(spectrum, self.transfer)
        q = musco_q(spectrum, r_m, self.transfer)
        return ThresholdResult(float(delta_c_of_q(q)), q, r_m * spectrum.k_peak)
