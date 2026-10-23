"""Primordial curvature power spectra P_zeta(k), with k in Mpc^-1.

Every spectrum is A * shape(k).  Variances are linear in A, so the abundance code
precomputes them at unit amplitude and rescales.
"""

from dataclasses import dataclass, replace

import numpy as np

SQRT2PI = np.sqrt(2.0 * np.pi)


class PowerSpectrum:
    A: float
    is_delta = False

    def shape(self, k):
        raise NotImplementedError

    def __call__(self, k):
        return self.A * self.shape(np.asarray(k, dtype=float))

    def with_amplitude(self, A):
        return replace(self, A=float(A))

    def support(self, rel=1e-16):
        """(k_lo, k_hi) outside which P / max(P) < rel."""
        raise NotImplementedError

    @property
    def lnk_scale(self):
        """Smallest scale in ln k on which the shape varies (sets the grid spacing)."""
        raise NotImplementedError

    @property
    def k_peak(self):
        raise NotImplementedError


@dataclass(frozen=True)
class LogNormal(PowerSpectrum):
    """Lognormal peak of width Delta.

    normalization="integrated": P = A / (sqrt(2 pi) Delta) exp(-ln^2(k/k*) / (2 Delta^2)),
    so that A = int P dln k (2306.17149, 2402.11033, the project plan).
    normalization="peak": P = A exp(...), A is the peak value (pbhstat).
    """

    A: float
    k_star: float
    width: float
    normalization: str = "integrated"

    def __post_init__(self):
        if self.normalization not in ("integrated", "peak"):
            raise ValueError("normalization must be 'integrated' or 'peak'")
        if self.width <= 0:
            raise ValueError("width must be positive; use DiracDelta for the monochromatic limit")

    @property
    def peak_value(self):
        if self.normalization == "peak":
            return self.A
        return self.A / (SQRT2PI * self.width)

    def shape(self, k):
        ln = np.log(k / self.k_star)
        norm = 1.0 if self.normalization == "peak" else 1.0 / (SQRT2PI * self.width)
        return norm * np.exp(-(ln**2) / (2.0 * self.width**2))

    def support(self, rel=1e-16):
        half = self.width * np.sqrt(2.0 * np.log(1.0 / rel))
        return self.k_star * np.exp(-half), self.k_star * np.exp(half)

    @property
    def lnk_scale(self):
        return self.width

    @property
    def k_peak(self):
        return self.k_star


@dataclass(frozen=True)
class BrokenPowerLaw(PowerSpectrum):
    """Broken power law of 2306.17149 eq. (1); A is the peak value at k = k*."""

    A: float
    k_star: float
    alpha: float = 4.0
    beta: float = 3.0
    gamma: float = 1.0

    def shape(self, k):
        x = k / self.k_star
        a, b, g = self.alpha, self.beta, self.gamma
        return (a + b) ** g / (b * x ** (-a / g) + a * x ** (b / g)) ** g

    def support(self, rel=1e-16):
        a, b, g = self.alpha, self.beta, self.gamma
        lo = self.k_star * (rel * (b / (a + b)) ** g) ** (1.0 / a)
        hi = self.k_star * (rel * (a / (a + b)) ** g) ** (-1.0 / b)
        return lo, hi

    @property
    def lnk_scale(self):
        return float(np.clip(self.gamma / max(self.alpha, self.beta), 0.05, 1.0))

    @property
    def k_peak(self):
        return self.k_star


@dataclass(frozen=True)
class FlatBand(PowerSpectrum):
    """Scale-invariant band, P = A for k_min <= k <= k_max (the broad spectrum of 2211.01728)."""

    A: float
    k_min: float
    k_max: float

    def shape(self, k):
        return np.where((k >= self.k_min) & (k <= self.k_max), 1.0, 0.0)

    def support(self, rel=1e-16):
        return self.k_min, self.k_max

    @property
    def lnk_scale(self):
        return np.log(self.k_max / self.k_min)

    @property
    def k_peak(self):
        return np.sqrt(self.k_min * self.k_max)


@dataclass(frozen=True)
class DiracDelta(PowerSpectrum):
    """Monochromatic spectrum, P = A delta(ln k - ln k*)."""

    A: float
    k_star: float
    is_delta = True

    def shape(self, k):
        raise TypeError("a Dirac delta has no pointwise shape; moments are evaluated analytically")

    def support(self, rel=1e-16):
        return self.k_star, self.k_star

    @property
    def lnk_scale(self):
        return 0.0

    @property
    def k_peak(self):
        return self.k_star
