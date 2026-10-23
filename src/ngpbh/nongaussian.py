"""Local non-Gaussian templates zeta = F(zeta_g), with zeta_g Gaussian.

Only F' enters the compaction function, C_l = F'(zeta_g) C_g (2402.11033 eq. 40), so a
constant shift in F (e.g. subtracting <zeta_g^2>) is irrelevant and is dropped.
"""

from dataclasses import dataclass

import numpy as np


class LocalNG:
    """Interface: F and its derivatives, and the open interval of zeta_g where F is defined."""

    def F(self, z):
        raise NotImplementedError

    def dF(self, z):
        raise NotImplementedError

    def d2F(self, z):
        raise NotImplementedError

    def d3F(self, z):
        raise NotImplementedError

    def bounds(self):
        return -np.inf, np.inf

    @property
    def is_gaussian(self):
        return False


@dataclass(frozen=True)
class Gaussian(LocalNG):
    def F(self, z):
        return np.asarray(z, dtype=float)

    def dF(self, z):
        return np.ones_like(np.asarray(z, dtype=float))

    def d2F(self, z):
        return np.zeros_like(np.asarray(z, dtype=float))

    def d3F(self, z):
        return np.zeros_like(np.asarray(z, dtype=float))

    @property
    def is_gaussian(self):
        return True


@dataclass(frozen=True)
class Quadratic(LocalNG):
    """zeta = zeta_g + (3/5) f_NL zeta_g^2."""

    f_nl: float

    def F(self, z):
        z = np.asarray(z, dtype=float)
        return z + 0.6 * self.f_nl * z**2

    def dF(self, z):
        return 1.0 + 1.2 * self.f_nl * np.asarray(z, dtype=float)

    def d2F(self, z):
        return np.full_like(np.asarray(z, dtype=float), 1.2 * self.f_nl)

    def d3F(self, z):
        return np.zeros_like(np.asarray(z, dtype=float))

    @property
    def is_gaussian(self):
        return self.f_nl == 0.0


@dataclass(frozen=True)
class Logarithmic(LocalNG):
    """Ultra-slow-roll template zeta = -mu ln(1 - zeta_g / mu) (2402.11033 eq. 60).

    mu > 0 gives positive non-Gaussianity (exponential tail), mu < 0 the negative branch
    of the logarithmic family; to quadratic order f_NL = 5 / (6 mu).  The map is defined
    for 1 - zeta_g / mu > 0; beyond it lies the trapped-vacuum regime, which is excluded.
    """

    mu: float

    @classmethod
    def from_fnl(cls, f_nl):
        return cls(5.0 / (6.0 * f_nl))

    @classmethod
    def from_beta(cls, beta):
        """2306.17149 eq. (12): zeta = -(2/beta) ln(1 - beta zeta_g / 2)."""
        return cls(2.0 / beta)

    @property
    def f_nl(self):
        return 5.0 / (6.0 * self.mu)

    def _u(self, z):
        return 1.0 - np.asarray(z, dtype=float) / self.mu

    def F(self, z):
        return -self.mu * np.log(self._u(z))

    def dF(self, z):
        return 1.0 / self._u(z)

    def d2F(self, z):
        return 1.0 / (self.mu * self._u(z) ** 2)

    def d3F(self, z):
        return 2.0 / (self.mu**2 * self._u(z) ** 3)

    def bounds(self):
        return (-np.inf, self.mu) if self.mu > 0 else (self.mu, np.inf)


@dataclass(frozen=True)
class Curvaton(LocalNG):
    """Sudden-decay curvaton with quadratic potential, zeta = ln X (2211.01728 eqs. 4-6).

    X = e^zeta is the positive root of X^4 - g X - (1 - Omega) = 0, with
    g = Omega (2 r + 3 zeta_g)^2 / (4 r^2) and Omega = 4 r / (3 + r) the curvaton density
    fraction at decay (r = r_dec).  This is the quartic behind the closed form of eq. (5); for
    r = 1 it gives zeta = (2/3) ln(1 + 3 zeta_g / 2) (eq. 8), and to quadratic order
    f_NL = (5/3)(3/(4r) - 1 - r/2) (eq. 9).
    """

    r_dec: float

    def __post_init__(self):
        if not 0.0 < self.r_dec <= 1.0:
            raise ValueError("r_dec must be in (0, 1]")

    @property
    def omega(self):
        return 4.0 * self.r_dec / (3.0 + self.r_dec)

    @property
    def f_nl(self):
        r = self.r_dec
        return 5.0 / 3.0 * (3.0 / (4.0 * r) - 1.0 - 0.5 * r)

    def _g(self, z):
        r = self.r_dec
        return self.omega * (2.0 * r + 3.0 * np.asarray(z, dtype=float)) ** 2 / (4.0 * r**2)

    def _X(self, z):
        g = self._g(z)
        c = 1.0 - self.omega
        x = 1.0 + np.cbrt(g) + c**0.25  # right of the root; f is convex on X > 0
        for _ in range(100):
            step = (x**4 - g * x - c) / (4.0 * x**3 - g)
            x = x - step
            if np.all(np.abs(step) <= 1e-15 * x):
                break
        return x

    def F(self, z):
        return np.log(self._X(z))

    def dF(self, z):
        r = self.r_dec
        g = self._g(z)
        dg = 3.0 * self.omega * (2.0 * r + 3.0 * np.asarray(z, dtype=float)) / (2.0 * r**2)
        return dg / (4.0 * self._X(z) ** 3 - g)

    def d2F(self, z, h=1e-6):
        return (self.dF(np.asarray(z) + h) - self.dF(np.asarray(z) - h)) / (2.0 * h)

    def d3F(self, z, h=1e-4):
        return (self.d2F(np.asarray(z) + h) - self.d2F(np.asarray(z) - h)) / (2.0 * h)

    def bounds(self):
        # for r = 1 the root X vanishes at zeta_g = -2/3 (zeta -> -infinity)
        return (-2.0 / 3.0, np.inf) if self.r_dec == 1.0 else (-np.inf, np.inf)
