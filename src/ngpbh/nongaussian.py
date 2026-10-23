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
