"""Cached quadrature rules (numpy recomputes Gauss-Legendre nodes by eigendecomposition)."""

from functools import lru_cache

import numpy as np


@lru_cache(maxsize=None)
def leggauss(n):
    return np.polynomial.legendre.leggauss(n)
