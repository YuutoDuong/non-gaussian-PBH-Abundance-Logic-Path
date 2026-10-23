"""Fourier-space window and transfer functions, written in terms of x = k r."""

import numpy as np

SQRT3 = np.sqrt(3.0)


def _as_array(x):
    x = np.asarray(x, dtype=float)
    return np.atleast_1d(x), x.ndim == 0


def tophat(x):
    """Real-space top-hat window, W(x) = 3 (sin x - x cos x) / x^3."""
    xa, scalar = _as_array(x)
    out = np.empty_like(xa)
    small = np.abs(xa) < 1e-2
    xs = xa[small]
    out[small] = 1.0 - xs**2 / 10.0 + xs**4 / 280.0
    xl = xa[~small]
    out[~small] = 3.0 * (np.sin(xl) - xl * np.cos(xl)) / xl**3
    return out[0] if scalar else out


def shell(x):
    """Spherical-shell window, W_s(x) = sin x / x."""
    xa, scalar = _as_array(x)
    out = np.sinc(xa / np.pi)
    return out[0] if scalar else out


def gaussian(x):
    """Gaussian window in the convention W(x) = exp(-x^2 / 4) (Young 1905.01230)."""
    xa, scalar = _as_array(x)
    out = np.exp(-(xa**2) / 4.0)
    return out[0] if scalar else out


def transfer_rd(x):
    """Linear radiation-era transfer function at horizon crossing of r: T = W(x / sqrt 3)."""
    return tophat(np.asarray(x, dtype=float) / SQRT3)


WINDOWS = {"tophat": tophat, "gaussian": gaussian}


def tophat_envelope(x):
    """Upper bound on |W(x)| used only to trim integration ranges."""
    x = np.asarray(x, dtype=float)
    return np.where(x < 3.0, 1.0, 3.0 / x**2 * (1.0 + 1.0 / x))


def shell_envelope(x):
    x = np.asarray(x, dtype=float)
    return np.where(x < 1.0, 1.0, 1.0 / x)
