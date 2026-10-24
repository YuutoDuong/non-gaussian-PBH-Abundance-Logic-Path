"""Primordial-NG corrections to the SIGW spectrum (ngpbh.sigw_ng)."""

import numpy as np
import pytest

from ngpbh.nongaussian import Curvaton, Gaussian, Logarithmic, Quadratic
from ngpbh.pta import DATA_SUBDIR, DEFAULT_DIR, FreeSpectrum, log_posterior_sigw, log_posterior_sigw_ng
from ngpbh.sigw import SIGWTemplate, grid_for, kernel, omega_rd
from ngpbh.sigw_ng import (NGSIGWTemplate, disconnected_terms, kernel_split, local_coefficients, mc_term,
                           one_loop_shape, vertex)
from ngpbh.spectra import BrokenPowerLaw, LogNormal

HAVE_NG15 = (DEFAULT_DIR / DATA_SUBDIR / "density.npy").exists()


def test_kernel_split_reproduces_kernel():
    rng = np.random.default_rng(3)
    near = np.sqrt(3) + rng.choice([-1, 1], 400) * 10 ** rng.uniform(-9, -1, 400)  # around the resonance
    t = np.concatenate([1 + 10 ** rng.uniform(-6, 5, 4000), near])
    s = rng.uniform(0, 1, t.size)
    a, ia, ib = kernel_split(t, s)
    assert np.allclose(a**2 * (ia**2 + ib**2) / 2, kernel(t, s), rtol=1e-11)


def test_vertex_matches_kernel_variables():
    """|z|^2 = uv a and the transfer functions of a vertex from its momentum vector."""
    rng = np.random.default_rng(4)
    q = rng.normal(size=(3000, 3)) * 10 ** rng.uniform(-2, 2, (3000, 1))
    k = 0.7
    zz, ia, ib = vertex(q, k)
    v = np.linalg.norm(q, axis=1) / k
    u = np.linalg.norm(np.array([0.0, 0.0, k]) - q, axis=1) / k
    a, ia2, ib2 = kernel_split(u + v, u - v)
    assert np.allclose(np.abs(zz), u * v * a, rtol=1e-8)
    assert np.allclose(ia, ia2, rtol=1e-8, atol=1e-14) and np.allclose(ib, ib2, rtol=1e-8, atol=1e-14)


@pytest.mark.parametrize("spec", [LogNormal(1.0, 1.0, 0.5), BrokenPowerLaw(1.0, 1.0, 4.0, 3.0, 1.0)])
def test_monte_carlo_gaussian_term_matches_quadrature(spec):
    x = np.array([0.01, 0.5, 1.5])
    det = omega_rd(spec, x, grid=grid_for(spec, x.min()))
    for xi, d in zip(x, det):
        m, e = mc_term(spec, "gaussian", xi, n_log2=12, n_rand=8)
        assert m == pytest.approx(d, rel=max(5 * e / d, 0.01))


def test_one_loop_spectrum_infrared_limit():
    """s2 -> 2 x^3 int dp s^2/p^4 = x^3 exp(9 D^2/4)/(sqrt(pi) D) for a lognormal of width D."""
    for width in (0.5, 1.0):
        s2 = one_loop_shape(LogNormal(1.0, 1.0, width))
        x = 1e-3
        assert s2.shape(x) == pytest.approx(x**3 * np.exp(9 * width**2 / 4) / (np.sqrt(np.pi) * width), rel=1e-3)


def test_disconnected_terms_match_monte_carlo():
    """hybrid = 2 Omega_bil[s, s2] and reducible = Omega_G[s2] against the 3- and 4-momentum MC."""
    spec = LogNormal(1.0, 1.0, 1.0)
    x = np.array([0.03, 1.0])
    det = disconnected_terms(spec, x)
    for term in ("hybrid", "reducible"):
        for xi, d in zip(x, det[term]):
            m, e = mc_term(spec, term, xi, n_log2=12, n_rand=8)
            assert m == pytest.approx(d, rel=max(5 * e / d, 0.015))


# Adshead, Lozanov & Weiner 2105.01659 Fig. 4, lognormal sigma = 1/5, A = F_NL = 1, read from the
# figure's vector paths (scripts/validate_sigw_ng.py): (term, x, value)
ALW_FIG4 = [("C", 0.45, 0.672), ("C", 1.0, 0.395), ("Z", 0.6, 0.053), ("Z", 1.0, -0.041), ("Z", 1.7, 0.136),
            ("planar", 0.6, 1.134), ("nonplanar", 0.6, -0.151), ("nonplanar", 1.7, -0.141)]


@pytest.mark.parametrize("term,x,ref", ALW_FIG4)
def test_connected_terms_match_adshead_et_al(term, x, ref):
    m, e = mc_term(LogNormal(1.0, 1.0, 0.2), term, x, n_log2=13, n_rand=8)
    scale = {"C": 0.67, "Z": 0.14, "planar": 1.13, "nonplanar": 0.15}[term]  # curve maxima
    assert abs(m - ref) < 0.04 * scale + 3 * e


def test_local_coefficients():
    assert local_coefficients(Gaussian()) == (0.0, 0.0)
    assert local_coefficients(Quadratic(-2.0)) == pytest.approx((-1.2, 0.0))
    mu = 2.0 / 3.0
    assert local_coefficients(Logarithmic(mu)) == pytest.approx((1 / (2 * mu), 1 / (3 * mu**2)), rel=1e-12)
    r = 0.9  # sudden-decay curvaton: f_NL and g_NL of Sasaki, Valiviita & Wands (astro-ph/0607627)
    f_nl = 5 / (4 * r) - 5 / 3 - 5 * r / 6
    g_nl = 25 / 54 * (-9 / r + 0.5 + 10 * r + 3 * r**2)
    F, G = local_coefficients(Curvaton(r))
    assert F == pytest.approx(0.6 * f_nl, rel=1e-6) and G == pytest.approx(0.36 * g_nl, rel=1e-4)


def _toy_template():
    t2 = SIGWTemplate(LogNormal(1.0, 1.0, 1.0), n_x=200)
    x = np.geomspace(1e-4, 30.0, 40)
    return NGSIGWTemplate(t2, x, 2.0 + 0.0 * x, 1.5 + 0.0 * x, 1.0)


def test_ng_template_series():
    tm = _toy_template()
    x = np.array([0.1, 1.0, 3.0])
    t2 = tm.t2(x)
    assert np.allclose(tm.omega(x, 0.1), 0.01 * t2)
    om = tm.omega(x, 0.1, F=1.2, G=0.5, order=4)
    assert np.allclose(om, t2 * (0.01 + 1e-3 * (1.44 * 2.0 + 12 * 0.5) + 1.2**4 * 1e-4 * 1.5))
    assert np.allclose(tm.omega(x, 0.1, F=1.2, order=3), t2 * (0.01 + 1e-3 * 1.44 * 2.0))


@pytest.mark.skipif(not HAVE_NG15, reason="NG15 free spectrum not downloaded (scripts/fetch_ng15.py)")
def test_ng_posterior_reduces_to_gaussian():
    tm = _toy_template()
    fs = FreeSpectrum(n_bins=14, smooth=0.05)
    lk, la = np.linspace(6.0, 9.0, 13), np.linspace(-2.5, 0.0, 26)
    assert np.allclose(log_posterior_sigw_ng(tm, fs, lk, la), log_posterior_sigw(tm.t2, fs, lk, la))
    # positive F^2 A^3 terms raise Omega: the posterior moves to smaller A
    lp = log_posterior_sigw_ng(tm, fs, lk, la, F=1.2)
    p0, p1 = np.exp(log_posterior_sigw(tm.t2, fs, lk, la)), np.exp(lp)
    assert (p1.sum(0) @ la) / p1.sum() < (p0.sum(0) @ la) / p0.sum()
