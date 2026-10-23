"""Phase 3b: maxima of the non-Gaussian compaction C_l = F'(zeta_r) C_g (docs/phase3b_derivation.tex)."""

from dataclasses import dataclass

import numpy as np
import pytest
from scipy.integrate import quad

from ngpbh import windows
from ngpbh.moments import MomentsConfig, compaction_moments
from ngpbh.nongaussian import Curvaton, LocalNG, Logarithmic, Quadratic
from ngpbh.spectra import LogNormal
from ngpbh.stats import (
    CollapseModel,
    NGPeakGeometry,
    PeaksConfig,
    QuadratureConfig,
    _table_2d,
    bbks_f,
    bbks_shape_moments,
    log_beta,
    log_G,
    log_G2,
    make_peak_geometry,
)
from ngpbh.validation import kac_rice_mc, semi_analytic_weight

CFG = MomentsConfig(transfer=True, x_cut=4.49)
C3A = PeaksConfig(constraint="C_g")
C3B = PeaksConfig(constraint="C_l")


def moments(width, rk=2.74, sigma_c=0.1):
    unit = compaction_moments(LogNormal(1.0, 1e7, width), rk / 1e7, CFG)
    return unit.scaled(sigma_c**2 / unit.sigma_c2)


@dataclass(frozen=True)
class Linear(LocalNG):
    """F = a zeta_g: non-Gaussian in name only (F'' = 0)."""

    a: float

    def F(self, z):
        return self.a * np.asarray(z, dtype=float)

    def dF(self, z):
        return np.full_like(np.asarray(z, dtype=float), self.a)

    def d2F(self, z):
        return np.zeros_like(np.asarray(z, dtype=float))

    def d3F(self, z):
        return np.zeros_like(np.asarray(z, dtype=float))


def test_shape_quadrature_reproduces_bbks_f():
    x = np.array([0.3, 0.7, 1.0, 2.0, 3.0, 5.0, 8.0, 12.0])
    f, f2 = bbks_shape_moments(x)
    assert np.allclose(f, bbks_f(x), rtol=1e-9, atol=0)
    # large x: every configuration is a peak, f2 = x^2/3 - E[3y^2 + z^2]/3 = x^2/3 - 1/3
    assert f2[-1] == pytest.approx(12.0**2 / 3 - 1 / 3, rel=1e-9)
    # first-order coefficient exceeds the spherical-peak estimate 9/x (Jensen)
    assert np.all(27 * f2 / f > 9 / x)


def test_table_2d_against_direct_quadrature():
    rng = np.random.default_rng(0)
    m = rng.uniform(1.0, 40.0, 3000)
    s = np.exp(rng.uniform(np.log(0.01), 0.0, 3000))
    log_g, r2 = _table_2d().lookup(m, s)
    assert np.max(np.abs(log_g - log_G(m, s))) < 5e-4
    r2_direct = 27.0 * np.exp(log_G2(m, s) - log_G(m, s))
    assert np.max(np.abs(r2 / r2_direct - 1)) < 2e-3


def test_new_moments_against_quad():
    spec = LogNormal(1.0, 1e7, 0.5)
    r = 2.74 / 1e7
    mom = compaction_moments(spec, r, CFG)
    phi = 2.0 / 3.0

    def integral(fn):
        g = lambda lnk: fn(np.exp(lnk) * r) * spec(np.exp(lnk)) * windows.transfer_rd(np.exp(lnk) * r) ** 2
        return quad(g, np.log(spec.k_star) - 12 * 0.5, np.log(4.49 / r), limit=400)[0]

    W, Ws = windows.tophat, windows.shell
    assert mom.s1r_2 == pytest.approx(integral(lambda x: x**2 * Ws(x) ** 2), rel=1e-6)
    assert mom.s2cr_2 == pytest.approx(integral(lambda x: 2 * phi / 3 * x**6 * W(x) * Ws(x)), rel=1e-6)
    assert mom.s2r_2 == pytest.approx(integral(lambda x: x**4 * Ws(x) ** 2), rel=1e-6)
    # Cauchy-Schwarz between C_g and zeta_r at each derivative order
    assert mom.scr1_2**2 < mom.s1_2 * mom.s1r_2
    assert mom.s2cr_2**2 < mom.s2_2 * mom.s2r_2


def test_gaussian_limit_is_phase_3a():
    mom = moments(1.0)
    col = CollapseModel(C_th=0.55)
    assert log_beta(mom, col, None, "peaks", C3B) == log_beta(mom, col, None, "peaks", C3A)


def test_linear_template_is_phase_3a():
    """F'' = 0: lam = 0 and kappa = 0 at every node, so 3b counts the same maxima as 3a
    (up to the 1D vs 2D table interpolation of G)."""
    mom = moments(1.0)
    col = CollapseModel(C_th=0.55)
    ng = Linear(1.3)
    lb_a = log_beta(mom, col, ng, "peaks", C3A)
    lb_b = log_beta(mom, col, ng, "peaks", C3B)
    assert lb_b == pytest.approx(lb_a, abs=5e-4)


@pytest.mark.parametrize(
    "width, ng, c",
    [(0.3, Quadratic(-2.0), 0.6), (1.0, Curvaton(0.3), 0.45), (1.0, Logarithmic(5 / 6), 0.6), (1.0, Quadratic(2.0), 0.6)],
)
def test_kac_rice_monte_carlo(width, ng, c):
    """The 20-variable Gaussian sampled independently of the derivation, with the exact Hessian
    of C_l: agrees with eq. (17) times the first-order rank-one factor."""
    mom = moments(width)
    z = mom.gamma_cr * (c / mom.sigma_c) * mom.sigma_r
    exact, no_r1, e_exact, e_no_r1 = kac_rice_mc(mom, ng, c, z, n_samples=400_000, seed=11)
    assert no_r1 == pytest.approx(semi_analytic_weight(mom, ng, c, z, rank_one=False), abs=4 * e_no_r1 + 5e-4 * no_r1)
    assert exact == pytest.approx(semi_analytic_weight(mom, ng, c, z, rank_one=True), abs=4 * e_exact + 5e-4 * exact)


def test_logarithmic_template_has_no_rank_one_term():
    mom = moments(1.0)
    geo = NGPeakGeometry.from_moments(mom, C3B)
    ng = Logarithmic(5 / 6)
    z = np.linspace(-0.3, 0.5, 9)
    c_l = 0.6
    args = (np.array(c_l / (mom.sigma_c * ng.dF(z))), z / mom.sigma_r, c_l, ng.dF(z), ng.d2F(z), ng.d3F(z), 1.0)
    assert np.allclose(geo.log_weight_ng(*args, rank_one=True), geo.log_weight_ng(*args, rank_one=False), atol=1e-12)


def test_narrow_spectra_approach_phase_3a():
    """For a monochromatic spectrum C_g, zeta_r and Y are proportional fields, so 3b = 3a."""
    col = CollapseModel(C_th=0.55)
    ng = Quadratic(1.0)
    diff = {}
    for width in (0.05, 1.0):
        mom = moments(width)
        diff[width] = abs(log_beta(mom, col, ng, "peaks", C3B) - log_beta(mom, col, ng, "peaks", C3A))
    assert diff[0.05] < 0.1 * diff[1.0]


def test_phase_3b_grid_refinement():
    mom = moments(1.0)
    col = CollapseModel(C_th=0.55)
    for ng in (Quadratic(-2.0), Logarithmic(5 / 6)):
        lb = log_beta(mom, col, ng, "peaks", C3B)
        lb_fine = log_beta(mom, col, ng, "peaks", C3B, QuadratureConfig().refined())
        assert lb_fine == pytest.approx(lb, abs=1e-4)


def test_geometry_amplitude_independent():
    """NGPeakGeometry built from unit moments and from scaled moments gives the same weights."""
    unit = compaction_moments(LogNormal(1.0, 1e7, 1.0), 2.74 / 1e7, CFG)
    A = 0.01 / unit.sigma_c2
    col = CollapseModel(C_th=0.55)
    ng = Quadratic(-1.0)
    lb_unit = log_beta(unit.scaled(A), col, ng, "peaks", C3B, geometry=make_peak_geometry(unit, C3B))
    lb_scaled = log_beta(unit.scaled(A), col, ng, "peaks", C3B)
    assert lb_unit == pytest.approx(lb_scaled, abs=1e-10)
