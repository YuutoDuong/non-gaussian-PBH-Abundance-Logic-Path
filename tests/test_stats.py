from dataclasses import dataclass

import numpy as np
import pytest
from scipy import integrate

from ngpbh.moments import MomentsConfig, compaction_moments
from ngpbh.nongaussian import Gaussian, LocalNG, Logarithmic, Quadratic
from ngpbh.spectra import LogNormal
from ngpbh.stats import (
    CollapseModel,
    PeakGeometry,
    PeaksConfig,
    QuadratureConfig,
    bbks_f,
    log_beta,
    log_G,
)

COL = CollapseModel(C_th=0.56)


@pytest.fixture(scope="module")
def unit_mom():
    return compaction_moments(LogNormal(1.0, 1.0, 0.5), 2.7, MomentsConfig(transfer=True, x_cut=4.49))


@dataclass(frozen=True)
class Flip(LocalNG):
    """zeta = -zeta_g: statistically identical to a Gaussian zeta, but F' = -1."""

    def F(self, z):
        return -np.asarray(z, dtype=float)

    def dF(self, z):
        return -np.ones_like(np.asarray(z, dtype=float))

    def d2F(self, z):
        return np.zeros_like(np.asarray(z, dtype=float))

    def d3F(self, z):
        return np.zeros_like(np.asarray(z, dtype=float))


def dblquad_log_beta(mom, col, ng):
    sc, sr, rho = mom.sigma_c, mom.sigma_r, mom.gamma_cr
    lo_b, hi_b = ng.bounds()
    z_lo = max(-16.0, lo_b / sr * (1 - 1e-9))
    z_hi = min(16.0, hi_b / sr * (1 - 1e-9))

    def f(cl, z):
        dF = float(ng.dF(sr * z))
        if dF == 0.0:
            return 0.0
        nu = cl / (sc * dF)
        lp = -(nu**2 - 2 * rho * nu * z + z**2) / (2 * (1 - rho**2)) - np.log(2 * np.pi * np.sqrt(1 - rho**2))
        C = cl - cl**2 / (4 * col.phi)
        if C <= col.C_th:
            return 0.0
        return np.exp(lp) / (sc * abs(dF)) * col.K * (C - col.C_th) ** col.gamma

    val, _ = integrate.dblquad(f, z_lo, z_hi, col.cl_lower, col.cl_upper, epsabs=0, epsrel=1e-9)
    return np.log(val)


def test_bbks_f_and_G():
    assert bbks_f(np.array([0.19999]))[0] == pytest.approx(bbks_f(np.array([0.20001]))[0], rel=1e-3)
    x = np.array([8.0, 12.0])
    assert np.allclose(bbks_f(x), x**3 - 3 * x, rtol=1e-12)
    assert np.exp(log_G(8.0, 0.5))[0] == pytest.approx(8**3 - 3 * 8 + 3 * 8 * 0.25, rel=1e-10)
    assert np.exp(log_G(1.5, 1e-3))[0] == pytest.approx(bbks_f(np.array([1.5]))[0], rel=1e-5)


@pytest.mark.parametrize("statistic", ["threshold", "peaks"])
def test_gaussian_limits(unit_mom, statistic):
    mom = unit_mom.scaled(0.02)
    ref = log_beta(mom, COL, None, statistic)
    # peaks: the NG path marginalises G(a nu + b z, s) over z numerically, ~1e-6 accurate
    tol = 1e-6 if statistic == "threshold" else 1e-5
    for ng in (Quadratic(0.0), Quadratic(1e-10), Logarithmic(1e10), Logarithmic(-1e10), Flip()):
        assert log_beta(mom, COL, ng, statistic) == pytest.approx(ref, abs=tol), ng


@pytest.mark.parametrize(
    "ng", [Quadratic(1.0), Quadratic(-1.0), Quadratic(-3.0), Logarithmic(2.5), Logarithmic(0.8), Logarithmic(-2.5)]
)
def test_threshold_matches_dblquad(unit_mom, ng):
    mom = unit_mom.scaled(0.03)
    assert log_beta(mom, COL, ng) == pytest.approx(dblquad_log_beta(mom, COL, ng), abs=2e-6)


@pytest.mark.parametrize("statistic", ["threshold", "peaks"])
@pytest.mark.parametrize("ng", [Quadratic(2.0), Quadratic(-2.0), Logarithmic(1.0)])
def test_grid_doubling_is_stable(unit_mom, statistic, ng):
    mom = unit_mom.scaled(0.02)
    base = log_beta(mom, COL, ng, statistic)
    fine = log_beta(mom, COL, ng, statistic, quad=QuadratureConfig().refined())
    assert fine == pytest.approx(base, abs=1e-5)


def test_threshold_monte_carlo(unit_mom):
    """Direct sampling of (C_g, zeta_r) at large variance, where beta ~ 1e-3."""
    mom = unit_mom.scaled(0.15)
    ng = Quadratic(1.0)
    rng = np.random.default_rng(1)
    n = 4_000_000
    rho = mom.gamma_cr
    z = rng.standard_normal(n)
    nu = rho * z + np.sqrt(1 - rho**2) * rng.standard_normal(n)
    c_l = ng.dF(mom.sigma_r * z) * mom.sigma_c * nu
    ok = (c_l > COL.cl_lower) & (c_l < COL.cl_upper)
    w = np.zeros(n)
    w[ok] = np.exp(COL.log_mass_ratio(c_l[ok]))
    beta_mc, err = w.mean(), w.std() / np.sqrt(n)
    beta = np.exp(log_beta(mom, COL, ng))
    assert abs(beta_mc - beta) < 4 * err
    assert err / beta < 0.03


def test_flip_symmetry_counts_troughs(unit_mom):
    """With F' < 0 a peak of C_l is a trough of C_g; Flip must reproduce the Gaussian peaks."""
    mom = unit_mom.scaled(0.02)
    geo = PeakGeometry.from_moments(mom, PeaksConfig())
    g = log_beta(mom, COL, None, "peaks", geometry=geo)
    assert log_beta(mom, COL, Flip(), "peaks", geometry=geo) == pytest.approx(g, abs=1e-5)


def test_high_peak_limit(unit_mom):
    mom = unit_mom.scaled(0.006)  # nu ~ 10
    full = log_beta(mom, COL, None, "peaks", PeaksConfig(density="bbks"))
    high = log_beta(mom, COL, None, "peaks", PeaksConfig(density="high_peak"))
    assert abs(full - high) < 0.05
    assert full > log_beta(mom, COL, None, "threshold")


def test_linear_estimator_press_schechter_dirac():
    """pbhstat-style PS: beta = 2 int_{g_c}^{4/3} K (g - g_c)^gamma N(g; sigma) dg."""
    from ngpbh.spectra import DiracDelta

    mom = compaction_moments(DiracDelta(0.03, 1.0), 2.74, MomentsConfig(transfer=False))
    col = CollapseModel(C_th=0.55, K=4.0, gamma=0.36, estimator="linear", C_l_th=0.77)
    s = mom.sigma_c
    ref, _ = integrate.quad(
        lambda g: 2 * 4.0 * (g - 0.77) ** 0.36 * np.exp(-0.5 * (g / s) ** 2) / (np.sqrt(2 * np.pi) * s),
        0.77,
        4.0 / 3.0,
        epsabs=0,
        epsrel=1e-12,
        limit=200,
    )
    assert log_beta(mom, col, quad=QuadratureConfig(ps_factor=2.0)) == pytest.approx(np.log(ref), abs=1e-7)
