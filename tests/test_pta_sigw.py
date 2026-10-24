"""Phase 5: scalar-induced GWs, thermal history, f_PBH conventions and the PTA likelihood."""

import numpy as np
import pytest
from scipy import integrate

from ngpbh import thermal
from ngpbh.abundance import AbundanceModel
from ngpbh.cosmology import HorizonMass
from ngpbh.moments import MomentsConfig
from ngpbh.pta import DEFAULT_DIR, DATA_SUBDIR, FreeSpectrum, credible_level, posterior_mass, sigma_1d, sigma_2d
from ngpbh.sigw import SIGWGrid, SIGWTemplate, grid_for, kernel, omega_rd, omega_rd_delta, omega_rd_ir
from ngpbh.spectra import BrokenPowerLaw, FlatBand, LogNormal
from ngpbh.stats import CollapseModel

HAVE_NG15 = (DEFAULT_DIR / DATA_SUBDIR / "density.npy").exists()


def test_flat_spectrum_gives_kohri_terada_constant():
    flat = FlatBand(A=1.0, k_min=1e-30, k_max=1e30)
    assert omega_rd(flat, [1.0], grid=SIGWGrid(t_max=1e12))[0] == pytest.approx(0.8222436, abs=2e-7)


def test_dirac_delta_closed_form():
    """Kohri & Terada: (3/64) ((4-x^2)/4)^2 x^2 (3x^2-2)^2 [pi^2 (3x^2-2)^2 Th(2 sqrt3 - 3x) + (4 + (3x^2-2) ln|1 - 4/(3x^2)|)^2]."""
    x = np.array([0.2, 0.7, 1.0, 1.1, 1.5, 1.9])
    b = 3 * x**2 - 2
    ref = 3 / 64 * ((4 - x**2) / 4) ** 2 * x**2 * b**2 * (
        np.pi**2 * b**2 * (2 * np.sqrt(3) > 3 * x) + (4 + b * np.log(np.abs(1 - 4 / (3 * x**2)))) ** 2)
    assert np.allclose(omega_rd_delta(x), ref, rtol=1e-12)


def test_kernel_corner_is_finite():
    """At t, s -> 1 the log bracket cancels to O(e^3); the kernel tends to 2 x projection^2 -> 0."""
    t = 1 + np.logspace(-9, -2, 8)
    k = kernel(t, 2 - t)
    assert np.all(np.isfinite(k)) and np.all(k >= 0) and k[0] < 1e-12


@pytest.mark.parametrize("spec", [LogNormal(1.0, 1.0, 1.0), BrokenPowerLaw(1.0, 1.0, 4.0, 3.0, 1.0)])
def test_quadrature_matches_dblquad(spec):
    x = 0.5
    f = lambda s, t: kernel(t, s) * spec(x * (t + s) / 2) * spec(x * (t - s) / 2) / 6.0
    t_hi = 2 * spec.support(1e-10)[1] / x
    edges = [1.0, np.sqrt(3)] + list(np.geomspace(2.0, t_hi, 16))
    ref = sum(integrate.dblquad(f, a, b, 0.0, 1.0, epsrel=1e-10)[0] for a, b in zip(edges[:-1], edges[1:]))
    assert omega_rd(spec, [x])[0] == pytest.approx(ref, rel=1e-8)


def test_grid_convergence_and_infrared_tail():
    spec = LogNormal(1.0, 1.0, 0.5)
    x = np.array([1e-3, 0.05, 0.5, 1.2, 3.0])
    g1 = grid_for(spec, x.min())
    g2 = SIGWGrid(t_max=10 * g1.t_max, h=g1.h / 2, n_gl=16, n_s=192)
    assert np.allclose(omega_rd(spec, x, grid=g1), omega_rd(spec, x, grid=g2), rtol=1e-9)
    assert omega_rd(spec, [1e-3])[0] / omega_rd_ir(spec, [1e-3])[0] == pytest.approx(1.0, abs=1e-4)


def test_template_scaling():
    tpl = SIGWTemplate(LogNormal(1.0, 1.0, 1.0), n_x=200)
    f = np.array([3e-9, 2e-8])
    a = tpl.h2_omega(f, A=0.1, k_star=1e8)
    b = tpl.h2_omega(f, A=0.2, k_star=1e8)
    assert np.allclose(b / a, 4.0)
    x = f / thermal.HZ_PER_INV_MPC / 1e8
    direct = thermal.H2_OMEGA_R * thermal.gw_dilution(f) * 0.01 * omega_rd(LogNormal(1.0, 1.0, 1.0), x)
    assert np.allclose(a, direct, rtol=1e-3)


def test_thermal_history_matches_table_frequencies():
    tab = thermal._TAB
    ratio = thermal.k_of_T(tab[:, 0]) * thermal.HZ_PER_INV_MPC / tab[:, 1]
    assert np.all(np.abs(ratio - 1) < 1e-3)
    assert thermal.horizon_mass(1.0) == pytest.approx(4.8e-2 * (thermal.g_rho(1.0) / 106.75) ** -0.5, rel=0.02)


def test_fpbh_conventions():
    m_h = np.array([1e-3, 1.0, 30.0])
    fivt, ent = HorizonMass(convention="fivt23"), HorizonMass(convention="entropy")
    ratio = np.exp(fivt.log_redshift_factor(m_h) - ent.log_redshift_factor(m_h))
    assert np.allclose(ratio, 12.74, rtol=2e-3)  # 1/(7.9e-10) vs entropy conservation (1.0e-8)
    for h in (fivt, ent):
        assert np.allclose(h.mh_of_r(h.r_of_mh(m_h)), m_h, rtol=1e-6)
    # 2306.17149 eqs. (4), (15) at T = 1 GeV
    T = 1.0
    k = 1.5e7 * (thermal.g_rho(T) / 106.75) ** 0.5 * (thermal.g_s(T) / 106.75) ** (-1 / 3) * T
    assert fivt.mh_of_r(1.0 / k) == pytest.approx(4.8e-2 * (thermal.g_rho(T) / 106.75) ** -0.5, rel=1e-4)
    # entropy conservation: f_PBH/beta = 3.76e8 (M_H/Msun)^-1/2 at g* = g*s = 106.75 (Carr et al. 2010 to 2%)
    w = ent.log_redshift_factor(thermal.horizon_mass(1e3)) + 0.5 * np.log(thermal.horizon_mass(1e3))
    g = thermal.g_rho(1e3)
    assert np.exp(w) * (g / 106.75) ** 0.25 == pytest.approx(3.76e8, rel=0.01)


def test_template_grid_cap_for_slow_uv_tails():
    """A BPL with beta = 0.25 has its support at k ~ 1e40 k*; capping t_max at 1e8 changes nothing."""
    spec = BrokenPowerLaw(1.0, 1.0, 4.0, 0.25, 1.0)
    x = np.array([1e-3, 0.1, 1.0])
    capped = omega_rd(spec, x, grid=grid_for(spec, 1e-4, rel=1e-10, t_max_cap=1e8))
    longer = omega_rd(spec, x, grid=grid_for(spec, 1e-4, rel=1e-10, t_max_cap=1e11))
    assert np.allclose(capped, longer, rtol=1e-8)
    assert SIGWTemplate.x_min * 1e10 <= 2e-9 / thermal.HZ_PER_INV_MPC  # lowest NG15 bin covered up to k* = 1e10


def test_amplitude_curve_matches_direct_solution():
    def model(k_star):
        return AbundanceModel(LogNormal(1.0, k_star, 1.0), CollapseModel(C_th=0.5, K=3.3, gamma=0.36),
                              moments_cfg=MomentsConfig(transfer=True, x_cut=4.49),
                              horizon=HorizonMass(convention="fivt23"), n_mh=40)
    curve = model(1e7).amplitude_curve([0.1, 30.0], f_targets=(1.0, 1e-3))
    assert curve[0, 1] == pytest.approx(model(3e8).solve_amplitude(1.0), rel=1e-5)
    assert curve[1, 0] == pytest.approx(model(1e6).solve_amplitude(1e-3), rel=1e-5)


@pytest.fixture
def synthetic_freespec(tmp_path):
    grid = np.linspace(-10.0, -4.0, 6001)
    centres = np.array([-6.5, -7.0])
    dens = -0.5 * ((grid[None, :] - centres[:, None]) / 0.1) ** 2 - np.log(0.1 * np.sqrt(2 * np.pi))
    np.save(tmp_path / "log10rhogrid.npy", grid)
    np.save(tmp_path / "density.npy", dens[None])
    np.save(tmp_path / "freqs.npy", np.array([2e-9, 4e-9]))
    return FreeSpectrum(tmp_path, n_bins=2), centres


def test_freespectrum_likelihood(synthetic_freespec):
    fs, centres = synthetic_freespec
    best = fs.loglike_rho(centres)
    assert best == pytest.approx(2 * (-np.log(0.1 * np.sqrt(2 * np.pi))) + 2 * np.log(fs._dx), rel=1e-6)
    assert fs.loglike_rho(centres + [0.1, 0.0]) == pytest.approx(best - 0.5, abs=1e-3)
    assert fs.loglike_rho(np.array([-30.0, -7.0])) == pytest.approx(fs.loglike_rho(np.array([-10.0, -7.0])))
    assert fs.loglike_rho(np.array([-3.0, -7.0])) == -np.inf
    # rho^2 = H100^2 h^2 Omega / (8 pi^4 f^5 T): Omega that puts bin 1 at its centre
    om = 10 ** (2 * centres) * 8 * np.pi**4 * fs.freqs**5 * fs.T / thermal.H100_HZ**2
    assert np.allclose(fs.log10_rho(om), centres)


def test_credible_level_on_gaussian():
    """For a 2D unit Gaussian, the HPD region first touching the half-plane x > d is the disc of radius d."""
    x = np.linspace(-8, 8, 800)  # no node on x = 0
    X, Y = np.meshgrid(x, x, indexing="ij")
    lp = -0.5 * (X**2 + Y**2)
    for d in (1.0, 2.0, 3.0):
        level = credible_level(lp, X > d)
        assert sigma_2d(level) == pytest.approx(d, abs=0.03)
    assert posterior_mass(lp, X > 0) == pytest.approx(0.5, abs=1e-3)
    assert sigma_1d(0.9545) == pytest.approx(2.0, abs=1e-3)


@pytest.mark.skipif(not HAVE_NG15, reason="needs data/ng15 (scripts/fetch_ng15.py)")
def test_ng15_powerlaw_amplitude():
    fs = FreeSpectrum(n_bins=14)
    la = np.linspace(-16, -13.5, 2501)
    lnl = fs.loglike_rho(fs.log10_rho_powerlaw(la, 13.0 / 3.0))
    p = np.exp(lnl - lnl.max())
    median = np.interp(0.5, np.cumsum(p) / p.sum(), la)
    assert 10**median == pytest.approx(2.4e-15, rel=0.1)  # NANOGrav 2306.16213
