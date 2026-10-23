import numpy as np
import pytest

from ngpbh import windows
from ngpbh.moments import MomentsConfig, compaction_moments
from ngpbh.spectra import BrokenPowerLaw, DiracDelta, FlatBand, LogNormal


def test_window_limits_and_zero():
    assert windows.tophat(0.0) == pytest.approx(1.0)
    assert windows.shell(0.0) == pytest.approx(1.0)
    assert windows.transfer_rd(0.0) == pytest.approx(1.0)
    # series and closed form agree across the switch at x = 1e-2
    assert windows.tophat(0.00999) == pytest.approx(windows.tophat(0.01001), rel=1e-5)
    # first zero of W: tan x = x
    assert abs(windows.tophat(4.493409457909064)) < 1e-12


def test_dirac_moments_match_young_2022():
    # Young 2201.13345 eq. (29): delta spectrum, r = 2.74/k_p, no transfer function.
    # (His sigma_r^2 = 0.141 A is a typo for W_s = 0.141; W_s^2 = 0.0204.)
    m = compaction_moments(DiracDelta(1.0, 1.0), 2.74, MomentsConfig(transfer=False))
    assert m.sigma_c2 == pytest.approx(2.01, rel=2e-3)
    assert m.sigma_cr2 == pytest.approx(0.200, rel=2e-2)
    assert m.sigma_r2 == pytest.approx(0.141**2, rel=3e-2)
    assert m.gamma_cr == pytest.approx(1.0, abs=1e-12)
    assert m.gamma_peak == pytest.approx(1.0, abs=1e-12)
    assert m.rstar_over_r == pytest.approx(np.sqrt(3.0) / 2.74, rel=1e-12)


@pytest.mark.parametrize("transfer", [False, True])
def test_narrow_lognormal_approaches_dirac(transfer):
    cfg = MomentsConfig(transfer=transfer)
    d = compaction_moments(DiracDelta(1.0, 1.0), 2.7, cfg)
    ln = compaction_moments(LogNormal(1.0, 1.0, 0.001), 2.7, cfg)  # finite-width error ~ Delta^2
    for name in ("sigma_c2", "sigma_cr2", "sigma_r2", "s1_2", "s2_2", "scr1_2"):
        assert getattr(ln, name) == pytest.approx(getattr(d, name), rel=2e-4)


def test_moments_are_linear_in_amplitude():
    cfg = MomentsConfig()
    m1 = compaction_moments(LogNormal(1.0, 1.0, 0.7), 3.0, cfg)
    m3 = compaction_moments(LogNormal(3.0, 1.0, 0.7), 3.0, cfg)
    assert m3.sigma_c2 == pytest.approx(3.0 * m1.sigma_c2, rel=1e-12)
    assert m1.scaled(3.0).sigma_r2 == pytest.approx(m3.sigma_r2, rel=1e-12)


@pytest.mark.parametrize(
    "spec",
    [LogNormal(1.0, 1.0, 1.0), LogNormal(1.0, 1.0, 0.1), BrokenPowerLaw(1.0, 1.0, 4.0, 3.0, 1.0), FlatBand(1.0, 1e-3, 1.0)],
)
def test_moment_grids_are_converged(spec):
    base = compaction_moments(spec, 3.0, MomentsConfig(transfer=True))
    fine = compaction_moments(spec, 3.0, MomentsConfig(transfer=True, points_per_period=48, shape_points=80))
    for name in ("sigma_c2", "sigma_cr2", "sigma_r2", "s1_2", "s2_2", "scr1_2"):
        assert getattr(fine, name) == pytest.approx(getattr(base, name), rel=5e-6)


def test_lognormal_normalisations():
    a = LogNormal(2.0, 1.0, 0.5)
    b = LogNormal(2.0, 1.0, 0.5, normalization="peak")
    assert a.peak_value == pytest.approx(2.0 / (np.sqrt(2 * np.pi) * 0.5))
    assert b(1.0) == pytest.approx(2.0)
    lnk = np.linspace(-6, 6, 20001)
    assert np.trapezoid(a(np.exp(lnk)), lnk) == pytest.approx(2.0, rel=1e-8)
