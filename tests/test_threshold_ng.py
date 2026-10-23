import numpy as np
import pytest
from scipy.optimize import brentq

from ngpbh.nongaussian import Gaussian, Logarithmic, Quadratic
from ngpbh.spectra import DiracDelta
from ngpbh.threshold import MuscoThreshold, cl_threshold, delta_c_of_q


def test_delta_c_limits():
    assert delta_c_of_q(0.01) == pytest.approx(0.4, abs=3e-3)
    assert delta_c_of_q(1e4) == pytest.approx(2.0 / 3.0, abs=1e-3)
    q = np.geomspace(0.02, 100, 50)
    assert np.all(np.diff(delta_c_of_q(q)) > 0)


def test_linear_threshold_inverts_compaction():
    phi = 2.0 / 3.0
    cl = cl_threshold(0.55, phi)
    assert cl - cl**2 / (4 * phi) == pytest.approx(0.55)
    assert cl == pytest.approx(0.7756, abs=1e-4)  # pbhstat's g_c = 0.77


@pytest.mark.parametrize("transfer", [False, True])
def test_musco_dirac(transfer):
    res = MuscoThreshold(transfer=transfer).evaluate(DiracDelta(1.0, 1.0))
    x_m = brentq(lambda x: (x**2 - 1) * np.sin(x) + x * np.cos(x), 2.0, 3.5)
    assert res.rm_kpeak == pytest.approx(x_m, rel=1e-9)
    assert res.C_th == pytest.approx(0.5947, abs=5e-4)  # ~0.59 quoted in 2011.03014


@pytest.mark.parametrize("ng", [Quadratic(1.3), Quadratic(-0.7), Logarithmic(2.5), Logarithmic(-1.5)])
def test_template_derivatives(ng):
    z = np.linspace(-0.3, 0.3, 7)
    h = 1e-5
    assert np.allclose(ng.dF(z), (ng.F(z + h) - ng.F(z - h)) / (2 * h), rtol=1e-7)
    assert np.allclose(ng.d2F(z), (ng.dF(z + h) - ng.dF(z - h)) / (2 * h), rtol=1e-6)
    assert np.allclose(ng.d3F(z), (ng.d2F(z + h) - ng.d2F(z - h)) / (2 * h), rtol=1e-5, atol=1e-9)


def test_log_template_expands_to_fnl():
    ng = Logarithmic.from_fnl(0.8)
    z = 1e-3
    assert ng.f_nl == pytest.approx(0.8)
    assert ng.F(z) == pytest.approx(z + 0.6 * 0.8 * z**2, rel=1e-6)
    assert Logarithmic.from_beta(3.0).f_nl == pytest.approx(5 * 3.0 / 12)  # 2406.20089
    assert Logarithmic(2.5).bounds() == (-np.inf, 2.5)
    assert Logarithmic(-2.5).bounds() == (-2.5, np.inf)
    assert Gaussian().is_gaussian and Quadratic(0.0).is_gaussian
