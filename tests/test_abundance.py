import numpy as np
import pytest
from scipy.optimize import brentq

from ngpbh.abundance import AbundanceModel
from ngpbh.cosmology import HorizonMass
from ngpbh.moments import MomentsConfig, compaction_moments
from ngpbh.nongaussian import Quadratic
from ngpbh.spectra import FlatBand, LogNormal
from ngpbh.stats import CollapseModel, QuadratureConfig, log_beta


def test_horizon_mass_relations_invert():
    for conv in ("ferrante", "pbhstat"):
        hz = HorizonMass(conv)
        m = np.array([1e-15, 1.0, 1e5])
        assert np.allclose(hz.mh_of_r(hz.r_of_mh(m)), m, rtol=1e-12)
    # 2211.01728 eq. (39): r_m = 2.4e-7 Mpc (g*/10.75)^(1/12) (M_H/Msun)^(1/2)
    hz = HorizonMass("ferrante", g_star=10.75)
    assert float(hz.r_of_mh(1.0)) == pytest.approx(2.4e-7, rel=0.02)


def test_ferrante_broad_spectrum_benchmark():
    """2211.01728 Fig. 8: flat P = A on [1e6, 1e15] Mpc^-1, non-linearities only, C_th = 0.56,
    K = 4.36: A = 9.2e-3 gives beta = 1e-16 at M_H = 1e-15 Msun."""
    hz = HorizonMass("ferrante", g_star=106.75)
    r = float(hz.r_of_mh(1e-15))
    mom = compaction_moments(FlatBand(1.0, 1e6, 1e15), r, MomentsConfig(transfer=True))
    col = CollapseModel(C_th=0.56, K=4.36, gamma=0.36)
    a = np.exp(brentq(lambda la: log_beta(mom.scaled(np.exp(la)), col) - np.log(1e-16), np.log(1e-3), np.log(0.1)))
    assert a == pytest.approx(9.2e-3, rel=0.02)


@pytest.fixture(scope="module")
def gaussian_model():
    return AbundanceModel(LogNormal(1.0, 1e7, 0.5), CollapseModel(C_th=0.56), n_mh=40)


def test_fpbh_increases_with_amplitude(gaussian_model):
    vals = [gaussian_model.log_fpbh(a) for a in (0.015, 0.02, 0.025, 0.03)]
    assert np.all(np.diff(vals) > 0)


def test_solve_amplitude(gaussian_model):
    a = gaussian_model.solve_amplitude(1e-2)
    assert gaussian_model.fpbh(a) == pytest.approx(1e-2, rel=1e-4)


@pytest.mark.parametrize("ng", [None, Quadratic(1.0)])
def test_amplitude_converged_under_grid_doubling(ng):
    kw = dict(ng=ng, mh_range=None)
    base = AbundanceModel(LogNormal(1.0, 1e7, 0.5), CollapseModel(C_th=0.56), n_mh=40, **kw)
    fine = AbundanceModel(
        LogNormal(1.0, 1e7, 0.5),
        CollapseModel(C_th=0.56),
        n_mh=80,
        quad=QuadratureConfig().refined(),
        moments_cfg=MomentsConfig(points_per_period=48, shape_points=80),
        **kw,
    )
    assert fine.solve_amplitude(1.0) == pytest.approx(base.solve_amplitude(1.0), rel=2e-4)


@pytest.mark.parametrize("r_dec, a_ref", [(0.5, 9.5e-3), (0.1, 2.8e-3)])
def test_ferrante_curvaton_benchmark(r_dec, a_ref):
    """2211.01728 Fig. 8 (right): curvaton NG, same setup as the Gaussian benchmark."""
    from ngpbh.nongaussian import Curvaton

    hz = HorizonMass("ferrante", g_star=106.75)
    mom = compaction_moments(FlatBand(1.0, 1e6, 1e15), float(hz.r_of_mh(1e-15)), MomentsConfig(transfer=True))
    col = CollapseModel(C_th=0.56, K=4.36, gamma=0.36)
    ng = Curvaton(r_dec)
    a = np.exp(brentq(lambda la: log_beta(mom.scaled(np.exp(la)), col, ng) - np.log(1e-16), np.log(1e-4), np.log(0.1)))
    assert a == pytest.approx(a_ref, rel=0.02)
