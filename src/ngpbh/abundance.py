"""f_PBH from beta(M_H), and the amplitude A_* that gives a target abundance.

    f_PBH = (1/Omega_DM) int dln M_H (M_eq/M_H)^(1/2) beta(M_H)       (2211.01728 eq. 61)

beta(M_H) is evaluated with all variances smoothed at the scale r(M_H) that crosses the horizon
when the horizon mass is M_H, and the transfer function taken at that time (tau = r).
"""

from dataclasses import dataclass, field
from typing import Optional

import numpy as np
from scipy.optimize import brentq
from scipy.special import logsumexp

from .cosmology import HorizonMass
from .moments import MomentsConfig, compaction_moments
from .nongaussian import LocalNG
from .spectra import PowerSpectrum
from .stats import CollapseModel, PeaksConfig, QuadratureConfig, log_beta, make_peak_geometry


@dataclass
class AbundanceModel:
    spectrum: PowerSpectrum  # only its shape is used; the amplitude is the free parameter
    collapse: CollapseModel
    ng: Optional[LocalNG] = None
    statistic: str = "threshold"
    moments_cfg: MomentsConfig = field(default_factory=MomentsConfig)
    horizon: HorizonMass = field(default_factory=HorizonMass)
    peaks: PeaksConfig = field(default_factory=PeaksConfig)
    quad: QuadratureConfig = field(default_factory=QuadratureConfig)
    mh_range: Optional[tuple] = None  # (M_H min, M_H max) [Msun]; automatic if None
    n_mh: int = 60
    rel_sigma: float = 0.05  # automatic range keeps sigma_c^2(r) >= rel_sigma * max

    def __post_init__(self):
        self._unit = self.spectrum.with_amplitude(1.0)
        if self.mh_range is not None:
            m_h = np.geomspace(*self.mh_range, self.n_mh)
            self.r = self.horizon.r_of_mh(m_h)
        else:
            self.r = self._auto_r_grid()
        self.m_h = self.horizon.mh_of_r(self.r)
        self.unit_moments = [compaction_moments(self._unit, float(r), self.moments_cfg) for r in self.r]
        self.geometries = None
        if self.statistic == "peaks":
            self.geometries = [make_peak_geometry(m, self.peaks) for m in self.unit_moments]
        ln_mh = np.log(self.m_h)
        w = np.zeros_like(ln_mh)
        w[1:] += 0.5 * np.diff(ln_mh)
        w[:-1] += 0.5 * np.diff(ln_mh)
        self._log_w = np.log(w) + self.horizon.log_redshift_factor(self.m_h)

    def _auto_r_grid(self):
        kp = self.spectrum.k_peak
        k_lo, k_hi = self._unit.support(1e-6) if not self._unit.is_delta else (kp, kp)
        r_coarse = np.geomspace(0.05 / k_hi, 200.0 / k_lo, 160)
        s2 = np.array([compaction_moments(self._unit, r, self.moments_cfg).sigma_c2 for r in r_coarse])
        good = np.where(s2 >= self.rel_sigma * s2.max())[0]
        lo = r_coarse[max(good[0] - 1, 0)]
        hi = r_coarse[min(good[-1] + 1, len(r_coarse) - 1)]
        return np.geomspace(lo, hi, self.n_mh)

    def log_beta_grid(self, A):
        geos = self.geometries or [None] * len(self.r)
        return np.array(
            [
                log_beta(m.scaled(A), self.collapse, self.ng, self.statistic, self.peaks, self.quad, g)
                for m, g in zip(self.unit_moments, geos)
            ]
        )

    def log_fpbh(self, A):
        return float(logsumexp(self.log_beta_grid(A) + self._log_w))

    def fpbh(self, A):
        return float(np.exp(self.log_fpbh(A)))

    def solve_amplitude(self, f_target=1.0, bracket=(1e-4, 1.0), xtol=1e-6):
        """A_* with f_PBH(A_*) = f_target, by root finding in ln A."""
        target = np.log(f_target)
        g = lambda lnA: max(self.log_fpbh(np.exp(lnA)), -1e300) - target  # brentq needs finite values
        lo, hi = np.log(bracket[0]), np.log(bracket[1])
        g_lo, g_hi = g(lo), g(hi)
        if not (g_lo < 0 < g_hi):
            raise RuntimeError(f"target not bracketed: log f_PBH = {g_lo + target:.3g}, {g_hi + target:.3g}")
        return float(np.exp(brentq(g, lo, hi, xtol=xtol)))
