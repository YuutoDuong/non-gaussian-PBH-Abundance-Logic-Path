"""Horizon-mass relations and the constants that convert beta(M_H) into f_PBH.

Four conventions; the comoving scale r crosses the horizon when k = 1/r in all of them.

* "ferrante": M_H = 17 Msun (g*/10.75)^(-1/6) (k_H / 1e6 Mpc^-1)^(-2) with k_H = 1/r,
  i.e. r_m = 2.4e-7 Mpc (g*/10.75)^(1/12) (M_H/Msun)^(1/2)  (2211.01728 eq. 39), and
  f_PBH = (1/Omega_DM) int dln M_H (M_eq/M_H)^(1/2) beta (eq. 61).
* "pbhstat": R = (M_H / M_eq)^(1/2) / k_eq with k_eq = 0.01 (Omega_m/0.31) Mpc^-1
  (2509.14326 eq. 4), same f_PBH formula.
* "fivt23" (2306.17149 eqs. 4, 14, 15), with g*(T) and g*s(T) of the Standard Model:
  k = 1.5e7 Mpc^-1 (g*/106.75)^(1/2) (g*s/106.75)^(-1/3) T/GeV,
  M_H = 4.8e-2 Msun (g*/106.75)^(-1/2) (T/GeV)^-2,
  f_PBH = (1/Omega_DM) int dln M_H (M_H/Msun)^(-1/2) (g*/106.75)^(3/4) (g*s/106.75)^(-1) beta/7.9e-10.
* "entropy": k = aH and M_H = 4 pi M_P^2/H exactly (`thermal`), and the PBH density redshifted
  with entropy conservation, rho_PBH,0 = (3/4) beta (g*/g*s) T s_0, so that
  f_PBH = (1/Omega_DM) int dln M_H (M_H/Msun)^(-1/2) (g*/106.75)^(3/4) (g*s/106.75)^(-1) beta/1.0e-8.
  This agrees with Carr et al. (0912.5297) to 2%. For the same beta, "ferrante" gives a 5.3 times
  larger f_PBH and "fivt23" a 12.7 times larger one (at g* = 106.75), i.e. they lower the
  amplitude A_* for f_PBH = 1 by about ln(5.3)/|ln beta| ~ 8% and ln(12.7)/|ln beta| ~ 12%.
"""

from dataclasses import dataclass
from functools import lru_cache

import numpy as np

from . import thermal

M_EQ = 2.8e17  # horizon mass at matter-radiation equality [Msun]
OMEGA_DM = 0.12 / 0.674**2  # 2211.01728, 2306.17149
CONVENTIONS = ("ferrante", "pbhstat", "fivt23", "entropy")


@lru_cache(maxsize=None)
def _thermal_table(convention, omega_dm):
    """ln k(T), ln M_H(T) and ln of the f_PBH weight per unit beta, on a grid in T (decreasing M_H)."""
    T = np.geomspace(1e-5, 1e6, 6001)
    gr, gs = thermal.g_rho(T), thermal.g_s(T)
    if convention == "fivt23":
        k = 1.5e7 * (gr / 106.75) ** 0.5 * (gs / 106.75) ** (-1.0 / 3.0) * T
        m_h = 4.8e-2 * (gr / 106.75) ** -0.5 * T**-2
        log_w = -0.5 * np.log(m_h) + 0.75 * np.log(gr / 106.75) - np.log(gs / 106.75) - np.log(7.9e-10 * omega_dm)
    else:
        k = thermal.k_of_T(T)
        m_h = thermal.horizon_mass(T)
        s0 = 2.0 * np.pi**2 / 45.0 * thermal.G_S0 * thermal.T_0**3
        rho_c_h2 = 3.0 * thermal.M_P**2 * (thermal.H100_HZ * thermal.HBARC_GEV_M / 299792458.0) ** 2
        log_w = np.log(0.75 * gr / gs * T * s0 / (rho_c_h2 * omega_dm * 0.674**2))
    return np.log(k), np.log(m_h), log_w


@dataclass(frozen=True)
class HorizonMass:
    convention: str = "ferrante"
    g_star: float = 106.75
    omega_m: float = 0.31
    m_eq: float = M_EQ
    omega_dm: float = OMEGA_DM

    def __post_init__(self):
        if self.convention not in CONVENTIONS:
            raise ValueError(f"convention must be one of {CONVENTIONS}")

    @property
    def _table(self):
        return _thermal_table(self.convention, self.omega_dm)

    def r_of_mh(self, m_h):
        """Comoving scale [Mpc] that crosses the horizon when the horizon mass is m_h [Msun]."""
        m_h = np.asarray(m_h, dtype=float)
        if self.convention == "ferrante":
            m17 = 17.0 * (self.g_star / 10.75) ** (-1.0 / 6.0)
            return np.sqrt(m_h / m17) / 1e6
        if self.convention == "pbhstat":
            k_eq = 0.01 * (self.omega_m / 0.31)
            return np.sqrt(m_h / self.m_eq) / k_eq
        ln_k, ln_m, _ = self._table
        return np.exp(-np.interp(np.log(m_h), ln_m[::-1], ln_k[::-1]))

    def mh_of_r(self, r):
        r = np.asarray(r, dtype=float)
        if self.convention == "ferrante":
            m17 = 17.0 * (self.g_star / 10.75) ** (-1.0 / 6.0)
            return m17 * (1e6 * r) ** 2
        if self.convention == "pbhstat":
            k_eq = 0.01 * (self.omega_m / 0.31)
            return self.m_eq * (k_eq * r) ** 2
        ln_k, ln_m, _ = self._table
        return np.exp(np.interp(-np.log(r), ln_k, ln_m))

    def log_redshift_factor(self, m_h):
        """ln of the factor that converts beta(M_H) into its f_PBH weight per unit ln M_H."""
        m_h = np.asarray(m_h, dtype=float)
        if self.convention in ("ferrante", "pbhstat"):
            return 0.5 * np.log(self.m_eq / m_h) - np.log(self.omega_dm)
        _, ln_m, log_w = self._table
        return np.interp(np.log(m_h), ln_m[::-1], log_w[::-1])
