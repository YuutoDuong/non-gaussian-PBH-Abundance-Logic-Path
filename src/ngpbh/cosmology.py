"""Horizon-mass relations and the constants that convert beta(M_H) into f_PBH.

Two conventions are provided because the project compares against both:

* "ferrante": M_H = 17 Msun (g*/10.75)^(-1/6) (k_H / 1e6 Mpc^-1)^(-2) with k_H = 1/r,
  i.e. r_m = 2.4e-7 Mpc (g*/10.75)^(1/12) (M_H/Msun)^(1/2)  (2211.01728 eq. 39).
* "pbhstat": R = (M_H / M_eq)^(1/2) / k_eq with k_eq = 0.01 (Omega_m/0.31) Mpc^-1
  (2509.14326 eq. 4).

In both, f_PBH = (1/Omega_DM) int dln M_H (M_eq/M_H)^(1/2) beta(M_H).
"""

from dataclasses import dataclass

import numpy as np

M_EQ = 2.8e17  # horizon mass at matter-radiation equality [Msun]
OMEGA_DM = 0.12 / 0.674**2  # 2211.01728, 2306.17149


@dataclass(frozen=True)
class HorizonMass:
    convention: str = "ferrante"
    g_star: float = 106.75
    omega_m: float = 0.31
    m_eq: float = M_EQ
    omega_dm: float = OMEGA_DM

    def __post_init__(self):
        if self.convention not in ("ferrante", "pbhstat"):
            raise ValueError("convention must be 'ferrante' or 'pbhstat'")

    def r_of_mh(self, m_h):
        """Comoving scale [Mpc] that crosses the horizon when the horizon mass is m_h [Msun]."""
        m_h = np.asarray(m_h, dtype=float)
        if self.convention == "ferrante":
            m17 = 17.0 * (self.g_star / 10.75) ** (-1.0 / 6.0)
            return np.sqrt(m_h / m17) / 1e6
        k_eq = 0.01 * (self.omega_m / 0.31)
        return np.sqrt(m_h / self.m_eq) / k_eq

    def mh_of_r(self, r):
        r = np.asarray(r, dtype=float)
        if self.convention == "ferrante":
            m17 = 17.0 * (self.g_star / 10.75) ** (-1.0 / 6.0)
            return m17 * (1e6 * r) ** 2
        k_eq = 0.01 * (self.omega_m / 0.31)
        return self.m_eq * (k_eq * r) ** 2

    def log_redshift_factor(self, m_h):
        """ln[(M_eq / M_H)^(1/2) / Omega_DM]: converts beta(M_H) into its f_PBH weight."""
        return 0.5 * np.log(self.m_eq / np.asarray(m_h, dtype=float)) - np.log(self.omega_dm)
