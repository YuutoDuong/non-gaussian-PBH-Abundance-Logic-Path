"""Standard-Model thermal history: g*(T), g*s(T), and when a comoving scale re-enters the horizon.

g*(T) and g*s(T) are the Saikawa & Shirai table (2005.03544, N_eff = 3.045) as distributed with
PTArcade 1.1.5, so present-day GW spectra use the same degrees of freedom as the NANOGrav
new-physics analyses. Radiation domination throughout:

    H(T) = pi sqrt(g*/90) T^2 / M_P,      a/a_0 = (g*s0 / g*s)^(1/3) T_0 / T,
    k(T) = a H  (comoving, Mpc^-1),       f = k / (2 pi)  (present-day frequency),
    M_H(T) = 4 pi M_P^2 / H              (horizon mass).
"""

from importlib.resources import files

import numpy as np

M_P = 2.435323e18  # reduced Planck mass [GeV]
T_0 = 2.7255 * 8.617333262e-14  # CMB temperature today [GeV]
HBARC_GEV_M = 1.973269804e-16  # hbar c [GeV m]
MPC_M = 3.0856775814913673e22
GEV_INV_MPC = MPC_M / HBARC_GEV_M  # 1 GeV in Mpc^-1
HZ_PER_INV_MPC = 299792458.0 / MPC_M / (2.0 * np.pi)  # f [Hz] = HZ_PER_INV_MPC * k [Mpc^-1]
H100_HZ = 1e5 / MPC_M  # 100 km/s/Mpc in s^-1
M_SUN_GEV = 1.98841e30 * 299792458.0**2 / 1.602176634e-10
H2_OMEGA_R = 4.2e-5  # photons + neutrinos today (2306.17149 eq. 3)

_TAB = np.loadtxt(files("ngpbh").joinpath("data/g_star.dat"))
_LN_T = np.log(_TAB[:, 0])
G_S0 = float(_TAB[0, 2])  # table starts after e+e- annihilation (T = 10 keV)
G_RHO0 = float(_TAB[0, 3])


def g_rho(T):
    """Energy degrees of freedom at temperature T [GeV] (constant outside the table)."""
    return np.interp(np.log(T), _LN_T, _TAB[:, 3])


def g_s(T):
    """Entropy degrees of freedom at temperature T [GeV]."""
    return np.interp(np.log(T), _LN_T, _TAB[:, 2])


def hubble(T):
    """H [GeV] at temperature T [GeV]."""
    T = np.asarray(T, dtype=float)
    return np.pi * np.sqrt(g_rho(T) / 90.0) * T**2 / M_P


def k_of_T(T):
    """Comoving wavenumber [Mpc^-1] equal to aH at temperature T [GeV]."""
    T = np.asarray(T, dtype=float)
    return (G_S0 / g_s(T)) ** (1.0 / 3.0) * T_0 / T * hubble(T) * GEV_INV_MPC


_LN_T_FINE = np.linspace(np.log(1e-6), np.log(1e7), 4001)
_LN_K_FINE = np.log(k_of_T(np.exp(_LN_T_FINE)))


def T_of_k(k):
    """Temperature [GeV] at which the comoving scale k [Mpc^-1] re-enters the horizon."""
    return np.exp(np.interp(np.log(k), _LN_K_FINE, _LN_T_FINE))


def T_of_f(f):
    """Temperature [GeV] at horizon re-entry of the mode with present-day frequency f [Hz]."""
    return T_of_k(np.asarray(f, dtype=float) / HZ_PER_INV_MPC)


def horizon_mass(T):
    """Horizon mass [Msun] at temperature T [GeV]."""
    return 4.0 * np.pi * M_P**2 / hubble(T) / M_SUN_GEV


def gw_dilution(f):
    """(g*/g*0)(g*s/g*s0)^(-4/3) at the horizon re-entry of frequency f [Hz] (2306.17149 eq. 3)."""
    T = T_of_f(f)
    return g_rho(T) / G_RHO0 * (g_s(T) / G_S0) ** (-4.0 / 3.0)
