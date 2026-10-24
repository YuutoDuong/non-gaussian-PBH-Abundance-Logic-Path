"""Primordial non-Gaussian corrections to the scalar-induced GW spectrum (radiation era).

For zeta = zeta_g + F (zeta_g^2 - <zeta_g^2>) the four-point function of zeta terminates at F^4 and
the spectrum splits into seven terms (Adshead, Lozanov & Weiner 2105.01659, eqs. 2.33-2.39):

    Omega = A^2 T2 + F^2 A^3 (hybrid + C + Z) + F^4 A^4 (reducible + planar + nonplanar)

for Delta^2_g = A s(k/k*), each term a function of x = k/k* only. All of them follow from

    Omega_X(k) = (k^3 / 12 pi^2) c_X int prod_loops d^3q/(2 pi)^3  Re[z1^2 conj(z2)^2] J12 prod_i P(m_i),

with z = (q_x + i q_y)/k the transverse momentum of a GW vertex (k along z), J12 = (I_A1 I_A2 +
I_B1 I_B2)/2 the oscillation average of the product of the two late-time transfer functions,
I = (I_A sin x + I_B cos x)/x (Kohri & Terada 1804.08577 eq. 26; J11 is the kernel of sigw.py),
and P(m) = 2 pi^2 Delta^2_g(m)/m^3 for each "spectral" momentum m_i, a linear combination of k and
the loop momenta. c_X = 2 (Gaussian), 8 (hybrid), 16 (C, Z), 8 (reducible), 32 (planar), 16
(nonplanar), per power of F; the disconnected terms have one vertex, q1 = q2.

The disconnected terms are the Gaussian functional of the full spectrum P_g + P_2 with
Delta^2_2 = F^2 A^2 s2(x) (one loop), so hybrid = 2 Omega_bil[s, s2] and reducible = Omega_G[s2];
they are computed by quadrature. The connected terms (5 and 8 dimensions) are computed by
randomised quasi-Monte Carlo over the spectral momenta: L of the L+1 spectral momenta are drawn
from the spectrum (ln m) with isotropic directions, the last is fixed by momentum conservation, and
the strategies (which momentum is fixed) are combined with the balance heuristic, which keeps the
weights bounded where one momentum becomes soft. For C and planar each sample is also averaged over
rotations about k (see _g_rotated), which removes the heavy-tailed variance of nearly independent
vertex pairs for broad spectra. The same engine evaluates the disconnected terms as a cross-check.
Validation against quadrature and against Adshead et al. Fig. 4: scripts/validate_sigw_ng.py.

Templates beyond quadratic order, zeta = zeta_g + F zeta_g^2 + G zeta_g^3 + ...: at O(A^3) only F
and G enter, the cubic term through the power renormalisation P -> (1 + 6 G sigma^2) P, i.e.
Omega += 12 G sigma^2 A^2 T2 (the g_NL-type connected trispectrum vanishes under the azimuthal
integral; Perna et al. 2403.06962 eq. 3.11). O(A^4) needs G^2, F^2 G, F H terms that are not
computed here; for those templates the series is truncated at O(A^3).
"""

from dataclasses import dataclass, replace

import numpy as np

from . import thermal
from .sigw import SQRT3, _G_SERIES, grid_for, omega_rd
from .spectra import PowerSpectrum


# --------------------------------------------------------------------------------------------
# transfer functions per vertex


def kernel_split(t, s, t2m3=None):
    """(a, I_A, I_B) with sigw.kernel(t, s) = a^2 (I_A^2 + I_B^2)/2.

    a = (t^2-1)(1-s^2)/(t^2-s^2) is the projection (|z|^2/uv), I_A and I_B the coefficients of
    sin x and cos x of the late-time transfer function times x (same variables as sigw.kernel).
    """
    t, s = np.broadcast_arrays(np.asarray(t, dtype=float), np.asarray(s, dtype=float))
    t2, s2 = t * t, s * s
    t2m3 = t2 - 3.0 if t2m3 is None else np.broadcast_to(t2m3, t.shape)
    y = 3.0 - s2
    xp = t2m3 + y
    a = (t2 - 1.0) * (1.0 - s2) / xp
    return (a,) + _transfer(xp, y, t2m3)


def _transfer(xp, y, t2m3):
    """I_A, I_B from xp = t^2 - s^2 = 4uv, y = 3 - s^2 and t2m3 = t^2 - 3 = xp - y."""
    xm = xp - 2.0 * y  # t^2 + s^2 - 6
    e = xp / y
    small = e < 0.05
    with np.errstate(divide="ignore", invalid="ignore"):
        log_term = np.log(np.abs(t2m3) / y)
        g_e3 = np.where(small, 0.0, (e + log_term - 0.5 * e * log_term) / e**3)
        if np.any(small):
            g_e3 = np.where(small, np.polynomial.polynomial.polyval(np.where(small, e, 0.0), _G_SERIES), g_e3)
        i_a = 24.0 * xm * g_e3 / y**2
        i_b = np.where(t2m3 > 0.0, 12.0 * np.pi * xm**2 / xp**3, 0.0)
    return i_a, i_b


def vertex(q, k):
    """z^2 = ((q_x + i q_y)/k)^2, I_A and I_B for internal momenta q (..., 3) and k = k z-hat."""
    qx, qy, qz = q[..., 0], q[..., 1], q[..., 2]
    perp2 = qx * qx + qy * qy
    qn = np.sqrt(perp2 + qz * qz)
    rn = np.sqrt(perp2 + (k - qz) ** 2)  # |k - q|
    s = (k - 2.0 * qz) / (rn + qn)  # (|k-q| - |q|)/k without cancellation
    xp = 4.0 * qn * rn / k**2
    y = 3.0 - s * s
    i_a, i_b = _transfer(xp, y, xp - y)
    z = (qx + 1j * qy) / k
    return z * z, i_a, i_b


# --------------------------------------------------------------------------------------------
# disconnected terms: one-loop spectrum and the bilinear Gaussian functional


@dataclass(frozen=True)
class Tabulated(PowerSpectrum):
    """A * shape tabulated in ln k (log-log interpolation); k^3 below the table, zero above."""

    A: float
    ln_k: np.ndarray
    ln_shape: np.ndarray
    scale: float  # ln k scale of the narrowest feature (for grid_for)

    def shape(self, k):
        lk = np.log(np.asarray(k, dtype=float))
        out = np.exp(np.interp(lk, self.ln_k, self.ln_shape))
        ir = lk < self.ln_k[0]
        out = np.where(ir, np.exp(self.ln_shape[0] + 3.0 * (lk - self.ln_k[0])), out)
        return np.where(lk > self.ln_k[-1], 0.0, out)

    def support(self, rel=1e-16):
        ok = self.ln_shape >= self.ln_shape.max() + np.log(rel)
        return float(np.exp(self.ln_k[ok][0])), float(np.exp(self.ln_k[ok][-1]))

    @property
    def lnk_scale(self):
        return self.scale

    @property
    def k_peak(self):
        return float(np.exp(self.ln_k[np.argmax(self.ln_shape)]))


def _unit(spectrum):
    unit = spectrum.with_amplitude(1.0)
    return replace(unit, k_star=1.0) if hasattr(unit, "k_star") else unit


def one_loop_shape(spectrum, n_y=12000, n_x=1600):
    """Delta^2_2 / (F^2 A^2) as a Tabulated spectrum, in units of k* (A = 1, k* = 1 internally).

    P_2(k) = 2 F^2 int d^3p/(2pi)^3 P(p) P(|k-p|) gives
    s2(x) = x^2 int dln p s(p)/p [C(ln(x+p)) - C(ln|x-p|)],  C(y) = int_-inf^y s(e^y') e^-y' dy'.
    """
    unit = _unit(spectrum)
    lo, hi = unit.support(1e-18)
    y = np.linspace(np.log(lo), np.log(hi), n_y)
    p = np.exp(y)
    sp = unit.shape(p)
    f = sp / p
    cum = np.concatenate([[0.0], np.cumsum(0.5 * (f[1:] + f[:-1]) * np.diff(y))])

    def big_c(z):
        with np.errstate(divide="ignore"):
            lz = np.log(z)
        return np.interp(lz, y, cum, left=0.0, right=cum[-1])

    xs = np.exp(np.linspace(np.log(lo) - 2.0, np.log(2.0 * hi), n_x))
    s2 = np.empty_like(xs)
    for i, x in enumerate(xs):
        integrand = f * (big_c(x + p) - big_c(np.abs(x - p)))
        s2[i] = x * x * np.trapezoid(integrand, y)
    keep = s2 > 0.0
    return Tabulated(1.0, np.log(xs[keep]), np.log(s2[keep]), unit.lnk_scale)


def omega_rd_bilinear(spec1, spec2, k, grid, chunk=8):
    """Omega_bil[P1, P2](k): the Gaussian functional with P(ku)P(kv) -> (P1(ku)P2(kv) + P2(ku)P1(kv))/2."""
    k = np.atleast_1d(np.asarray(k, dtype=float))
    out = np.empty_like(k)
    for i in range(0, k.size, chunk):
        kk = k[i : i + chunk, None, None]
        ku, kv = kk * grid.u, kk * grid.v
        p = 0.5 * (spec1(ku) * spec2(kv) + spec2(ku) * spec1(kv))
        out[i : i + chunk] = np.sum(grid.weight * p, axis=(1, 2))
    return out


def disconnected_terms(spectrum, x, s2=None, grid_kw=None):
    """Per unit A and F: {"gaussian": T2, "hybrid": 2 Omega_bil[s, s2], "reducible": Omega_G[s2]}."""
    unit = _unit(spectrum)
    s2 = one_loop_shape(spectrum) if s2 is None else s2
    x = np.atleast_1d(np.asarray(x, dtype=float))
    kw = grid_kw or {}
    g1 = grid_for(unit, x.min(), **kw)
    g2 = grid_for(s2, x.min(), **kw)
    gmix = g1 if g1.weight.size >= g2.weight.size else g2
    return {"gaussian": omega_rd(unit, x, grid=g1),
            "hybrid": 2.0 * omega_rd_bilinear(unit, s2, x, gmix),
            "reducible": omega_rd(s2, x, grid=g2)}


# --------------------------------------------------------------------------------------------
# Monte Carlo over the spectral momenta


@dataclass(frozen=True)
class Term:
    """rows[i] = coefficients of spectral momentum m_i over (k, loop momenta); c per power of F."""

    c: float
    rows: tuple
    connected: bool  # two GW vertices (q1, q2) or one (q1 = q2 = q, the first loop momentum)
    order: int  # power of F


TERMS = {
    "gaussian": Term(2.0, ((0, 1), (1, -1)), False, 0),
    "hybrid": Term(8.0, ((1, -1, 0), (0, 0, 1), (0, 1, -1)), False, 2),
    "C": Term(16.0, ((0, 0, 1), (1, 0, -1), (0, 1, -1)), True, 2),
    "Z": Term(16.0, ((0, 1, 0), (0, 0, 1), (1, -1, -1)), True, 2),
    "reducible": Term(8.0, ((0, 0, 1, 0), (0, 1, -1, 0), (0, 0, 0, 1), (1, -1, 0, -1)), False, 4),
    "planar": Term(32.0, ((0, 0, 0, 1), (1, 0, 0, -1), (0, 1, 0, -1), (0, 0, 1, -1)), True, 4),
    "nonplanar": Term(16.0, ((0, 0, 0, 1), (0, 1, 0, -1), (0, 0, 1, -1), (1, -1, -1, 1)), True, 4),
}


def _strategies(rows):
    """Omitted momentum j and the inverse of the loop-coefficient matrix of the others, if invertible."""
    rows = np.asarray(rows, dtype=float)
    out = []
    for j in range(rows.shape[0]):
        a = np.delete(rows, j, axis=0)[:, 1:]
        if abs(np.linalg.det(a)) > 0.5:
            out.append((j, np.linalg.inv(a)))
    return out


class LnSampler:
    """ln m drawn from (1 - eps) s(m)/sbar + eps U[a, b] (defensive log-uniform component)."""

    def __init__(self, spectrum, eps=0.1, extend_to=None, n=40000):
        self.unit = _unit(spectrum)
        lo, hi = self.unit.support(1e-16)
        y = np.linspace(np.log(lo), np.log(hi), n)
        s = self.unit.shape(np.exp(y))
        cdf = np.concatenate([[0.0], np.cumsum(0.5 * (s[1:] + s[:-1]) * np.diff(y))])
        self.sbar = cdf[-1]
        self._y, self._cdf = y, cdf / cdf[-1]
        lo8, hi8 = self.unit.support(1e-8)
        self.a = np.log(lo8)
        self.b = np.log(max(hi8, extend_to or 0.0))
        self.eps = eps

    def ppf(self, u):
        e = self.eps
        if e <= 0.0:
            return np.interp(u, self._cdf, self._y)
        flat = u < e
        y_s = np.interp(np.where(flat, 0.5, (u - e) / (1.0 - e)), self._cdf, self._y)
        return np.where(flat, self.a + (self.b - self.a) * u / e, y_s)

    def pdf(self, y):
        """Density in ln m."""
        dens = (1.0 - self.eps) * self.unit.shape(np.exp(y)) / self.sbar
        if self.eps > 0.0:
            dens = dens + np.where((y >= self.a) & (y <= self.b), self.eps / (self.b - self.a), 0.0)
        return dens


def _sobol(dim, n_log2, seed):
    """2^n_log2 scrambled Sobol points; seed is a sequence of ints (one stream per term/strategy)."""
    from scipy.stats import qmc

    rng = np.random.default_rng(np.random.SeedSequence(list(seed)))
    try:
        sob = qmc.Sobol(dim, scramble=True, rng=rng)
    except TypeError:  # SciPy < 1.15
        sob = qmc.Sobol(dim, scramble=True, seed=rng)
    return sob.random_base2(n_log2)


def _g(z1, a1, b1, z2, a2, b2):
    return 0.5 * np.real(z1 * np.conj(z2)) * (a1 * a2 + b1 * b2)


def _rotz(v, alpha):
    c, s = np.cos(alpha), np.sin(alpha)
    return np.stack([c * v[..., 0] - s * v[..., 1], s * v[..., 0] + c * v[..., 1], v[..., 2]], axis=-1)


def _g_rotated(name, loops, x, n_phi):
    """g averaged over rigid rotations about k of spectral momenta, n_phi angles each.

    C: q2 (and with it k - q2) with p = q1 - q2 fixed. Planar: p1 = q1 - q3 and p2 = q2 - q3
    independently, q3 fixed; g is bilinear in the two vertices, so the n_phi^2 average factorises.
    The rotations preserve the measure and every spectral magnitude, so the average is an unbiased
    estimator with the same weights; it cancels the cos 2(phi1 - phi2) oscillation of pairs of
    vertices that are nearly independent (one pair much harder or softer than k), which otherwise
    dominates the variance for broad spectra.
    """
    alphas = 2.0 * np.pi * np.arange(n_phi) / n_phi
    if name == "C":
        z2, a2, b2 = vertex(loops[1], x)
        p = loops[0] - loops[1]
        g = 0.0
        for al in alphas:
            g = g + _g(*vertex(p + _rotz(loops[1], al), x), z2 * np.exp(2j * al), a2, b2)
        return g / n_phi
    amp = []
    for i in (0, 1):
        p = loops[i] - loops[2]
        sa = sb = 0.0
        for al in alphas:
            z, a, b = vertex(loops[2] + _rotz(p, al), x)
            sa, sb = sa + z * a, sb + z * b
        amp.append((sa / n_phi, sb / n_phi))
    (a1, b1), (a2, b2) = amp
    return 0.5 * np.real(a1 * np.conj(a2) + b1 * np.conj(b2))


def mc_term(spectrum, name, x, n_log2=15, n_rand=4, seed=12345, eps=0.1, sampler=None, n_phi=8):
    """Term `name` per unit A^(L+1) F^order at x = k/k* by randomised QMC; returns (value, std. error).

    T_X(x) = (x^3 c_X / 6) sum_j E_j[ g prod_i s(m_i) / sum_l m_l^3 prod_(i != l) nu(ln m_i) ],
    with g = Re[z1^2 conj(z2)^2] J12, j and l over the valid strategies and nu the sampling density.
    For C and planar, g is averaged over n_phi rotations about k (see _g_rotated).
    """
    term = TERMS[name]
    rows = np.asarray(term.rows, dtype=float)
    n_mom, n_loop = rows.shape[0], rows.shape[1] - 1
    strategies = _strategies(rows)
    valid = [j for j, _ in strategies]
    smp = sampler if sampler is not None else LnSampler(spectrum, eps=eps, extend_to=3.0 * x)
    kvec = np.array([0.0, 0.0, x])
    est = np.zeros(n_rand)
    for r in range(n_rand):
        for j, ainv in strategies:
            u = _sobol(3 * n_loop, n_log2, [seed, r, j, n_mom])
            sampled = [i for i in range(n_mom) if i != j]
            mom = np.empty((n_mom, u.shape[0], 3))
            for c, i in enumerate(sampled):
                mag = np.exp(smp.ppf(u[:, 3 * c]))
                cos = 2.0 * u[:, 3 * c + 1] - 1.0
                sin = np.sqrt(np.maximum(1.0 - cos * cos, 0.0))
                phi = 2.0 * np.pi * u[:, 3 * c + 2]
                mom[i] = mag[:, None] * np.column_stack([sin * np.cos(phi), sin * np.sin(phi), cos])
            rhs = mom[sampled] - rows[sampled, 0][:, None, None] * kvec
            loops = np.einsum("ab,bnc->anc", ainv, rhs)
            mom[j] = rows[j, 0] * kvec + np.einsum("b,bnc->nc", rows[j, 1:], loops)
            if term.connected and name in ("C", "planar") and n_phi > 1:
                g = _g_rotated(name, loops, x, n_phi)
            elif term.connected:
                g = _g(*vertex(loops[0], x), *vertex(loops[1], x))
            else:
                z1, a1, b1 = vertex(loops[0], x)
                g = 0.5 * np.abs(z1) ** 2 * (a1 * a1 + b1 * b1)
            mags = np.linalg.norm(mom, axis=-1)
            with np.errstate(divide="ignore", over="ignore", under="ignore", invalid="ignore"):
                lnm = np.log(mags)
                num = g * np.prod(smp.unit.shape(mags), axis=0)
                nu = smp.pdf(lnm)
                den = np.zeros_like(num)
                for om in valid:
                    den += mags[om] ** 3 * np.prod(np.delete(nu, om, axis=0), axis=0)
                w = np.where(den > 0.0, num / den, 0.0)
            est[r] += np.mean(np.nan_to_num(w, nan=0.0, posinf=0.0, neginf=0.0))
    est *= x**3 * term.c / 6.0
    return float(est.mean()), float(est.std(ddof=1) / np.sqrt(n_rand)) if n_rand > 1 else np.nan


# --------------------------------------------------------------------------------------------
# local coefficients and the spectrum template


def local_coefficients(ng):
    """(F, G) of zeta = zeta_g + F zeta_g^2 + G zeta_g^3 + ... for a nongaussian.LocalNG template.

    A linear coefficient f1 != 1 is absorbed into zeta_g (F -> F/f1^2, G -> G/f1^3, A -> f1^2 A);
    it is 1 for every template in ngpbh.
    """
    if ng is None or ng.is_gaussian:
        return 0.0, 0.0
    z = np.array([0.0])
    f1 = float(ng.dF(z)[0])
    if abs(f1 - 1.0) > 1e-8:
        raise ValueError("templates with F'(0) != 1 need the amplitude rescaled; not supported")
    return float(ng.d2F(z)[0]) / 2.0, float(ng.d3F(z)[0]) / 6.0


@dataclass
class NGSIGWTemplate:
    """Omega_GW,RD(k) = A^2 T2 + A^3 (F^2 T3 + 12 G sbar T2) + F^4 A^4 T4, tabulated in x = k/k*.

    t2 is the fine Gaussian template (sigw.SIGWTemplate); T3 = hybrid + C + Z and T4 = reducible +
    planar + nonplanar are given on a coarser grid x_ng and interpolated as ratios to T2 (log x).
    """

    t2: object  # sigw.SIGWTemplate
    x_ng: np.ndarray
    r3: np.ndarray  # T3/T2 on x_ng
    r4: np.ndarray  # T4/T2 on x_ng
    sbar: float  # sigma^2 / A

    @classmethod
    def from_terms(cls, t2, terms):
        """From one shape's rows of results/sigw_ng_terms.csv (see load_terms)."""
        t3 = terms["hybrid"] + terms["C"] + terms["Z"]
        t4 = terms["reducible"] + terms["planar"] + terms["nonplanar"]
        return cls(t2, terms["x"], t3 / terms["T2"], t4 / terms["T2"], float(terms["sbar"][0]))

    def terms(self, x):
        x = np.asarray(x, dtype=float)
        t2 = self.t2(x)
        lx = np.log(np.clip(x, self.x_ng[0], self.x_ng[-1]))
        lg = np.log(self.x_ng)
        return t2, t2 * np.interp(lx, lg, self.r3), t2 * np.interp(lx, lg, self.r4)

    def omega(self, x, A, F=0.0, G=0.0, order=4):
        """Omega_GW,RD; A may broadcast against x (e.g. A[:, None], x[None, :])."""
        t2, t3, t4 = self.terms(x)
        out = A**2 * t2
        if order >= 3:
            out = out + A**3 * (F**2 * t3 + 12.0 * G * self.sbar * t2)
        if order >= 4:
            out = out + F**4 * A**4 * t4
        return out

    def h2_omega(self, f, A, k_star, F=0.0, G=0.0, order=4):
        x = np.asarray(f, dtype=float) / thermal.HZ_PER_INV_MPC / k_star
        return thermal.H2_OMEGA_R * thermal.gw_dilution(f) * self.omega(x, A, F, G, order)


def load_terms(path):
    """{shape: {column: array}} from a CSV written by scripts/sigw_ng_terms.py."""
    import csv

    out = {}
    with open(path) as fh:
        for r in csv.DictReader(fh):
            d = out.setdefault(r["shape"], {})
            for key, val in r.items():
                if key != "shape":
                    d.setdefault(key, []).append(float(val))
    return {s: {k: np.array(v) for k, v in d.items()} for s, d in out.items()}


__all__ = ["kernel_split", "vertex", "one_loop_shape", "omega_rd_bilinear", "disconnected_terms",
           "TERMS", "LnSampler", "mc_term", "local_coefficients", "NGSIGWTemplate", "load_terms", "SQRT3"]
