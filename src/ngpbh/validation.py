"""Independent checks of the Phase 3b peak density (docs/phase3b_derivation.tex, Sec. 9).

kac_rice_mc evaluates the Kac-Rice integrand for maxima of C_l = F'(zeta_r) C_g at one node
(C_g, zeta_r) = (c, z) by Monte Carlo, without using the reduction to Y = C_g + lam zeta_r:

* the 14 scalars-and-Hessian variables are drawn from their joint Gaussian conditional on
  (c, z), built component by component from the isotropic correlators;
* grad zeta_r is drawn from a proposal and weighted by the joint density of
  (grad C_g, grad zeta_r) at grad C_g = -lam grad zeta_r (the delta function of grad C_l);
* the Hessian of C_l is evaluated exactly, rank-one term included, and its determinant and
  definiteness are computed numerically.

It returns n(c, z) / P_G(c, z) in units of r^-3, with and without the rank-one term, to be
compared with NGPeakGeometry.log_weight_ng minus log V.
"""

import numpy as np

# order of the six independent components of a symmetric 3x3 matrix
_PAIRS = ((0, 0), (1, 1), (2, 2), (0, 1), (0, 2), (1, 2))


def _hess_cov(i, j, k, l):
    """(delta_ij delta_kl + delta_ik delta_jl + delta_il delta_jk) / 15."""
    d = lambda a, b: 1.0 if a == b else 0.0
    return (d(i, j) * d(k, l) + d(i, k) * d(j, l) + d(i, l) * d(j, k)) / 15.0


def _blocks(mom):
    """Second moments in units r = 1: sig0, sig1, sig2 as 2x2 matrices over (C_g, zeta_r)."""
    sig0 = np.array([[mom.sigma_c2, mom.sigma_cr2], [mom.sigma_cr2, mom.sigma_r2]])
    sig1 = np.array([[mom.s1_2, mom.scr1_2], [mom.scr1_2, mom.s1r_2]])
    sig2 = np.array([[mom.s2_2, mom.s2cr_2], [mom.s2cr_2, mom.s2r_2]])
    return sig0, sig1, sig2


def scalar_hessian_covariance(mom):
    """Covariance of (C_g, zeta_r, H_c[6], H_r[6]) from the isotropic correlators."""
    sig0, sig1, sig2 = _blocks(mom)
    n = 14
    cov = np.zeros((n, n))
    cov[:2, :2] = sig0
    for A in range(2):
        for B in range(2):
            for p, (i, j) in enumerate(_PAIRS):
                # <A d_i d_j B> = -delta_ij sig1_AB / 3
                cov[A, 2 + 6 * B + p] = cov[2 + 6 * B + p, A] = -(1.0 if i == j else 0.0) * sig1[A, B] / 3.0
                for q, (k, l) in enumerate(_PAIRS):
                    cov[2 + 6 * A + p, 2 + 6 * B + q] = sig2[A, B] * _hess_cov(i, j, k, l)
    return cov


def _sym(v):
    """(N, 6) components -> (N, 3, 3) symmetric matrices."""
    out = np.empty(v.shape[:-1] + (3, 3))
    for p, (i, j) in enumerate(_PAIRS):
        out[..., i, j] = v[..., p]
        out[..., j, i] = v[..., p]
    return out


def kac_rice_mc(mom, ng, c, z, n_samples=1_000_000, seed=0, proposal_var=None):
    """Monte Carlo n(c, z) / P_G(c, z) for maxima of C_l; returns (exact, without rank-one term,
    standard errors of both).  mom: CompactionMoments at the actual amplitude (r = 1 units)."""
    rng = np.random.default_rng(seed)
    d1, d2, d3 = float(ng.dF(z)), float(ng.d2F(z)), float(ng.d3F(z))
    lam = c * d2 / d1

    # Hessians conditional on the scalars
    cov = scalar_hessian_covariance(mom)
    s_idx, h_idx = slice(0, 2), slice(2, 14)
    css = cov[s_idx, s_idx]
    chs = cov[h_idx, s_idx]
    chh = cov[h_idx, h_idx]
    gain = chs @ np.linalg.inv(css)
    mean_h = gain @ np.array([c, z])
    cov_h = chh - gain @ chs.T
    L = np.linalg.cholesky(cov_h + 1e-14 * np.eye(12) * np.trace(cov_h))
    h = mean_h + rng.standard_normal((n_samples, 12)) @ L.T
    Hc, Hr = _sym(h[:, :6]), _sym(h[:, 6:])

    # gradients: per component, (g_c, g_r) Gaussian with covariance sig1 / 3, components independent
    _, sig1, _ = _blocks(mom)
    cg = sig1 / 3.0
    inv_cg = np.linalg.inv(cg)
    det_cg = np.linalg.det(cg)
    pv = cg[1, 1] if proposal_var is None else proposal_var
    g = rng.standard_normal((n_samples, 3)) * np.sqrt(pv)
    gc = -lam * g
    # log of the 6-dim joint density at (gc, g), minus log of the proposal density
    quad = inv_cg[0, 0] * gc**2 + 2.0 * inv_cg[0, 1] * gc * g + inv_cg[1, 1] * g**2
    log_p = -0.5 * quad.sum(axis=1) - 3.0 * np.log(2.0 * np.pi) - 1.5 * np.log(det_cg)
    log_q = -0.5 * (g**2).sum(axis=1) / pv - 1.5 * np.log(2.0 * np.pi * pv)
    w = np.exp(log_p - log_q)

    outer_cr = gc[:, :, None] * g[:, None, :]
    rank_one = d2 * (outer_cr + np.swapaxes(outer_cr, 1, 2)) + d3 * c * g[:, :, None] * g[:, None, :]
    base = d1 * Hc + d2 * c * Hr
    results = []
    for M in (base + rank_one, base):
        ev = np.linalg.eigvalsh(M)
        is_max = np.all(ev < 0.0, axis=1)
        val = np.where(is_max, np.abs(np.prod(ev, axis=1)), 0.0) * w / abs(d1) ** 3
        results.append((val.mean(), val.std() / np.sqrt(n_samples)))
    (exact, e_exact), (no_r1, e_no_r1) = results
    return exact, no_r1, e_exact, e_no_r1


def semi_analytic_weight(mom, ng, c, z, rank_one=True):
    """n(c, z) / P_G(c, z) as in ngpbh.stats (units r^-3, no volume factor): eq. (17) of the
    derivation, times the rank-one factor of Sec. 7 if rank_one."""
    from .stats import NGPeakGeometry, PeaksConfig

    geo = NGPeakGeometry.from_moments(mom, PeaksConfig(constraint="C_l"))
    d1, d2, d3 = (np.array(float(f(z))) for f in (ng.dF, ng.d2F, ng.d3F))
    c_l = c * d1
    nu, zz = np.array(c / mom.sigma_c), np.array(z / mom.sigma_r)
    # geo holds the moments at this amplitude, so A^(1/2) relative to them is 1
    log_w = geo.log_weight_ng(nu, zz, c_l, d1, d2, d3, 1.0, rank_one)
    return float(np.exp(log_w - geo.log_volume))
