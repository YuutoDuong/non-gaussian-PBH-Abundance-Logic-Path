"""Phase 1 exit check: required amplitude A_* from ngpbh and from pbhstat, same conventions.

pbhstat conventions (2509.14326, code v0.0.6): peak-normalised lognormal, real-space top hat
with the k < 4.49/R cut and no transfer function, linear estimator g with g_c = 0.77, K = 4,
gamma = 0.36, high-peak number density with V = 4 pi R^3/3, R = (M_H/M_eq)^(1/2)/k_eq,
Omega_CDM = 0.26.

Two normalisations in the v0.0.6 code differ from the paper and from standard BBKS:
  * Press-Schechter has no factor 2 (the paper's eq. 7 writes beta = 2 int ...);
  * peaks integrates n(nu) over dg without dnu/dg = 1/sigma_0 and with an extra 1/2
    (paper eq. 17), i.e. beta_pbhstat = (sigma_0 / 2) beta_BBKS.
--mode pbhstat applies both to ngpbh so that the numerics can be compared;
--mode standard keeps ngpbh's conventions (factor 1 in PS, BBKS normalisation).

Usage: python scripts/compare_pbhstat.py [--fpbh 1e-3] [--kstar 1e6] [--mode pbhstat]
"""

import argparse
import csv
import time
import warnings
from pathlib import Path

import numpy as np
from scipy.optimize import brentq

from ngpbh.abundance import AbundanceModel
from ngpbh.cosmology import HorizonMass
from ngpbh.moments import MomentsConfig
from ngpbh.spectra import LogNormal
from ngpbh.stats import CollapseModel, PeaksConfig, QuadratureConfig

warnings.filterwarnings("ignore")
WIDTHS = (0.1, 0.3, 0.5, 1.0, 2.0)


class PbhstatPeaksNormalisation(AbundanceModel):
    """Peaks beta multiplied by sigma_0/2, reproducing pbhstat v0.0.6's PeaksTheoryModel."""

    def log_beta_grid(self, A):
        lb = super().log_beta_grid(A)
        sigma0 = np.array([np.sqrt(A * m.sigma_c2) for m in self.unit_moments])
        return lb + np.log(sigma0 / 2.0)


def ngpbh_model(width, kstar, statistic, mode):
    from pbhstat import constants

    cls = PbhstatPeaksNormalisation if (mode == "pbhstat" and statistic == "peaks") else AbundanceModel
    return cls(
        LogNormal(1.0, kstar, width, normalization="peak"),
        CollapseModel(C_th=0.55, K=4.0, gamma=0.36, estimator="linear", C_l_th=0.77),
        statistic=statistic,
        moments_cfg=MomentsConfig(window="tophat", transfer=False, x_cut=4.49),
        horizon=HorizonMass("pbhstat", omega_m=constants.Omega_m, omega_dm=constants.OmegaCDM),
        peaks=PeaksConfig(density="high_peak", volume="tophat"),
        quad=QuadratureConfig(ps_factor=1.0),
        n_mh=120,
        rel_sigma=1e-3,
    )


def pbhstat_fpbh(A, width, kstar, statistic):
    from pbhstat.mass_function import MassFunction
    from pbhstat.mass_variance import MassVariance
    from pbhstat.power_spectrum import PowerSpectrum

    stat = "press" if statistic == "threshold" else "peaks"
    ps = PowerSpectrum(shape="lognormal", amplitude=A, k_star=kstar, sigma_ln=width)
    # pbhstat returns NaN when k_values extend far beyond the peak; |ln k/k*| < 2 + 4 Delta is
    # stable (results change by < 1e-4 when this range is widened by 50%)
    half = 2.0 + 4.0 * width
    k = kstar * np.exp(np.linspace(-half, half, 2000))
    mv = MassVariance("realtophat", ps, statistics=stat, cutoff=True)
    mf = MassFunction(mass_variance=mv, statistics=stat, K=4, g_c=0.77, gamma=0.36)
    m, f = mf.evaluate(k, mpbh_vals=101)
    if not np.all(np.isfinite(f)):
        raise RuntimeError("pbhstat returned a non-finite mass function")
    return mf.fpbh(f, m)


def pbhstat_amplitude(width, kstar, statistic, f_target, guess):
    g = lambda lnA: np.log(pbhstat_fpbh(np.exp(lnA), width, kstar, statistic)) - np.log(f_target)
    lo, hi = np.log(guess / 1.3), np.log(guess * 1.3)
    return float(np.exp(brentq(g, lo, hi, xtol=1e-5)))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--fpbh", type=float, default=1e-3)
    ap.add_argument("--kstar", type=float, default=1e6)
    ap.add_argument("--mode", choices=("pbhstat", "standard"), default="pbhstat")
    ap.add_argument("--out", default=None)
    args = ap.parse_args()
    out = args.out or f"results/phase1_pbhstat_comparison_{args.mode}.csv"

    rows = []
    for width in WIDTHS:
        for statistic in ("threshold", "peaks"):
            t0 = time.time()
            model = ngpbh_model(width, args.kstar, statistic, args.mode)
            a_ours = model.solve_amplitude(args.fpbh, bracket=(1e-3, 3.0))
            t1 = time.time()
            a_ref = pbhstat_amplitude(width, args.kstar, statistic, args.fpbh, a_ours)
            t2 = time.time()
            diff = a_ours / a_ref - 1.0
            rows.append(dict(width=width, statistic=statistic, A_ngpbh=a_ours, A_pbhstat=a_ref, rel_diff=diff))
            print(
                f"Delta={width:3.1f} {statistic:9s} A_ngpbh={a_ours:.5e} A_pbhstat={a_ref:.5e} "
                f"diff={100 * diff:+.3f}%  ({t1 - t0:.1f}s / {t2 - t1:.1f}s)",
                flush=True,
            )
    Path(out).parent.mkdir(parents=True, exist_ok=True)
    with open(out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)
    worst = max(abs(r["rel_diff"]) for r in rows)
    print(f"max |A_ngpbh/A_pbhstat - 1| = {100 * worst:.3f}%  ->  {'PASS' if worst < 0.01 else 'FAIL'} (1% target)")


if __name__ == "__main__":
    main()
