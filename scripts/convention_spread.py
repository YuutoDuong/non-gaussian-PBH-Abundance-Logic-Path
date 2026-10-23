"""How much of R = A_*^peaks / A_*^threshold is convention (Gaussian, pbhstat physics setup).

Same spectra, window, estimator and constants as scripts/compare_pbhstat.py; only the
normalisations differ:
  pbhstat : PS factor 1, peaks beta x sigma_0/2       (pbhstat v0.0.6 code)
  bbks    : PS factor 1, BBKS peaks, V = 4 pi r^3 / 3  (ngpbh default)
  ps2     : PS factor 2, BBKS peaks, V = 4 pi r^3 / 3  (pbhstat paper eq. 7)
  gaussV  : PS factor 1, BBKS peaks, V = (2 pi)^{3/2} r^3 (Young 2201.13345)
"""

import csv
import warnings
from pathlib import Path

import numpy as np

from ngpbh.abundance import AbundanceModel
from ngpbh.cosmology import HorizonMass
from ngpbh.moments import MomentsConfig
from ngpbh.spectra import LogNormal
from ngpbh.stats import CollapseModel, PeaksConfig, QuadratureConfig

warnings.filterwarnings("ignore")
WIDTHS = (0.1, 0.3, 0.5, 1.0, 2.0)


class SigmaHalf(AbundanceModel):
    def log_beta_grid(self, A):
        lb = super().log_beta_grid(A)
        return lb + np.log(np.array([np.sqrt(A * m.sigma_c2) for m in self.unit_moments]) / 2.0)


def amplitude(width, statistic, ps_factor=1.0, volume="tophat", sigma_half=False, kstar=1e6, fpbh=1e-3):
    cls = SigmaHalf if sigma_half else AbundanceModel
    model = cls(
        LogNormal(1.0, kstar, width, normalization="peak"),
        CollapseModel(C_th=0.55, K=4.0, gamma=0.36, estimator="linear", C_l_th=0.77),
        statistic=statistic,
        moments_cfg=MomentsConfig(window="tophat", transfer=False, x_cut=4.49),
        horizon=HorizonMass("pbhstat", omega_m=0.315, omega_dm=0.26),
        peaks=PeaksConfig(density="high_peak", volume=volume),
        quad=QuadratureConfig(ps_factor=ps_factor),
        n_mh=120,
        rel_sigma=1e-3,
    )
    return model.solve_amplitude(fpbh, bracket=(1e-3, 3.0))


def main():
    rows = []
    for w in WIDTHS:
        ps1 = amplitude(w, "threshold")
        ps2 = amplitude(w, "threshold", ps_factor=2.0)
        pk = amplitude(w, "peaks")
        pk_half = amplitude(w, "peaks", sigma_half=True)
        pk_gv = amplitude(w, "peaks", volume="gaussian")
        row = dict(width=w, R_pbhstat=pk_half / ps1, R_bbks=pk / ps1, R_ps2=pk / ps2, R_gaussV=pk_gv / ps1)
        rows.append(row)
        print(" ".join(f"{k}={v:.4f}" for k, v in row.items()), flush=True)
    out = Path("results/phase1_convention_spread.csv")
    with open(out, "w", newline="") as fh:
        wr = csv.DictWriter(fh, fieldnames=list(rows[0]))
        wr.writeheader()
        wr.writerows(rows)


if __name__ == "__main__":
    main()
