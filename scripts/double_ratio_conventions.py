"""Does the double ratio R(NG)/R(Gaussian) cancel the normalisation conventions?

Physics setup of scripts/scan_amplitudes.py (Musco threshold with transfer function, k < 4.49/r,
K = 3.3, gamma = 0.36); only the Press-Schechter factor and the peak volume change:
  bbks   : PS factor 1, V = 4 pi r^3 / 3   (ngpbh default)
  ps2    : PS factor 2, V = 4 pi r^3 / 3
  gaussV : PS factor 1, V = (2 pi)^{3/2} r^3
"""

import warnings
from concurrent.futures import ProcessPoolExecutor

from ngpbh.abundance import AbundanceModel
from ngpbh.moments import MomentsConfig
from ngpbh.nongaussian import Logarithmic, Quadratic
from ngpbh.spectra import LogNormal
from ngpbh.stats import CollapseModel, PeaksConfig, QuadratureConfig
from ngpbh.threshold import MuscoThreshold

warnings.filterwarnings("ignore")
TEMPLATES = {"gaussian": None, "fNL=+1": Quadratic(1.0), "fNL=-2": Quadratic(-2.0), "USR mu=5/6": Logarithmic(5.0 / 6.0)}
CASES = ((0.5, ("fNL=+1", "fNL=-2")), (2.0, ("USR mu=5/6",)))
CONVENTIONS = {"bbks": (1.0, "tophat"), "ps2": (2.0, "tophat"), "gaussV": (1.0, "gaussian")}


def solve(job):
    width, label, statistic, ps_factor, volume = job
    spec = LogNormal(1.0, 1e7, width)
    th = MuscoThreshold(transfer=True).evaluate(spec)
    model = AbundanceModel(
        spec,
        CollapseModel(C_th=th.C_th, K=3.3, gamma=0.36),
        ng=TEMPLATES[label],
        statistic=statistic,
        moments_cfg=MomentsConfig(transfer=True, x_cut=4.49),
        peaks=PeaksConfig(volume=volume),
        quad=QuadratureConfig(ps_factor=ps_factor),
        n_mh=50,
    )
    return job, model.solve_amplitude(1.0, bracket=(1e-4, 2.0))


def main():
    jobs = {
        (w, t, s, ps, vol)
        for w, labels in CASES
        for t in ("gaussian",) + labels
        for ps, vol in CONVENTIONS.values()
        for s in ("threshold", "peaks")
    }
    with ProcessPoolExecutor(4) as pool:
        A = dict(pool.map(solve, sorted(jobs)))

    def R(w, t, conv):
        ps, vol = CONVENTIONS[conv]
        return A[(w, t, "peaks", ps, vol)] / A[(w, t, "threshold", ps, vol)]

    for w, labels in CASES:
        print(f"Delta={w} R_gauss: " + " ".join(f"{c}={R(w, 'gaussian', c):.4f}" for c in CONVENTIONS))
        for t in labels:
            print(f"  {t:11s} R/R_gauss: " + " ".join(f"{c}={R(w, t, c) / R(w, 'gaussian', c):.4f}" for c in CONVENTIONS))


if __name__ == "__main__":
    main()
