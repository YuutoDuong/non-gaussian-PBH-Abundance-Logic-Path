"""Phase 5: how strongly PBH overproduction disfavours a SIGW explanation of NANOGrav 15, per statistic.

Stage "lines": the overproduction line A_*(k*) (f_PBH = 1 and 1e-3) for every spectral shape, NG
template and statistic (threshold, peaks 3a, peaks 3b), on a grid in k*. All statistics share the
threshold, window, transfer function, k < 4.49/r cut, K = 3.3 and gamma = 0.36. Variants:
  base    - Musco threshold without the transfer function and the f_PBH normalisation of
            2306.17149 (eqs. 4, 14, 15): their setup, which our bands reproduce (phase5_validation);
  T       - Musco threshold with the transfer function (the Phase 1-4 default);
  entropy - base, with f_PBH from entropy conservation (12.7x smaller f_PBH for the same beta).
Stage "tension": the NG15 posterior in (log10 k*, log10 A) for each shape (uniform priors on
[5, 10] x [-3, 0], 14-bin HD free spectrum) and, for each line, the HPD credible level at which
the allowed region A < A_*(k*) is first reached (2306.17149's "excluded at n sigma"), the
posterior mass in it, and the profile-likelihood Delta chi^2.

Usage: python scripts/pta_tension.py [lines|tension|all] [--workers 4]
Writes results/phase5_lines.csv, results/phase5_tension.csv, results/phase5_summary.md.
"""

import argparse
import csv
import itertools
import time
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np
from scipy.interpolate import CubicSpline, RegularGridInterpolator

from ngpbh.abundance import AbundanceModel
from ngpbh.cosmology import HorizonMass
from ngpbh.moments import MomentsConfig
from ngpbh.nongaussian import Curvaton, Logarithmic, Quadratic
from ngpbh.pta import FreeSpectrum, credible_level, log_posterior_sigw, posterior_mass, sigma_1d
from ngpbh.sigw import SIGWTemplate
from ngpbh.spectra import BrokenPowerLaw, LogNormal
from ngpbh.stats import CollapseModel, PeaksConfig
from ngpbh.threshold import MuscoThreshold

K_REF = 1e7
SHAPES = {
    "BPL(4,3,1)": BrokenPowerLaw(1.0, K_REF, 4.0, 3.0, 1.0),
    "LN D=0.5": LogNormal(1.0, K_REF, 0.5),
    "LN D=1": LogNormal(1.0, K_REF, 1.0),
    "LN D=2": LogNormal(1.0, K_REF, 2.0),
}
TEMPLATES = {
    "gaussian": None,
    "fNL=-2": Quadratic(-2.0),
    "fNL=-1": Quadratic(-1.0),
    "fNL=-0.5": Quadratic(-0.5),
    "fNL=+0.5": Quadratic(0.5),
    "fNL=+1": Quadratic(1.0),
    "fNL=+2": Quadratic(2.0),
    "USR beta=3": Logarithmic.from_beta(3.0),  # 2306.17149's quasi-inflection-point case (f_NL = 5/4)
    "USR mu=5/6": Logarithmic(5.0 / 6.0),
    "USR mu=5/2": Logarithmic(2.5),
    "curvaton r=0.9": Curvaton(0.9),
}
STATISTICS = {
    "threshold": ("threshold", PeaksConfig()),
    "peaks_3a": ("peaks", PeaksConfig(constraint="C_g")),
    "peaks_3b": ("peaks", PeaksConfig(constraint="C_l")),
}
VARIANTS = {
    "base": dict(transfer=False, convention="fivt23"),
    "T": dict(transfer=True, convention="fivt23"),
    "entropy": dict(transfer=False, convention="entropy"),
}
ENTROPY_TEMPLATES = ("gaussian", "fNL=-2", "fNL=+2", "USR beta=3")
LOG10_KSTAR_LINE = np.arange(5.0, 10.001, 0.25)
F_TARGETS = (1.0, 1e-3)
LOG10_KSTAR = np.linspace(5.0, 10.0, 501)
LOG10_A = np.linspace(-3.0, 0.0, 601)


def line_job(job):
    shape, template, statistic, variant = job
    cfg = VARIANTS[variant]
    spec = SHAPES[shape]
    th = MuscoThreshold(transfer=cfg["transfer"]).evaluate(spec)
    stat, peaks = STATISTICS[statistic]
    model = AbundanceModel(
        spec,
        CollapseModel(C_th=th.C_th, K=3.3, gamma=0.36),
        ng=TEMPLATES[template],
        statistic=stat,
        moments_cfg=MomentsConfig(transfer=True, x_cut=4.49),
        horizon=HorizonMass(convention=cfg["convention"]),
        peaks=peaks,
        n_mh=50,
    )
    t0 = time.time()
    try:
        a = model.amplitude_curve(10**LOG10_KSTAR_LINE / K_REF, f_targets=F_TARGETS, bracket=(1e-4, 2.0))
    except RuntimeError as err:
        print(f"  {job}: {err}", flush=True)
        a = np.full((len(F_TARGETS), LOG10_KSTAR_LINE.size), np.nan)
    rows = []
    for i, f in enumerate(F_TARGETS):
        for lk, val in zip(LOG10_KSTAR_LINE, a[i]):
            rows.append(dict(shape=shape, template=template, statistic=statistic, variant=variant, C_th=round(th.C_th, 5),
                             f_pbh=f, log10_kstar=round(lk, 3), log10_A_star=np.log10(val)))
    return rows, time.time() - t0


def stage_lines(workers, out="results/phase5_lines.csv"):
    jobs = [(s, t, st, "base") for s, t, st in itertools.product(SHAPES, TEMPLATES, STATISTICS)]
    jobs += [(s, t, st, "T") for s, t, st in itertools.product(SHAPES, TEMPLATES, STATISTICS)]
    jobs += [(s, t, st, "entropy") for s, t, st in itertools.product(SHAPES, ENTROPY_TEMPLATES, STATISTICS)]
    jobs.sort(key=lambda j: j[2] != "peaks_3b")  # slowest first
    rows = []
    with ProcessPoolExecutor(workers) as pool:
        for job, (r, dt) in zip(jobs, pool.map(line_job, jobs)):
            rows += r
            mid = [x["log10_A_star"] for x in r if x["f_pbh"] == 1.0 and x["log10_kstar"] == 8.0][0]
            print(f"{job[0]:10s} {job[1]:15s} {job[2]:9s} {job[3]:7s} log10 A*(k*=1e8) = {mid:.4f} ({dt:.0f}s)", flush=True)
    Path(out).parent.mkdir(parents=True, exist_ok=True)
    with open(out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)


def load_lines(path="results/phase5_lines.csv"):
    lines = {}
    with open(path) as fh:
        for r in csv.DictReader(fh):
            key = (r["shape"], r["template"], r["statistic"], r["variant"], float(r["f_pbh"]))
            lines.setdefault(key, []).append((float(r["log10_kstar"]), float(r["log10_A_star"])))
    return {k: np.array(sorted(v)) for k, v in lines.items()}


def shape_posteriors():
    fs = FreeSpectrum(n_bins=14)
    return {name: log_posterior_sigw(SIGWTemplate(sp.with_amplitude(1.0)), fs, LOG10_KSTAR, LOG10_A)
            for name, sp in SHAPES.items()}


def tension_metrics(logl, line):
    """Credible level, posterior mass and profile Delta chi^2 of the allowed region A < A_*(k*)."""
    ok = np.isfinite(line[:, 1])
    a_star = CubicSpline(line[ok, 0], line[ok, 1])(LOG10_KSTAR)
    allowed = LOG10_A[None, :] < a_star[:, None]
    level = credible_level(logl, allowed)
    mass = posterior_mass(logl, allowed)
    dchi2 = 2.0 * (logl.max() - logl[allowed].max()) if allowed.any() else np.inf
    return level, mass, dchi2


def stage_tension(out="results/phase5_tension.csv"):
    lines = load_lines()
    posts = shape_posteriors()
    np.savez_compressed("results/phase5_posteriors.npz", log10_kstar=LOG10_KSTAR, log10_A=LOG10_A,
                        **{name: lp for name, lp in posts.items()})
    rows = []
    for (shape, template, statistic, variant, f), line in sorted(lines.items()):
        level, mass, dchi2 = tension_metrics(posts[shape], line)
        rows.append(dict(shape=shape, template=template, statistic=statistic, variant=variant, f_pbh=f,
                         log10_A_star_k8=float(np.interp(8.0, line[:, 0], line[:, 1])),
                         credible_level=level, sigma=sigma_1d(level), p_allowed=mass, dchi2=dchi2))
    with open(out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)
    write_summary(rows)


def write_summary(rows, out="results/phase5_summary.md"):
    idx = {(r["shape"], r["template"], r["statistic"], r["variant"], r["f_pbh"]): r for r in rows}
    md = ["# Phase 5: NG15 tension per statistic", "",
          "sigma = Gaussian-equivalent (two-sided 1D) of the HPD credible level at which the region "
          "A < A_*(k*) is first reached; 2306.17149 quote 'excluded at over 2 sigma' for their Gaussian band.",
          "Posterior: NG15 HD free spectrum (14 bins), uniform in log10 k* in [5,10] and log10 A in [-3,0], same shape "
          "as the line.", ""]
    for variant, f in [("base", 1.0), ("base", 1e-3), ("T", 1.0), ("entropy", 1.0)]:
        templates = ENTROPY_TEMPLATES if variant == "entropy" else TEMPLATES
        md += [f"## variant {variant}, f_PBH = {f:g}", "",
               "| shape | template | threshold | peaks 3a | peaks 3b | 3b - threshold | NG shift (3b) |",
               "|---|---|---|---|---|---|---|"]
        for shape in SHAPES:
            g = idx.get((shape, "gaussian", "peaks_3b", variant, f))
            for t in templates:
                s = [idx.get((shape, t, st, variant, f)) for st in STATISTICS]
                if None in s:
                    continue
                sig = [x["sigma"] for x in s]
                ng_shift = sig[2] - g["sigma"] if g else np.nan
                md.append(f"| {shape} | {t} | {sig[0]:.2f} | {sig[1]:.2f} | {sig[2]:.2f} | {sig[2] - sig[0]:+.2f} | "
                          f"{ng_shift:+.2f} |")
        md.append("")
    Path(out).write_text("\n".join(md))
    print("\n".join(md))


def stage_compare(out="results/phase5_compare.md"):
    """2306.17149's own lines (digitised, with QCD) against ours, on our posteriors.

    Their Fig. 3 overlays the BPL-marginalised NG15 posterior on bands for BPL(4,3,1); the
    marginalised posterior comes from scripts/pta_bpl_marginal.py.
    """
    dig = defaultdict(list)
    with open("results/digitized_2306.17149_fig3.csv") as fh:
        for r in csv.DictReader(fh):
            dig[(r["curve"], r["edge"])].append((float(r["log10_kstar"]), float(r["log10_A"])))
    lines = load_lines()
    posts = dict(np.load("results/phase5_posteriors.npz"))
    cases = {"BPL(4,3,1), fixed shape": posts["BPL(4,3,1)"]}
    marg = Path("results/phase5_posterior_bpl_marginal.npz")
    if marg.exists():
        cases["BPL marginalised (their Fig. 1/3)"] = np.load(marg)["logpost"]
    md = ["# Phase 5: 2306.17149's lines vs ours", "",
          "sigma of the region below each line (HPD level, 1D-equivalent). Their lines include QCD effects "
          "on the formation parameters; ours do not. Their lines extend over 6 < log10 k* < 9 only; beyond, "
          "they are continued flat.", "",
          "| posterior | template | f_PBH | their line | ngpbh threshold (base) | ngpbh peaks 3b (base) |",
          "|---|---|---|---|---|---|"]
    for pname, lp in cases.items():
        for t in ("gaussian", "fNL=-2", "USR beta=3", "curvaton r=0.9"):
            for f, edge in ((1.0, "f=1"), (1e-3, "f=1e-3")):
                d = np.array(sorted(dig[(t, edge)]))
                _, keep = np.unique(d[:, 0], return_index=True)
                d = d[keep]
                theirs = np.column_stack([LOG10_KSTAR, np.interp(LOG10_KSTAR, d[:, 0], d[:, 1])])
                sig = [sigma_1d(tension_metrics(lp, ln)[0]) for ln in
                       (theirs, lines[("BPL(4,3,1)", t, "threshold", "base", f)],
                        lines[("BPL(4,3,1)", t, "peaks_3b", "base", f)])]
                md.append(f"| {pname} | {t} | {f:g} | {sig[0]:.2f} | {sig[1]:.2f} | {sig[2]:.2f} |")
    if marg.exists():
        md += ["", "## Their NG15 contours inside our BPL-marginalised posterior", "",
               "HPD level of our posterior along their digitised contours (median, 10-90%), excluding the "
               "parts on the A = 1 prior edge and the plot edge at log10 k* = 9:", "",
               "| their contour | our HPD level | 1D convention | 2D convention |", "|---|---|---|---|"]
        lp = cases["BPL marginalised (their Fig. 1/3)"]
        p = np.exp(lp - lp.max())
        p /= p.sum()
        order = np.argsort(p.ravel())[::-1]
        lev = np.empty(p.size)
        lev[order] = np.cumsum(p.ravel()[order])
        interp = RegularGridInterpolator((LOG10_KSTAR, LOG10_A), lev.reshape(p.shape))
        for edge, c1, c2 in (("1sigma", 0.683, 0.393), ("2sigma", 0.954, 0.865), ("3sigma", 0.997, 0.989)):
            pts = np.array(dig[("NG15", edge)])
            pts = pts[(pts[:, 1] < -0.05) & (pts[:, 0] < 8.95)]
            v = interp(pts)
            md.append(f"| {edge} | {np.median(v):.3f} ({np.percentile(v, 10):.3f}-{np.percentile(v, 90):.3f}) | {c1} | {c2} |")
        q = lambda w, x: np.interp([0.0228, 0.5, 0.9772], np.cumsum(w) / w.sum(), x)
        k2, a2 = q(p.sum(1), LOG10_KSTAR), q(p.sum(0), LOG10_A)
        md += ["", f"1D marginals (median, 2-sigma): log10 k* = {k2[1]:.2f} (+{k2[2] - k2[1]:.2f}/-{k2[1] - k2[0]:.2f}), "
               f"log10 A = {a2[1]:.2f} (+{a2[2] - a2[1]:.2f}/-{a2[1] - a2[0]:.2f}); 2306.17149 Fig. 1: "
               "7.92 (+0.89/-0.62) and -0.86 (+0.79/-0.65)."]
    Path(out).write_text("\n".join(md) + "\n")
    print("\n".join(md))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("stage", nargs="?", default="all", choices=["lines", "tension", "compare", "all"])
    ap.add_argument("--workers", type=int, default=4)
    args = ap.parse_args()
    if args.stage in ("lines", "all"):
        stage_lines(args.workers)
    if args.stage in ("tension", "all"):
        stage_tension()
    if args.stage in ("compare", "all"):
        stage_compare()


if __name__ == "__main__":
    main()
