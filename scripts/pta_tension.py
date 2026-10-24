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

Stage "ngsigw": the same with the primordial-NG corrections to the SIGW spectrum (Adshead et al.
2105.01659; ngpbh.sigw_ng, terms from scripts/sigw_ng_terms.py): every NG template gets its own
posterior, exact to A^4 for the quadratic template and to A^3 (F and G) for the curvaton and USR
templates. The Gaussian-SIGW posterior of stage "tension" is what 2306.17149 used.

Usage: python scripts/pta_tension.py [lines|tension|compare|ngsigw|all] [--workers 4]
Writes results/phase5_lines.csv, results/phase5_tension.csv, results/phase5_summary.md and, for
stage ngsigw, results/phase5_ngsigw_tension.csv and results/phase5_ngsigw_summary.md.
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
from ngpbh.pta import (FreeSpectrum, credible_level, log_posterior_sigw, log_posterior_sigw_ng, posterior_mass,
                       sigma_1d, sigma_one_sided)
from ngpbh.sigw import SIGWTemplate
from ngpbh.sigw_ng import NGSIGWTemplate, load_terms, local_coefficients
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
SMOOTH = 0.05  # dex; regularises the KDE tails (ngpbh.pta.FreeSpectrum), bulk unchanged
SIGMA_RELIABLE = 4.0  # beyond this the free-spectrum KDEs carry no information (see docs/phase5_pta.md)
NG_TERMS = "results/sigw_ng_terms.csv"
NG_FIGURE_CASES = [("BPL(4,3,1)", "fNL=-2"), ("BPL(4,3,1)", "fNL=+2"), ("LN D=0.5", "fNL=-2"), ("LN D=0.5", "fNL=+2")]


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


def shape_posteriors(smooth=SMOOTH):
    fs = FreeSpectrum(n_bins=14, smooth=smooth)
    return {name: log_posterior_sigw(SIGWTemplate(sp.with_amplitude(1.0)), fs, LOG10_KSTAR, LOG10_A)
            for name, sp in SHAPES.items()}


def compact(logpost, decimals):
    """ln L relative to its maximum, clipped at -200 (beyond any quoted level), as float32."""
    return np.round(np.maximum(logpost - logpost.max(), -200.0), int(-np.log10(decimals))).astype(np.float32)


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
    exact = shape_posteriors(smooth=None)  # ceffyl's KDEs as distributed
    # for figures only: every second grid point, ln L relative to its maximum, float32
    np.savez_compressed("results/phase5_posteriors.npz", log10_kstar=LOG10_KSTAR[::2], log10_A=LOG10_A[::2],
                        **{name: compact(lp[::2, ::2], 1e-3) for name, lp in posts.items()})
    rows = []
    for (shape, template, statistic, variant, f), line in sorted(lines.items()):
        level, mass, dchi2 = tension_metrics(posts[shape], line)
        level_x, mass_x, _ = tension_metrics(exact[shape], line)
        rows.append(dict(shape=shape, template=template, statistic=statistic, variant=variant, f_pbh=f,
                         log10_A_star_k8=float(np.interp(8.0, line[:, 0], line[:, 1])),
                         credible_level=level, sigma=sigma_1d(level), p_allowed=mass,
                         sigma_p=sigma_one_sided(mass), dchi2=dchi2,
                         sigma_ceffyl=sigma_1d(level_x), sigma_p_ceffyl=sigma_one_sided(mass_x)))
    with open(out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)
    write_summary(rows)


def fmt_sigma(x):
    return f">{SIGMA_RELIABLE:.0f}" if x > SIGMA_RELIABLE else f"{x:.2f}"


def write_summary(rows, out="results/phase5_summary.md"):
    idx = {(r["shape"], r["template"], r["statistic"], r["variant"], r["f_pbh"]): r for r in rows}
    md = ["# Phase 5: NG15 tension per statistic", "",
          "Each cell: sigma_HPD / sigma_P. sigma_HPD is the two-sided 1D Gaussian equivalent of the HPD credible "
          "level at which the region A < A_*(k*) is first reached (2306.17149's 'excluded at n sigma'); sigma_P is "
          "the one-sided Gaussian equivalent of the posterior mass in that region. Posterior: NG15 HD free spectrum "
          f"(14 bins, KDE tails smoothed by {SMOOTH} dex), uniform in log10 k* in [5,10] and log10 A in [-3,0], same "
          f"shape as the line. Values above {SIGMA_RELIABLE:.0f} sigma lie where the free-spectrum KDEs have no samples.",
          ""]
    for variant, f in [("base", 1.0), ("base", 1e-3), ("T", 1.0), ("entropy", 1.0)]:
        templates = ENTROPY_TEMPLATES if variant == "entropy" else TEMPLATES
        md += [f"## variant {variant}, f_PBH = {f:g}", "",
               "| shape | template | log10 A* (k*=1e8): threshold, 3b | threshold | peaks 3a | peaks 3b | "
               "3b - threshold (P) | NG shift, 3b (P) |",
               "|---|---|---|---|---|---|---|---|"]
        for shape in SHAPES:
            g = idx.get((shape, "gaussian", "peaks_3b", variant, f))
            for t in templates:
                s = [idx.get((shape, t, st, variant, f)) for st in STATISTICS]
                if None in s:
                    continue
                cells = [f"{fmt_sigma(x['sigma'])} / {fmt_sigma(x['sigma_p'])}" for x in s]
                d_stat = s[2]["sigma_p"] - s[0]["sigma_p"]
                d_ng = s[2]["sigma_p"] - g["sigma_p"] if g else np.nan
                md.append(f"| {shape} | {t} | {s[0]['log10_A_star_k8']:.3f}, {s[2]['log10_A_star_k8']:.3f} | "
                          f"{cells[0]} | {cells[1]} | {cells[2]} | {d_stat:+.2f} | {d_ng:+.2f} |")
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
    fs = FreeSpectrum(n_bins=14, smooth=SMOOTH)
    bpl = log_posterior_sigw(SIGWTemplate(SHAPES["BPL(4,3,1)"]), fs, LOG10_KSTAR, LOG10_A)
    cases = {"BPL(4,3,1), fixed shape": bpl}
    marg = Path("results/phase5_posterior_bpl_marginal.npz")
    if marg.exists():
        cases["BPL marginalised (their Fig. 1/3)"] = np.load(marg)["logpost"]
    md = ["# Phase 5: 2306.17149's lines vs ours", "",
          "sigma_HPD / sigma_P of the region below each line (as in phase5_summary.md). Their lines include QCD "
          "effects on the formation parameters; ours do not. Their lines extend over 6 < log10 k* < 9 only; "
          "beyond, they are continued flat.", "",
          "| posterior | template | f_PBH | their line | ngpbh threshold (base) | ngpbh peaks 3b (base) |",
          "|---|---|---|---|---|---|"]
    for pname, lp in cases.items():
        for t in ("gaussian", "fNL=-2", "USR beta=3", "curvaton r=0.9"):
            for f, edge in ((1.0, "f=1"), (1e-3, "f=1e-3")):
                d = np.array(sorted(dig[(t, edge)]))
                _, keep = np.unique(d[:, 0], return_index=True)
                d = d[keep]
                theirs = np.column_stack([LOG10_KSTAR, np.interp(LOG10_KSTAR, d[:, 0], d[:, 1])])
                cells = []
                for ln in (theirs, lines[("BPL(4,3,1)", t, "threshold", "base", f)],
                           lines[("BPL(4,3,1)", t, "peaks_3b", "base", f)]):
                    level, mass, _ = tension_metrics(lp, ln)
                    cells.append(f"{fmt_sigma(sigma_1d(level))} / {fmt_sigma(sigma_one_sided(mass))}")
                md.append(f"| {pname} | {t} | {f:g} | {cells[0]} | {cells[1]} | {cells[2]} |")
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


def ng_templates(path=NG_TERMS):
    terms = load_terms(path)
    return {name: NGSIGWTemplate.from_terms(SIGWTemplate(sp.with_amplitude(1.0)), terms[name])
            for name, sp in SHAPES.items()}


def ng_settings(template):
    """(F, G, order, alternative order) of the NG-SIGW series for a Phase 5 template.

    The quadratic template is exact at A^4; its O(A^3) truncation is the alternative. The
    curvaton and USR maps are truncated at A^3 (F, G); adding the F^4 A^4 terms, without the
    G^2, F^2 G and F H ones, is the alternative.
    """
    ng = TEMPLATES[template]
    F, G = local_coefficients(ng)
    return (F, G, 4, 3) if isinstance(ng, Quadratic) else (F, G, 3, 4)


def stage_ngsigw(out="results/phase5_ngsigw_tension.csv"):
    lines = load_lines()
    with open("results/phase5_tension.csv") as fh:
        gauss = {(r["shape"], r["template"], r["statistic"], r["variant"], float(r["f_pbh"])): r for r in csv.DictReader(fh)}
    tmpl = ng_templates()
    fs = FreeSpectrum(n_bins=14, smooth=SMOOTH)
    fsx = FreeSpectrum(n_bins=14, smooth=None)
    rows, figs, med = [], {}, {}
    for shape in SHAPES:
        lp_g = log_posterior_sigw(tmpl[shape].t2, fs, LOG10_KSTAR, LOG10_A)
        med[(shape, "gaussian")] = posterior_medians(lp_g)
        for template in TEMPLATES:
            if template == "gaussian":
                continue
            F, G, order, alt = ng_settings(template)
            t0 = time.time()
            lp = log_posterior_sigw_ng(tmpl[shape], fs, LOG10_KSTAR, LOG10_A, F, G, order)
            lpx = log_posterior_sigw_ng(tmpl[shape], fsx, LOG10_KSTAR, LOG10_A, F, G, order)
            lpa = log_posterior_sigw_ng(tmpl[shape], fs, LOG10_KSTAR, LOG10_A, F, G, alt)
            med[(shape, template)] = posterior_medians(lp)
            if (shape, template) in NG_FIGURE_CASES:
                figs[f"{shape}|{template}"] = compact(lp[::2, ::2], 1e-3)
            for key, line in sorted(lines.items()):
                if key[:2] != (shape, template):
                    continue
                level, mass, dchi2 = tension_metrics(lp, line)
                _, mass_x, _ = tension_metrics(lpx, line)
                _, mass_a, _ = tension_metrics(lpa, line)
                g = gauss[key]
                rows.append(dict(shape=shape, template=template, statistic=key[2], variant=key[3], f_pbh=key[4],
                                 F=round(F, 6), G=round(G, 6), order=order,
                                 sigma_p_gsigw=float(g["sigma_p"]), sigma_p=sigma_one_sided(mass), sigma=sigma_1d(level),
                                 p_allowed=mass, dchi2=dchi2, sigma_p_ceffyl=sigma_one_sided(mass_x),
                                 sigma_p_alt=sigma_one_sided(mass_a), alt_order=alt))
            print(f"{shape:10s} {template:15s} F={F:+.3f} G={G:.3f} A^{order} ({time.time() - t0:.0f}s)", flush=True)
    np.savez_compressed("results/phase5_ngsigw_posteriors.npz", log10_kstar=LOG10_KSTAR[::2], log10_A=LOG10_A[::2], **figs)
    with open(out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)
    write_ngsigw_summary(rows, med)


def posterior_medians(logpost):
    """1D posterior medians of log10 k* and log10 A."""
    p = np.exp(logpost - logpost.max())

    def median(w, x):
        return float(np.interp(0.5, np.cumsum(w) / w.sum(), x))

    return median(p.sum(1), LOG10_KSTAR), median(p.sum(0), LOG10_A)


def write_ngsigw_summary(rows, med, out="results/phase5_ngsigw_summary.md"):
    idx = {(r["shape"], r["template"], r["statistic"], r["variant"], r["f_pbh"]): r for r in rows}
    md = ["# Phase 5 with the primordial-NG corrections to the SIGW spectrum", "",
          "Each NG template gets its own NG15 posterior: Omega_GW = A^2 T2 + A^3 [F^2 T3 + 12 G sbar T2] + F^4 A^4 T4 "
          "(exact for the quadratic template; truncated at A^3 for the curvaton and USR maps). Cells: sigma_P with the "
          "Gaussian SIGW (as in phase5_summary.md) -> with the NG SIGW. 'alt': the A^3 truncation for the quadratic "
          "template, the A^3 series plus the F^4 A^4 terms for the others. Values above "
          f"{SIGMA_RELIABLE:.0f} sigma lie beyond the free-spectrum KDE samples.", "",
          "## Posterior medians (1D marginals)", "",
          "| shape | template | F | G | log10 k*: Gaussian SIGW -> NG SIGW | log10 A |", "|---|---|---|---|---|---|"]
    for shape in SHAPES:
        k0, a0 = med[(shape, "gaussian")]
        for t in TEMPLATES:
            if t == "gaussian":
                continue
            r = next(x for x in rows if x["shape"] == shape and x["template"] == t)
            k1, a1 = med[(shape, t)]
            md.append(f"| {shape} | {t} | {r['F']:+.3f} | {r['G']:.3f} | {k0:.3f} -> {k1:.3f} | "
                      f"{a0:.3f} -> {a1:.3f} ({a1 - a0:+.3f}) |")
    md.append("")
    for variant, f in [("base", 1.0), ("base", 1e-3), ("T", 1.0), ("entropy", 1.0)]:
        templates = [t for t in (ENTROPY_TEMPLATES if variant == "entropy" else TEMPLATES) if t != "gaussian"]
        md += [f"## variant {variant}, f_PBH = {f:g}", "",
               "| shape | template | threshold | peaks 3a | peaks 3b | NG-SIGW shift (thr, 3b) | "
               "3b - threshold: G-SIGW -> NG-SIGW | alt (thr, 3b) |",
               "|---|---|---|---|---|---|---|---|"]
        for shape in SHAPES:
            for t in templates:
                s = [idx.get((shape, t, st, variant, f)) for st in STATISTICS]
                if None in s:
                    continue
                cells = [f"{fmt_sigma(x['sigma_p_gsigw'])} -> {fmt_sigma(x['sigma_p'])}" for x in s]
                d_thr = s[0]["sigma_p"] - s[0]["sigma_p_gsigw"]
                d_3b = s[2]["sigma_p"] - s[2]["sigma_p_gsigw"]
                stat_g = s[2]["sigma_p_gsigw"] - s[0]["sigma_p_gsigw"]
                stat_n = s[2]["sigma_p"] - s[0]["sigma_p"]
                alt = f"{fmt_sigma(s[0]['sigma_p_alt'])}, {fmt_sigma(s[2]['sigma_p_alt'])}"
                md.append(f"| {shape} | {t} | {cells[0]} | {cells[1]} | {cells[2]} | {d_thr:+.2f}, {d_3b:+.2f} | "
                          f"{stat_g:+.2f} -> {stat_n:+.2f} | {alt} |")
        md.append("")
    Path(out).write_text("\n".join(md))
    print("\n".join(md))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("stage", nargs="?", default="all", choices=["lines", "tension", "compare", "ngsigw", "all"])
    ap.add_argument("--workers", type=int, default=4)
    args = ap.parse_args()
    if args.stage in ("lines", "all"):
        stage_lines(args.workers)
    if args.stage in ("tension", "all"):
        stage_tension()
    if args.stage in ("compare", "all"):
        stage_compare()
    if args.stage in ("ngsigw", "all"):
        stage_ngsigw()


if __name__ == "__main__":
    main()
