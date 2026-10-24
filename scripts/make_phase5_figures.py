"""Phase 5 figures from results/phase5_*.{csv,npz} and the digitised 2306.17149 Fig. 3.

fig_phase5_tension.png    NG15 posterior with the overproduction lines per statistic (BPL and LN D=1)
fig_phase5_sigma.png      exclusion significance per NG template and statistic, per spectral shape
fig_phase5_validation.png our threshold-statistics bands and posterior against 2306.17149 Fig. 3

Same palette rules as make_figures.py: f_NL sign on a diverging scale with a gray Gaussian
midpoint; the USR/curvaton templates take categorical slots; statistics are told apart by line
style and marker, so colour never carries the statistic alone.
"""

import csv
from collections import defaultdict
from pathlib import Path

import numpy as np

from make_figures import AXIS, DIVERGING, GRID, INK2, MUTED, SURFACE, plt

NG_COLOURS = dict(DIVERGING)
NG_COLOURS.update({"USR beta=3": "#1baf7a", "USR mu=5/6": "#eda100", "USR mu=5/2": "#eb6834", "curvaton r=0.9": "#7c5cc4"})
NG_LABELS = {
    "gaussian": "Gaussian", "fNL=-2": r"$f_{NL}=-2$", "fNL=-1": r"$f_{NL}=-1$", "fNL=-0.5": r"$f_{NL}=-0.5$",
    "fNL=+0.5": r"$f_{NL}=+0.5$", "fNL=+1": r"$f_{NL}=+1$", "fNL=+2": r"$f_{NL}=+2$",
    "USR beta=3": r"USR $\beta=3$", "USR mu=5/6": r"USR $\mu_*=5/6$", "USR mu=5/2": r"USR $\mu_*=5/2$",
    "curvaton r=0.9": r"curvaton $r=0.9$",
}
STAT_STYLE = {"threshold": ("-", "o", "threshold"), "peaks_3a": (":", "^", "peaks 3a"), "peaks_3b": ("--", "s", "peaks 3b")}
POST_FILL = ["#d6e6f8", "#9cc3ef", "#5b9be2"]  # 99.7, 95.4, 68.3% HPD, one sequential hue


def read_lines(path="results/phase5_lines.csv"):
    out = defaultdict(list)
    with open(path) as fh:
        for r in csv.DictReader(fh):
            out[(r["shape"], r["template"], r["statistic"], r["variant"], float(r["f_pbh"]))].append(
                (float(r["log10_kstar"]), float(r["log10_A_star"])))
    return {k: np.array(sorted(v)) for k, v in out.items()}


def read_tension(path="results/phase5_tension.csv"):
    with open(path) as fh:
        return {(r["shape"], r["template"], r["statistic"], r["variant"], float(r["f_pbh"])): r
                for r in csv.DictReader(fh)}


def hpd_levels(logpost, levels=(0.997, 0.954, 0.683)):
    p = np.exp(logpost - logpost.max())
    p /= p.sum()
    s = np.sort(p.ravel())[::-1]
    c = np.cumsum(s)
    return [s[np.searchsorted(c, lv)] for lv in levels]


def posterior_panel(ax, lk, la, logpost):
    p = np.exp(logpost - logpost.max())
    p /= p.sum()
    lv = hpd_levels(logpost)
    ax.contourf(lk, la, p.T, levels=[lv[0], lv[1], lv[2], p.max() * 1.01], colors=POST_FILL)
    ax.contour(lk, la, p.T, levels=lv, colors=["#2a78d6"], linewidths=0.6)


def axes_style(ax):
    ax.grid(True, color=GRID, linewidth=0.6)
    ax.set_axisbelow(True)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    ax.set_xlabel(r"$\log_{10}(k_*/{\rm Mpc}^{-1})$")
    ax.set_ylabel(r"$\log_{10} A$")


def fig_tension(lines, tension, post, out):
    lk, la = post["log10_kstar"], post["log10_A"]
    fig, axes = plt.subplots(1, 2, figsize=(9.0, 3.9), sharey=True)
    for ax, shape in zip(axes, ["BPL(4,3,1)", "LN D=1"]):
        posterior_panel(ax, lk, la, post[shape])
        for t in ["fNL=-2", "gaussian", "USR beta=3"]:
            for st in ["threshold", "peaks_3b"]:
                ls, _, slab = STAT_STYLE[st]
                line = lines[(shape, t, st, "base", 1.0)]
                sig = float(tension[(shape, t, st, "base", 1.0)]["sigma"])
                ax.plot(line[:, 0], line[:, 1], ls, color=NG_COLOURS[t], linewidth=1.5,
                        label=f"{NG_LABELS[t]}, {slab}: {sig:.1f}$\\sigma$")
        axes_style(ax)
        ax.set_xlim(6.0, 9.5)
        ax.set_ylim(-2.6, 0.0)
        ax.set_title(f"{shape.replace('D=', 'Δ=')}: NG15 posterior and $f_{{PBH}}=1$ lines", loc="left")
        ax.legend(loc="lower right", fontsize=6.3, handlelength=2.6)
    if axes[0].get_legend() is not None:
        axes[0].text(0.02, 0.97, "fill: 68.3 / 95.4 / 99.7% HPD", transform=axes[0].transAxes, fontsize=7,
                     color=MUTED, va="top")
    fig.tight_layout()
    fig.savefig(out, dpi=200)
    plt.close(fig)


def fig_sigma(tension, out, variant="base", f=1.0):
    shapes = ["BPL(4,3,1)", "LN D=0.5", "LN D=1", "LN D=2"]
    templates = ["fNL=-2", "fNL=-1", "fNL=-0.5", "gaussian", "fNL=+0.5", "fNL=+1", "fNL=+2",
                 "curvaton r=0.9", "USR mu=5/2", "USR mu=5/6", "USR beta=3"]
    fig, axes = plt.subplots(1, 4, figsize=(10.0, 3.8), sharey=True)
    y = np.arange(len(templates))[::-1]
    for ax, shape in zip(axes, shapes):
        for yy, t in zip(y, templates):
            vals = [float(tension[(shape, t, st, variant, f)]["sigma"]) for st in STAT_STYLE]
            ax.plot([min(vals), max(vals)], [yy, yy], color=AXIS, linewidth=1.0, zorder=1)
            for (st, (_, mk, slab)), v in zip(STAT_STYLE.items(), vals):
                ax.plot(v, yy, mk, color=NG_COLOURS[t], markeredgecolor=SURFACE if st != "peaks_3a" else INK2,
                        markeredgewidth=0.6, zorder=2, label=slab if yy == y[0] else None,
                        markerfacecolor=NG_COLOURS[t] if st != "peaks_3a" else SURFACE)
        for s in (2.0, 3.0):
            ax.axvline(s, color=MUTED, linewidth=0.8, linestyle=(0, (3, 3)), zorder=0)
        ax.grid(True, axis="x", color=GRID, linewidth=0.6)
        for side in ("top", "right"):
            ax.spines[side].set_visible(False)
        ax.set_title(shape.replace("D=", "Δ="), loc="left")
        ax.set_xlabel(r"exclusion of $f_{PBH}<1$ [$\sigma$]")
    axes[0].set_yticks(y)
    axes[0].set_yticklabels([NG_LABELS[t] for t in templates], fontsize=7.5)
    handles = [plt.Line2D([], [], marker=mk, linestyle="", color=INK2, markerfacecolor=INK2 if st != "peaks_3a" else SURFACE,
                          label=slab) for st, (_, mk, slab) in STAT_STYLE.items()]
    fig.legend(handles=handles, loc="upper right", ncol=3, fontsize=7.5)
    fig.suptitle("NG15 exclusion of the PBH-allowed region, per statistic (2306.17149 setup)", x=0.01, ha="left",
                 fontsize=9.5, color=INK2)
    fig.tight_layout(rect=(0, 0, 1, 0.93))
    fig.savefig(out, dpi=200)
    plt.close(fig)


def fig_validation(lines, post, marg, digitized, out):
    fig, axes = plt.subplots(1, 2, figsize=(9.0, 3.9))
    ax = axes[0]
    for t in ["fNL=-2", "curvaton r=0.9", "gaussian", "USR beta=3"]:
        c = NG_COLOURS[t]
        for f, ls in [(1.0, "-"), (1e-3, "-")]:
            ours = lines[("BPL(4,3,1)", t, "threshold", "base", f)]
            theirs = digitized[(t, "f=1" if f == 1.0 else "f=1e-3")]
            ax.plot(theirs[:, 0], theirs[:, 1], color=c, linewidth=3.2, alpha=0.35, solid_capstyle="butt")
            ax.plot(ours[:, 0], ours[:, 1], "--", color=c, linewidth=1.1,
                    label=NG_LABELS[t] if f == 1.0 else None)
    axes_style(ax)
    ax.set_xlim(6.0, 9.0)
    ax.set_ylim(-2.35, -1.25)
    ax.set_title("BPL(4,3,1), threshold statistics: $f_{PBH}=1$ and $10^{-3}$", loc="left")
    ax.text(0.02, 0.03, "wide: 2306.17149 Fig. 3 (with QCD)   dashed: ngpbh (no QCD)", transform=ax.transAxes,
            fontsize=6.8, color=MUTED)
    ax.legend(loc="upper right", fontsize=7)
    ax = axes[1]
    lk, la = post["log10_kstar"], post["log10_A"]
    if marg is not None:
        posterior_panel(ax, lk, la, marg)
        title = "NG15 posterior, BPL shape marginalised"
    else:
        posterior_panel(ax, lk, la, post["BPL(4,3,1)"])
        title = "NG15 posterior, BPL(4,3,1)"
    for edge, lw in [("1sigma", 1.0), ("2sigma", 1.0), ("3sigma", 1.0)]:
        d = digitized[("NG15", edge)]
        ax.plot(d[:, 0], d[:, 1], color="#1d6b2f", linewidth=lw, label="2306.17149 Fig. 3" if edge == "1sigma" else None)
    axes_style(ax)
    ax.set_xlim(6.0, 9.5)
    ax.set_ylim(-2.6, 0.0)
    ax.set_title(title, loc="left")
    ax.legend(loc="lower right", fontsize=7)
    fig.tight_layout()
    fig.savefig(out, dpi=200)
    plt.close(fig)


def read_digitized(path="results/digitized_2306.17149_fig3.csv"):
    out = defaultdict(list)
    with open(path) as fh:
        for r in csv.DictReader(fh):
            out[(r["curve"], r["edge"])].append((float(r["log10_kstar"]), float(r["log10_A"])))
    return {k: np.array(v) for k, v in out.items()}


def main():
    lines = read_lines()
    tension = read_tension()
    post = dict(np.load("results/phase5_posteriors.npz"))
    marg_path = Path("results/phase5_posterior_bpl_marginal.npz")
    marg = np.load(marg_path)["logpost"] if marg_path.exists() else None
    fig_tension(lines, tension, post, "results/fig_phase5_tension.png")
    fig_sigma(tension, "results/fig_phase5_sigma.png")
    fig_validation(lines, post, marg, read_digitized(), "results/fig_phase5_validation.png")


if __name__ == "__main__":
    main()
