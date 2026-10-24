"""Phase 5 figures from results/phase5_*.{csv,npz} and the digitised 2306.17149 Fig. 3.

fig_phase5_tension.png    NG15 posterior with the overproduction lines per statistic (BPL and LN D=1)
fig_phase5_sigma.png      exclusion significance per NG template and statistic, per spectral shape
fig_phase5_validation.png our threshold-statistics bands and posterior against 2306.17149 Fig. 3
fig_sigw_ng_terms.png     NG corrections to the SIGW spectrum per spectral shape, T3/T2 and T4/T2
fig_phase5_ngsigw.png     the NG15 posterior with the NG-corrected SIGW, and the change of sigma_P

Same palette rules as make_figures.py: f_NL sign on a diverging scale with a gray Gaussian
midpoint; the USR/curvaton templates take categorical slots; statistics are told apart by line
style and marker, so colour never carries the statistic alone.
"""

import csv
from collections import defaultdict
from pathlib import Path

import numpy as np

from make_figures import AXIS, CATEGORICAL, DIVERGING, GRID, INK2, MUTED, SURFACE, plt

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
                sig = float(tension[(shape, t, st, "base", 1.0)]["sigma_p"])
                sig_txt = ">4" if sig > 4 else f"{sig:.1f}"
                ax.plot(line[:, 0], line[:, 1], ls, color=NG_COLOURS[t], linewidth=1.5,
                        label=f"{NG_LABELS[t]}, {slab}: {sig_txt}$\\sigma$")
        axes_style(ax)
        ax.set_xlim(6.0, 9.5)
        ax.set_ylim(-2.6, 0.0)
        ax.set_title(f"{shape.replace('D=', 'Δ=')}: NG15 posterior and $f_{{PBH}}=1$ lines", loc="left")
        ax.legend(loc="upper left", bbox_to_anchor=(0.0, 0.93), fontsize=6.3, handlelength=2.6)
        ax.text(0.02, 0.98, "NG15 posterior: 68.3 / 95.4 / 99.7% HPD;  $\\sigma=\\sigma_P$", transform=ax.transAxes,
                fontsize=6.8, color=MUTED, va="top")
    fig.tight_layout()
    fig.savefig(out, dpi=200)
    plt.close(fig)


def fig_sigma(tension, out, variant="base", f=1.0, cap=4.0):
    shapes = ["BPL(4,3,1)", "LN D=0.5", "LN D=1", "LN D=2"]
    templates = ["fNL=-2", "fNL=-1", "fNL=-0.5", "gaussian", "fNL=+0.5", "fNL=+1", "fNL=+2",
                 "curvaton r=0.9", "USR mu=5/2", "USR mu=5/6", "USR beta=3"]
    fig, axes = plt.subplots(1, 4, figsize=(10.0, 3.8), sharey=True)
    y = np.arange(len(templates))[::-1]
    for ax, shape in zip(axes, shapes):
        for yy, t in zip(y, templates):
            vals = [min(float(tension[(shape, t, st, variant, f)]["sigma_p"]), cap + 0.15) for st in STAT_STYLE]
            ax.plot([min(vals), max(vals)], [yy, yy], color=AXIS, linewidth=1.0, zorder=1)
            for (st, (_, mk, slab)), v in zip(STAT_STYLE.items(), vals):
                ax.plot(v, yy, mk, color=NG_COLOURS[t], markeredgecolor=SURFACE if st != "peaks_3a" else INK2,
                        markeredgewidth=0.6, zorder=2, label=slab if yy == y[0] else None,
                        markerfacecolor=NG_COLOURS[t] if st != "peaks_3a" else SURFACE)
        for s in (2.0, 3.0):
            ax.axvline(s, color=MUTED, linewidth=0.8, linestyle=(0, (3, 3)), zorder=0)
        ax.axvspan(cap, cap + 0.4, color=GRID, zorder=0)
        ax.set_xlim(0.8, cap + 0.4)
        ax.text(cap + 0.2, y[-1] - 0.9, f">{cap:.0f}", fontsize=7, color=MUTED, ha="center")
        ax.grid(True, axis="x", color=GRID, linewidth=0.6)
        for side in ("top", "right"):
            ax.spines[side].set_visible(False)
        ax.set_title(shape.replace("D=", "Δ="), loc="left")
        ax.set_xlabel(r"exclusion of $f_{PBH}<1$, $\sigma_P$")
    axes[0].set_yticks(y)
    axes[0].set_yticklabels([NG_LABELS[t] for t in templates], fontsize=7.5)
    handles = [plt.Line2D([], [], marker=mk, linestyle="", color=INK2, markerfacecolor=INK2 if st != "peaks_3a" else SURFACE,
                          label=slab) for st, (_, mk, slab) in STAT_STYLE.items()]
    fig.legend(handles=handles, loc="upper right", ncol=3, fontsize=7.5)
    fig.suptitle("NG15 exclusion of the PBH-allowed region per statistic (posterior mass; 2306.17149 conventions)",
                 x=0.01, ha="left", fontsize=9.5, color=INK2)
    fig.tight_layout(rect=(0, 0, 1, 0.95))
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
    ax.legend(loc="center right", bbox_to_anchor=(1.0, 0.36), fontsize=7)
    ax = axes[1]
    if marg is not None:
        posterior_panel(ax, marg["log10_kstar"], marg["log10_A"], marg["logpost"])
        title = "NG15 posterior, BPL shape marginalised"
    else:
        posterior_panel(ax, post["log10_kstar"], post["log10_A"], post["BPL(4,3,1)"])
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


SHAPE_COLOURS = dict(zip(["BPL(4,3,1)", "LN D=0.5", "LN D=1", "LN D=2"], CATEGORICAL))


def fig_ngsigw_terms(out, path="results/sigw_ng_terms.csv"):
    from ngpbh.sigw_ng import load_terms

    terms = load_terms(path)
    fig, axes = plt.subplots(1, 2, figsize=(9.0, 3.6), sharey=True)
    band = (1.28e6 / 10**7.8, 1.79e7 / 10**7.8)  # the 14 NG15 bins for k* = 10^7.8 Mpc^-1
    for ax, (key, title) in zip(axes, [("T3", r"$O(F^2A^3)$: (hybrid + C + Z) / $T_2$"),
                                       ("T4", r"$O(F^4A^4)$: (reducible + planar + nonplanar) / $T_2$")]):
        ax.axvspan(*band, color=GRID, alpha=0.6, zorder=0, linewidth=0)
        for shape, d in terms.items():
            num = d["hybrid"] + d["C"] + d["Z"] if key == "T3" else d["reducible"] + d["planar"] + d["nonplanar"]
            ax.plot(d["x"], num / d["T2"], color=SHAPE_COLOURS[shape], label=shape.replace("D=", "Δ="))
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.grid(True, color=GRID, linewidth=0.6)
        ax.set_axisbelow(True)
        for side in ("top", "right"):
            ax.spines[side].set_visible(False)
        ax.set_xlabel(r"$k/k_*$")
        ax.set_title(title, loc="left")
    axes[0].set_ylabel(r"ratio to the Gaussian term $T_2$")
    axes[0].text(band[0] * 1.1, 0.06, "NG15 bins,\n$k_*=10^{7.8}$", fontsize=7, color=MUTED)
    axes[0].legend(loc="upper left", fontsize=7)
    fig.suptitle(r"$\Omega_{GW}=A^2T_2+F^2A^3T_3+F^4A^4T_4$ for "
                 r"$\zeta=\zeta_g+F(\zeta_g^2-\langle\zeta_g^2\rangle)$, $F=3f_{NL}/5$",
                 x=0.01, ha="left", fontsize=9.5, color=INK2)
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    fig.savefig(out, dpi=200)
    plt.close(fig)


def fig_ngsigw(lines, post, out, path="results/phase5_ngsigw_tension.csv", npz="results/phase5_ngsigw_posteriors.npz",
               cap=4.0):
    ng = dict(np.load(npz))
    lk, la = ng["log10_kstar"], ng["log10_A"]
    fig, axes = plt.subplots(1, 2, figsize=(9.6, 4.5), gridspec_kw={"width_ratios": [1.0, 1.15]})
    ax = axes[0]
    posterior_panel(ax, post["log10_kstar"], post["log10_A"], post["BPL(4,3,1)"])
    p = np.exp(ng["BPL(4,3,1)|fNL=-2"] - ng["BPL(4,3,1)|fNL=-2"].max())
    p /= p.sum()
    ax.contour(lk, la, p.T, levels=hpd_levels(ng["BPL(4,3,1)|fNL=-2"]), colors=[INK2], linewidths=0.9,
               linestyles=[(0, (4, 2))])
    for st in ("threshold", "peaks_3b"):
        ls, _, slab = STAT_STYLE[st]
        line = lines[("BPL(4,3,1)", "fNL=-2", st, "base", 1.0)]
        ax.plot(line[:, 0], line[:, 1], ls, color=NG_COLOURS["fNL=-2"], linewidth=1.5, label=f"$f_{{PBH}}=1$, {slab}")
    ax.plot([], [], color=INK2, linewidth=0.9, linestyle=(0, (4, 2)), label="NG15 posterior with the NG SIGW")
    axes_style(ax)
    ax.set_xlim(6.0, 9.5)
    ax.set_ylim(-2.6, 0.0)
    ax.set_title(r"BPL(4,3,1), $f_{NL}=-2$: Gaussian SIGW (fill) vs NG SIGW", loc="left")
    ax.legend(loc="lower right", fontsize=6.8)
    ax = axes[1]
    rows = [r for r in csv.DictReader(open(path)) if r["variant"] == "base" and float(r["f_pbh"]) == 1.0]
    templates = ["fNL=-2", "fNL=-1", "fNL=-0.5", "fNL=+0.5", "fNL=+1", "fNL=+2",
                 "curvaton r=0.9", "USR mu=5/2", "USR mu=5/6", "USR beta=3"]
    y = np.arange(len(templates))[::-1]
    offsets = dict(zip(SHAPE_COLOURS, (-0.27, -0.09, 0.09, 0.27)))
    for r in rows:
        if r["statistic"] not in ("threshold", "peaks_3b"):
            continue
        yy = y[templates.index(r["template"])] + offsets[r["shape"]]
        d = float(r["sigma_p"]) - float(r["sigma_p_gsigw"])
        reliable = max(float(r["sigma_p"]), float(r["sigma_p_gsigw"])) <= cap
        mk = "o" if r["statistic"] == "threshold" else "s"
        col = SHAPE_COLOURS[r["shape"]]
        ax.plot(d, yy, mk, markersize=4.2, color=col, markerfacecolor=col if reliable else SURFACE,
                markeredgecolor=col if not reliable else SURFACE, markeredgewidth=0.8)
    ax.axvline(0.0, color=AXIS, linewidth=0.9, zorder=0)
    ax.set_yticks(y)
    ax.set_yticklabels([NG_LABELS[t] for t in templates], fontsize=7.5)
    ax.grid(True, axis="x", color=GRID, linewidth=0.6)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    ax.set_xlabel(r"$\Delta\sigma_P$ = NG SIGW $-$ Gaussian SIGW")
    ax.set_title(r"change of the exclusion, $f_{PBH}=1$ (2306.17149 conventions)", loc="left")
    handles = [plt.Line2D([], [], marker="s", linestyle="", color=c, label=s.replace("D=", "Δ="))
               for s, c in SHAPE_COLOURS.items()]
    handles += [plt.Line2D([], [], marker="o", linestyle="", color=INK2, label="threshold"),
                plt.Line2D([], [], marker="s", linestyle="", color=INK2, label="peaks 3b"),
                plt.Line2D([], [], marker="o", linestyle="", color=INK2, markerfacecolor=SURFACE,
                           label=f"beyond {cap:.0f}σ")]
    ax.legend(handles=handles, loc="upper center", bbox_to_anchor=(0.45, -0.16), fontsize=6.8, ncol=4)
    fig.tight_layout()
    fig.savefig(out, dpi=200)
    plt.close(fig)


def main():
    lines = read_lines()
    tension = read_tension()
    post = dict(np.load("results/phase5_posteriors.npz"))
    marg_path = Path("results/phase5_posterior_bpl_marginal.npz")
    marg = dict(np.load(marg_path)) if marg_path.exists() else None
    fig_tension(lines, tension, post, "results/fig_phase5_tension.png")
    fig_sigma(tension, "results/fig_phase5_sigma.png")
    fig_validation(lines, post, marg, read_digitized(), "results/fig_phase5_validation.png")
    if Path("results/sigw_ng_terms.csv").exists():
        fig_ngsigw_terms("results/fig_sigw_ng_terms.png")
    if Path("results/phase5_ngsigw_tension.csv").exists():
        fig_ngsigw(lines, post, "results/fig_phase5_ngsigw.png")


if __name__ == "__main__":
    main()
