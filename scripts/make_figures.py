"""Figures from results/*.csv: PNGs for the repository and the note, and vector PDFs at print size
(the JCAP text width) for the paper, in docs/paper/figures.

Colour: the sign of f_NL is a polarity, so the quadratic family uses a diverging scale
(blue arm f_NL < 0, red arm f_NL > 0, gray Gaussian midpoint; each arm validated as an ordinal
ramp); the logarithmic family and the convention comparison use the first categorical slots.
Light steps sit below 3:1 contrast, so lines are direct-labelled as well as in the legend.
"""

import csv
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.ticker import FormatStrFormatter  # noqa: E402

INK, INK2, MUTED, GRID, AXIS = "#0b0b0b", "#52514e", "#898781", "#e1e0d9", "#c3c2b7"
SURFACE = "#ffffff"
TEXTWIDTH = 6.09  # inches: the JCAP text width, 440 pt
PAPER_DIR = Path("docs/paper/figures")
GAUSS = INK2
DIVERGING = {  # f_NL -> colour
    "fNL=-2": "#184f95",
    "fNL=-1": "#2a78d6",
    "fNL=-0.5": "#6da7ec",
    "gaussian": GAUSS,
    "fNL=+0.5": "#f08a85",
    "fNL=+1": "#e34948",
    "fNL=+2": "#a8302f",
}
CATEGORICAL = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100"]
LOG_FAMILY = ["USR mu=5/6", "USR mu=5/2", "USR mu=10", "log mu=-5/2"]
LABELS = {
    "gaussian": "Gaussian",
    "fNL=-2": r"$f_{\rm NL}=-2$",
    "fNL=-1": r"$f_{\rm NL}=-1$",
    "fNL=-0.5": r"$f_{\rm NL}=-0.5$",
    "fNL=+0.5": r"$f_{\rm NL}=+0.5$",
    "fNL=+1": r"$f_{\rm NL}=+1$",
    "fNL=+2": r"$f_{\rm NL}=+2$",
    "USR mu=5/6": r"$\mu_*=5/6$",
    "USR mu=5/2": r"$\mu_*=5/2$",
    "USR mu=10": r"$\mu_*=10$",
    "log mu=-5/2": r"$\mu_*=-5/2$",
}

plt.rcParams.update(
    {
        "figure.facecolor": SURFACE,
        "axes.facecolor": SURFACE,
        "axes.edgecolor": AXIS,
        "axes.labelcolor": INK2,
        "axes.titlecolor": INK,
        "axes.titlesize": 8.5,
        "axes.labelsize": 8.5,
        "axes.linewidth": 0.6,
        "xtick.color": MUTED,
        "ytick.color": MUTED,
        "xtick.labelsize": 7.5,
        "ytick.labelsize": 7.5,
        "xtick.major.width": 0.6,
        "ytick.major.width": 0.6,
        "xtick.minor.width": 0.4,
        "ytick.minor.width": 0.4,
        "legend.fontsize": 6.8,
        "legend.frameon": False,
        "font.family": "sans-serif",
        "lines.linewidth": 1.3,
        "lines.markersize": 3.6,
        "pdf.fonttype": 42,  # embed TrueType, not Type 3, in the vector figures
    }
)


def read(path):
    with open(path) as fh:
        return list(csv.DictReader(fh))


def save(fig, out, paper=True):
    """The PNG for the repository and, for a figure of the paper, the vector PDF at print size."""
    fig.savefig(out, dpi=200)
    if paper:
        PAPER_DIR.mkdir(parents=True, exist_ok=True)
        fig.savefig(PAPER_DIR / (Path(out).stem + ".pdf"))
    plt.close(fig)


def style(ax):
    ax.grid(True, color=GRID, linewidth=0.6)
    ax.set_axisbelow(True)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    ax.set_xscale("log")
    ax.set_xticks([0.1, 0.3, 0.5, 1.0, 2.0])
    ax.set_xticklabels(["0.1", "0.3", "0.5", "1", "2"])
    ax.set_xlabel(r"lognormal width $\Delta$")


def end_labels(ax, items, min_gap_pt=9.0):
    """Direct labels at line ends, pushed apart vertically so that none overlap.

    items: (x, y, text) in data coordinates; offsets are solved in display points."""
    fig = ax.figure
    fig.canvas.draw()
    to_pt = 72.0 / fig.dpi
    ys = [ax.transData.transform((x, y))[1] * to_pt for x, y, _ in items]
    order = sorted(range(len(items)), key=lambda i: ys[i])
    placed = {}
    last = -1e9
    for i in order:
        placed[i] = max(ys[i], last + min_gap_pt)
        last = placed[i]
    shift = (max(placed[i] - ys[i] for i in placed) / 2.0) if placed else 0.0
    for i, (x, y, text) in enumerate(items):
        dy = placed[i] - ys[i] - shift if len(items) > 1 else 0.0
        ax.annotate(text, (x, y), xytext=(5, dy), textcoords="offset points", fontsize=7, color=INK2, va="center")


def series(rows, template, statistic, cut="4.49"):
    pts = sorted(
        (float(r["width"]), float(r["A_star"]))
        for r in rows
        if r["template"] == template and r["statistic"] == statistic and r["x_cut"] == cut
    )
    return [p[0] for p in pts], [p[1] for p in pts]


def ratio(rows, template, peaks="peaks"):
    w, a_pk = series(rows, template, peaks)
    _, a_th = series(rows, template, "threshold")
    return w, [p / t for p, t in zip(a_pk, a_th)]


def correction(rows, template):
    """A_* of Phase 3b peaks over Phase 3a peaks."""
    w, a_b = series(rows, template, "peaks_3b")
    _, a_a = series(rows, template, "peaks")
    return w, [b / a for b, a in zip(a_b, a_a)]


def fig_amplitudes(rows, out):
    fig, axes = plt.subplots(1, 2, figsize=(9.0, 3.8), sharey=True)
    ax = axes[0]
    for t, c in DIVERGING.items():
        w, a = series(rows, t, "threshold")
        ax.plot(w, a, "-o", color=c, label=LABELS[t], markeredgecolor=SURFACE, markeredgewidth=0.8)
    end_labels(ax, [(series(rows, t, "threshold")[0][-1], series(rows, t, "threshold")[1][-1], LABELS[t])
                    for t in ("fNL=-2", "gaussian", "fNL=+2")])
    ax.set_title("Quadratic template", loc="left")
    ax.set_ylabel(r"$A_*$ for $f_{PBH}=1$ (threshold statistics)")
    ax.legend(loc="upper left", ncol=2)
    ax = axes[1]
    w, a = series(rows, "gaussian", "threshold")
    ax.plot(w, a, "-o", color=GAUSS, label="Gaussian", markeredgecolor=SURFACE, markeredgewidth=0.8)
    labels = []
    for t, c in zip(LOG_FAMILY, CATEGORICAL):
        w, a = series(rows, t, "threshold")
        ax.plot(w, a, "-o", color=c, label=LABELS[t], markeredgecolor=SURFACE, markeredgewidth=0.8)
        labels.append((w[-1], a[-1], LABELS[t]))
    ax.set_title(r"Logarithmic (USR) template $\zeta=-\mu_*\ln(1-\zeta_g/\mu_*)$", loc="left")
    ax.legend(loc="upper left")
    for a_ in axes:
        style(a_)
        a_.set_yscale("log")
        a_.set_xlim(0.085, 3.4)
    end_labels(ax, labels)
    fig.suptitle(
        r"Phase 2: required amplitude, lognormal $P_\zeta$ at $k_*=10^7\,$Mpc$^{-1}$, Musco threshold, transfer function, $k<4.49/r$",
        fontsize=9,
        color=INK2,
        x=0.01,
        ha="left",
    )
    fig.tight_layout()
    save(fig, out, paper=False)


def fig_ratio(rows, out):
    """R for maxima of C_l (solid) and of C_g (dashed, same colour), and the ratio of the two amplitudes."""
    fig, axes = plt.subplots(1, 3, figsize=(TEXTWIDTH, 2.6), layout="constrained")
    families = [
        (axes[0], list(DIVERGING.items()), "quadratic"),
        (axes[1], [("gaussian", GAUSS)] + list(zip(LOG_FAMILY, CATEGORICAL)), "logarithmic (USR)"),
    ]
    for ax, members, title in families:
        for t, c in members:
            w, r3a = ratio(rows, t, "peaks")
            ax.plot(w, r3a, "--", color=c, linewidth=0.8, alpha=0.7)
            w, r3b = ratio(rows, t, "peaks_3b")
            ax.plot(w, r3b, "-o", color=c, label=LABELS[t], markeredgecolor=SURFACE, markeredgewidth=0.6)
        ax.set_title(title, loc="left")
        style(ax)
        ax.set_ylim(0.72, 0.94)
        ax.set_yticks([0.75, 0.80, 0.85, 0.90])
        ax.yaxis.set_major_formatter(FormatStrFormatter("%.2f"))
        ax.set_xlim(0.085, 2.4)
        ax.legend(loc="upper left", ncol=2 if len(members) > 5 else 1, fontsize=6.2, handlelength=1.5,
                  columnspacing=0.8, borderaxespad=0.3, labelspacing=0.3)
    axes[0].set_ylabel(r"$R=A_*^{\rm peaks}/A_*^{\rm threshold}$")
    # the line-style key sits apart from the colour key, in a free corner of the middle panel
    key = [plt.Line2D([], [], color=MUTED, linewidth=1.1, label=r"maxima of $C_\ell$"),
           plt.Line2D([], [], color=MUTED, linewidth=0.8, linestyle="--", label=r"maxima of $C_{\rm g}$")]
    axes[1].add_artist(axes[1].get_legend())
    axes[1].legend(handles=key, loc="lower right", fontsize=6.2, handlelength=2.0, borderaxespad=0.3)
    ax = axes[2]
    # colours as in the other panels; the logarithmic family is dashed with square markers here
    for t, c in DIVERGING.items():
        if t != "gaussian":
            w, q = correction(rows, t)
            ax.plot(w, q, "-o", color=c, markeredgecolor=SURFACE, markeredgewidth=0.6)
    for t, c in zip(LOG_FAMILY, CATEGORICAL):
        w, q = correction(rows, t)
        ax.plot(w, q, "--s", color=c, markeredgecolor=SURFACE, markeredgewidth=0.6)
    ax.axhline(1.0, color=GAUSS, linewidth=0.8)
    style(ax)
    ax.set_xlim(0.085, 2.4)
    ax.set_yticks([0.985, 0.990, 0.995, 1.000])
    ax.yaxis.set_major_formatter(FormatStrFormatter("%.3f"))
    ax.set_title(r"$C_\ell$ vs $C_{\rm g}$ maxima", loc="left")
    ax.set_ylabel(r"$A_*^{(C_\ell)}/A_*^{(C_{\rm g})}$")
    save(fig, out)


def fig_lattice(rows, out):
    """Maxima of C_l counted on 256^3 Gaussian random fields over the Kac-Rice predictions."""
    cases = [("gaussian", GAUSS, "Gaussian (grid calibration)"), ("fNL=-2", DIVERGING["fNL=-2"], LABELS["fNL=-2"]),
             ("fNL=+2", DIVERGING["fNL=+2"], LABELS["fNL=+2"]), ("mu=5/6", CATEGORICAL[2], r"USR $\mu_*=5/6$")]
    fig, axes = plt.subplots(1, 2, figsize=(TEXTWIDTH, 2.5), sharey=True, layout="constrained")
    for ax, key, title in ((axes[0], "3b", "prediction: maxima of $C_\\ell$"),
                           (axes[1], "3a", "prediction: maxima of $C_{\\rm g}$")):
        for i, (case, c, label) in enumerate(cases):
            sel = [r for r in rows if r["case"] == case]
            nu = [float(r["nu"]) + 0.03 * (i - 1.5) for r in sel]
            val = [float(r[f"lattice_over_{key}"]) for r in sel]
            err = [float(r["err_rel"]) * v for r, v in zip(sel, val)]
            ax.errorbar(nu, val, yerr=err, fmt="-o", color=c, label=label, capsize=1.5, elinewidth=0.8,
                        markeredgecolor=SURFACE, markeredgewidth=0.6)
        ax.axhline(1.0, color=AXIS, linewidth=0.8)
        ax.grid(True, color=GRID, linewidth=0.5)
        ax.set_axisbelow(True)
        for side in ("top", "right"):
            ax.spines[side].set_visible(False)
        ax.set_xticks([1.5, 2.0, 2.5, 3.0])
        ax.set_xlabel(r"threshold $C_\ell > \nu\,\sigma_c$")
        ax.set_title(title, loc="left")
    axes[0].set_ylabel("lattice count / prediction")
    axes[0].set_ylim(0.88, 1.22)
    axes[0].legend(loc="upper right", fontsize=6.4)
    save(fig, out)


def fig_conventions(conv_rows, out):
    fig, ax = plt.subplots(figsize=(4.6, 4.1))
    conv = [
        ("R_pbhstat", "pbhstat v0.0.6 code"),
        ("R_ps2", r"PS $\times 2$, BBKS, $V=4\pi r^3/3$"),
        ("R_bbks", r"PS $\times 1$, BBKS, $V=4\pi r^3/3$"),
        ("R_gaussV", r"PS $\times 1$, BBKS, $V=(2\pi)^{3/2}r^3$"),
    ]
    for (key, label), c in zip(conv, CATEGORICAL):
        w = [float(r["width"]) for r in conv_rows]
        v = [float(r[key]) for r in conv_rows]
        ax.plot(w, v, "-o", color=c, label=label, markeredgecolor=SURFACE, markeredgewidth=0.8)
    style(ax)
    ax.set_ylim(0.72, 0.94)
    ax.set_title("Gaussian R: normalisation conventions only", loc="left")
    ax.text(0.02, 0.02, "pbhstat physics setup (linear estimator, no transfer)", transform=ax.transAxes,
            fontsize=7, color=MUTED)
    ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.17), ncol=1, fontsize=7)
    fig.tight_layout()
    save(fig, out, paper=False)


def main():
    rows = read("results/phase2_amplitudes.csv")
    conv = read("results/phase1_convention_spread.csv")
    Path("results").mkdir(exist_ok=True)
    fig_amplitudes(rows, "results/fig_phase2_threshold_amplitudes.png")
    fig_ratio(rows, "results/fig_phase3_ratio.png")
    fig_conventions(conv, "results/fig_phase1_conventions.png")
    if Path("results/phase3b_lattice.csv").exists():
        fig_lattice(read("results/phase3b_lattice.csv"), "results/fig_phase3b_lattice.png")


if __name__ == "__main__":
    main()
