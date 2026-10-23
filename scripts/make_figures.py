"""Figures from results/*.csv (static PNGs for the repo and the note).

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

INK, INK2, MUTED, GRID, AXIS = "#0b0b0b", "#52514e", "#898781", "#e1e0d9", "#c3c2b7"
SURFACE = "#fcfcfb"
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
    "fNL=-2": r"$f_{NL}=-2$",
    "fNL=-1": r"$f_{NL}=-1$",
    "fNL=-0.5": r"$f_{NL}=-0.5$",
    "fNL=+0.5": r"$f_{NL}=+0.5$",
    "fNL=+1": r"$f_{NL}=+1$",
    "fNL=+2": r"$f_{NL}=+2$",
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
        "axes.titlesize": 10,
        "axes.labelsize": 9,
        "xtick.color": MUTED,
        "ytick.color": MUTED,
        "xtick.labelsize": 8,
        "ytick.labelsize": 8,
        "legend.fontsize": 7.5,
        "legend.frameon": False,
        "font.family": "sans-serif",
        "lines.linewidth": 1.6,
        "lines.markersize": 4.5,
    }
)


def read(path):
    with open(path) as fh:
        return list(csv.DictReader(fh))


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


def ratio(rows, template):
    w, a_pk = series(rows, template, "peaks")
    _, a_th = series(rows, template, "threshold")
    return w, [p / t for p, t in zip(a_pk, a_th)]


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
    fig.savefig(out, dpi=200)
    plt.close(fig)


def fig_ratio(rows, conv_rows, out):
    fig, axes = plt.subplots(1, 3, figsize=(12.5, 4.3))
    ax = axes[0]
    for t, c in DIVERGING.items():
        w, r = ratio(rows, t)
        ax.plot(w, r, "-o", color=c, label=LABELS[t], markeredgecolor=SURFACE, markeredgewidth=0.8)
    quad_labels = [(ratio(rows, t)[0][-1], ratio(rows, t)[1][-1], LABELS[t]) for t in ("fNL=-2", "gaussian", "fNL=+2")]
    ax.set_title("Quadratic template", loc="left")
    ax.set_ylabel(r"$R=A_*^{\rm peaks}/A_*^{\rm threshold}$")
    ax.legend(loc="upper left", ncol=2)
    ax = axes[1]
    w, r = ratio(rows, "gaussian")
    ax.plot(w, r, "-o", color=GAUSS, label="Gaussian", markeredgecolor=SURFACE, markeredgewidth=0.8)
    log_labels = [(ratio(rows, "gaussian")[0][-1], ratio(rows, "gaussian")[1][-1], "Gaussian")]
    for t, c in zip(LOG_FAMILY, CATEGORICAL):
        w, r = ratio(rows, t)
        ax.plot(w, r, "-o", color=c, label=LABELS[t], markeredgecolor=SURFACE, markeredgewidth=0.8)
        log_labels.append((w[-1], r[-1], LABELS[t]))
    ax.set_title("Logarithmic (USR) template", loc="left")
    ax.legend(loc="upper left")
    for a_ in axes[:2]:
        style(a_)
        a_.set_ylim(0.72, 0.94)
        a_.set_xlim(0.085, 3.4)
    end_labels(axes[0], quad_labels)
    end_labels(axes[1], log_labels)
    ax = axes[2]
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
    ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.17), ncol=2, fontsize=7)
    fig.suptitle(
        "Phase 3a preview: peaks in the Gaussian-peak approximation (not the final Phase 3 result); "
        "left two panels share every ingredient except the peak constraint",
        fontsize=9,
        color=INK2,
        x=0.01,
        ha="left",
    )
    fig.tight_layout()
    fig.savefig(out, dpi=200)
    plt.close(fig)


def main():
    rows = read("results/phase2_amplitudes.csv")
    conv = read("results/phase1_convention_spread.csv")
    Path("results").mkdir(exist_ok=True)
    fig_amplitudes(rows, "results/fig_phase2_threshold_amplitudes.png")
    fig_ratio(rows, conv, "results/fig_phase3a_ratio_preview.png")


if __name__ == "__main__":
    main()
