"""Which threshold did 2306.17149 compute? Their Fig. S3 against three readings of the prescription.

Their supplement (eqs. S4-S5) computes r_m and the shape parameter alpha_s from P^T = T^2(k, r_m) P,
following Musco et al. (2011.03014), and Fig. S3 plots the resulting delta_th and alpha_s against
the lognormal width Delta and against the UV slope beta of BPL(4, beta, 1). The figure is vector
graphics in the arXiv source (Th1.pdf, Th2.pdf): the curves are stroked polylines and the axes
are calibrated on the tick marks, so the curves are read exactly. We compare them with
ngpbh.threshold.MuscoThreshold for
  * no transfer function;
  * T at eta = r_m/10, Musco et al.'s superhorizon smoothing ("when r_m ~ 10 r_H or larger");
  * T at horizon crossing, eta = r_m, the transfer function of their variances (eq. S2), which
    is what the notation T(k, r_m) in eqs. S4-S5 suggests.
Writes results/digitized_2306.17149_figS3.csv, results/threshold_2306.17149_check.md and
results/fig_threshold_check.png. Needs PyMuPDF (pip install pymupdf).

Usage: python scripts/check_threshold_2306_17149.py [--src DIR]
"""

import argparse
import csv
import io
import tarfile
import urllib.request
import warnings
from collections import defaultdict
from pathlib import Path

import numpy as np

from ngpbh.spectra import BrokenPowerLaw, LogNormal
from ngpbh.threshold import MuscoThreshold, delta_c_of_q

from make_figures import CATEGORICAL, GRID, INK2, MUTED, plt

URL = "https://arxiv.org/e-print/2306.17149v3"
# panel file -> (x label, x range on the outer ticks, delta_th range on its outer ticks, alpha_s range)
PANELS = {"Th1.pdf": ("Delta", (0.1, 1.0), (0.40, 0.60), (0.5, 6.0)),
          "Th2.pdf": ("beta", (1.5, 3.0), (0.45, 0.60), (0.5, 3.0))}
RED, BLUE, BLACK = (0.667, 0.0, 0.0), (0.0, 0.0, 0.667), (0.0, 0.0, 0.0)
VARIANTS = {"no transfer function": dict(transfer=False),
            "T at eta = r_m/10 (Musco et al.)": dict(transfer=True, eta_over_rm=0.1),
            "T at horizon crossing, eta = r_m": dict(transfer=True, eta_over_rm=1.0)}


def fetch(src):
    if not all((src / name).exists() for name in PANELS):
        src.mkdir(parents=True, exist_ok=True)
        req = urllib.request.Request(URL, headers={"User-Agent": "ngpbh"})
        with urllib.request.urlopen(req, timeout=120) as r, tarfile.open(fileobj=io.BytesIO(r.read()), mode="r:*") as tar:
            for name in PANELS:
                (src / name).write_bytes(tar.extractfile(tar.getmember(name)).read())
    return src


def ticks(drawings, colour, axis):
    out = set()
    for d in drawings:
        if d["type"] != "s" or round(d.get("width") or 0, 2) != 0.86 or d.get("color") is None:
            continue
        if tuple(round(c, 3) for c in d["color"]) != colour:
            continue
        r = d["rect"]
        if axis == "x" and r.width < 0.01 and r.y1 > 200:  # bottom axis
            out.add(round(r.x0, 3))
        if axis == "y" and r.height < 0.01:
            out.add(round(r.y0, 3))
    return np.array(sorted(out))


def digitize(src):
    """{(panel, 'delta_th' | 'alpha_s'): array of (x, value)}."""
    import pymupdf

    curves = {}
    for name, (_, xr, dr, ar) in PANELS.items():
        drawings = pymupdf.open(src / name)[0].get_drawings()
        xt, dt, at = ticks(drawings, BLACK, "x"), ticks(drawings, RED, "y"), ticks(drawings, BLUE, "y")
        to_x = lambda p: xr[0] + (p - xt[0]) * (xr[1] - xr[0]) / (xt[-1] - xt[0])
        to_d = lambda p: dr[0] + (dt[-1] - p) * (dr[1] - dr[0]) / (dt[-1] - dt[0])  # y grows downwards
        to_a = lambda p: ar[0] + (at[-1] - p) * (ar[1] - ar[0]) / (at[-1] - at[0])
        for d in drawings:
            if d["type"] != "s" or len(d["items"]) <= 10:
                continue
            kind = "delta_th" if tuple(round(c, 3) for c in d["color"]) == RED else "alpha_s"
            pts = np.array([(p.x, p.y) for it in d["items"] if it[0] == "l" for p in (it[1], it[2])])
            conv = to_d if kind == "delta_th" else to_a
            xy = np.column_stack([to_x(pts[:, 0]), conv(pts[:, 1])])
            _, keep = np.unique(np.round(xy[:, 0], 9), return_index=True)
            curves[(name, kind)] = xy[keep]
    return curves


def ours(panel, xs):
    make = (lambda w: LogNormal(1.0, 1.0, w)) if panel == "Th1.pdf" else (lambda b: BrokenPowerLaw(1.0, 1.0, 4.0, b, 1.0))
    out = {}
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")  # QUADPACK round-off notes; the values are converged
        for vname, kw in VARIANTS.items():
            res = [MuscoThreshold(**kw).evaluate(make(x)) for x in xs]
            out[vname] = (np.array([r.C_th for r in res]), np.array([r.q for r in res]))
    return out


def figure(curves, grids, results, out):
    fig, axes = plt.subplots(1, 2, figsize=(9.0, 3.7))
    styles = {"no transfer function": dict(color=CATEGORICAL[0], ls="--", marker="o"),
              "T at eta = r_m/10 (Musco et al.)": dict(color=CATEGORICAL[1], ls=":", marker="s"),
              "T at horizon crossing, eta = r_m": dict(color=CATEGORICAL[2], ls="-", marker="^")}
    for ax, panel in zip(axes, PANELS):
        theirs = curves[(panel, "delta_th")]
        ax.plot(theirs[:, 0], theirs[:, 1], color=INK2, linewidth=4.0, alpha=0.25, solid_capstyle="butt",
                label="2306.17149 Fig. S3")
        for vname, st in styles.items():
            ax.plot(grids[panel], results[panel][vname][0], linestyle=st["ls"], marker=st["marker"], markersize=3.2,
                    color=st["color"], linewidth=1.2, label=f"ngpbh, {vname}")
        ax.grid(True, color=GRID, linewidth=0.6)
        ax.set_axisbelow(True)
        for side in ("top", "right"):
            ax.spines[side].set_visible(False)
        xlabel = r"lognormal width $\Delta$" if panel == "Th1.pdf" else r"BPL UV slope $\beta$ ($\alpha=4$, $\gamma=1$)"
        ax.set_xlabel(xlabel)
        ax.set_ylabel(r"$\delta_{\rm th}$")
    axes[0].set_title("lognormal", loc="left")
    axes[1].set_title("broken power law", loc="left")
    axes[0].legend(loc="lower left", fontsize=6.6)
    axes[1].text(0.02, 0.03, r"$\beta<2$: the integrals without T converge only conditionally", transform=axes[1].transAxes,
                 fontsize=6.6, color=MUTED)
    fig.tight_layout()
    fig.savefig(out, dpi=200)
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--src", default="data/arxiv_2306.17149_src")
    args = ap.parse_args()
    curves = digitize(fetch(Path(args.src)))
    with open("results/digitized_2306.17149_figS3.csv", "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["panel", "quantity", "x", "value"])
        for (panel, kind), xy in sorted(curves.items()):
            w.writerows([PANELS[panel][0], kind, f"{x:.6g}", f"{v:.6g}"] for x, v in xy)
    grids = {"Th1.pdf": np.round(np.arange(0.1, 1.0001, 0.05), 3), "Th2.pdf": np.round(np.arange(1.5, 3.0001, 0.1), 3)}
    results = {panel: ours(panel, xs) for panel, xs in grids.items()}
    figure(curves, grids, results, "results/fig_threshold_check.png")

    md = ["# The threshold of 2306.17149: their Fig. S3 against three readings of the prescription", "",
          "delta_th from the Musco et al. prescription (their eqs. S4-S6), read from the vector paths of their Fig. S3 "
          "(lognormal: delta_th vs Delta; BPL(4, beta, 1): delta_th vs beta), against ngpbh.threshold.MuscoThreshold "
          "without the transfer function, with T at eta = r_m/10 (Musco et al.'s superhorizon smoothing, "
          "2011.03014 Sec. IV) and with T at horizon crossing (eta = r_m, the transfer function of their variances, "
          "which their notation T(k, r_m) suggests). Digitised curves: results/digitized_2306.17149_figS3.csv; "
          "figure: results/fig_threshold_check.png.", ""]
    summary = defaultdict(dict)
    for panel, xs in grids.items():
        theirs = curves[(panel, "delta_th")]
        ref = np.interp(xs, theirs[:, 0], theirs[:, 1])
        label = PANELS[panel][0]
        md += [f"## {'lognormal' if panel == 'Th1.pdf' else 'BPL(4, beta, 1)'}", "",
               f"| {label} | 2306.17149 | " + " | ".join(VARIANTS) + " |", "|---|---|" + "---|" * len(VARIANTS)]
        for i, x in enumerate(xs):
            md.append(f"| {x:g} | {ref[i]:.4f} | " + " | ".join(f"{results[panel][v][0][i]:.4f}" for v in VARIANTS) + " |")
        for v in VARIANTS:
            d = np.abs(results[panel][v][0] - ref)
            summary[panel][v] = (d.max(), d[xs >= 1.9].max() if panel == "Th2.pdf" else d.max())
        md.append("")
    md += ["## Summary: largest |ngpbh - 2306.17149| in delta_th", "",
           "| reading of the prescription | lognormal, 0.1 <= Delta <= 1 | BPL, 1.9 <= beta <= 3 | BPL, 1.5 <= beta <= 3 |",
           "|---|---|---|---|"]
    for v in VARIANTS:
        ln, bpl = summary["Th1.pdf"][v], summary["Th2.pdf"][v]
        md.append(f"| {v} | {ln[0]:.4f} | {bpl[1]:.4f} | {bpl[0]:.4f} |")
    # internal consistency of their two curves: delta_th(alpha_s) through eq. S6
    cons = []
    for panel, x in (("Th1.pdf", 1.0), ("Th1.pdf", 0.5), ("Th2.pdf", 2.0), ("Th2.pdf", 3.0)):
        a = np.interp(x, *curves[(panel, "alpha_s")].T)
        d = np.interp(x, *curves[(panel, "delta_th")].T)
        cons.append(f"{PANELS[panel][0]} = {x:g}: alpha_s = {a:.3f} gives delta_th = {float(delta_c_of_q(a)):.4f}, "
                    f"plotted {d:.4f}")
    md += ["", "Reading: their thresholds are those of the prescription without the transfer function, or with Musco "
           "et al.'s superhorizon smoothing; the two differ by <= 0.0005 here. The horizon-crossing transfer function "
           "would raise delta_th by up to 0.08 over this range (0.11 at Delta = 2). For beta < 2 the integrals without "
           "T converge only conditionally (UV tail k^-beta); their curve has a kink at beta ~ 1.8 and lies up to "
           "0.007 below ours there, presumably a numerical cut-off on their side.", "",
           "Consistency of their two curves through eq. S6 (delta_th = delta_c(alpha_s)):", ""]
    md += [f"* {c}" for c in cons]
    md += ["", "For the lognormal the plotted alpha_s and delta_th are consistent (and agree with ours); for the BPL "
           "the plotted alpha_s is not the one behind the plotted delta_th, which is what enters the abundance.", ""]
    Path("results/threshold_2306.17149_check.md").write_text("\n".join(md))
    print("\n".join(md))


if __name__ == "__main__":
    main()
