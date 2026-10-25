"""Validation of the primordial-NG corrections to the SIGW spectrum (ngpbh.sigw_ng).

1. Internal: the Monte Carlo engine against quadrature for the terms that have both, the Gaussian
   one and the disconnected NG ones (hybrid, reducible), for the four Phase 5 shapes.
2. External: all six NG terms against Adshead, Lozanov & Weiner (2105.01659) Fig. 4 (lognormal
   spectrum, sigma = 1/5 and 1/10, A = F_NL = 1). The figure is vector graphics in the arXiv
   source (lognormal-vs-monochromatic-all-contributions-linear.pdf): the curves are stroked
   polylines, the axes are calibrated on the major tick marks and their labels, so the curves are
   read exactly (no pixel digitising).

Writes results/digitized_2105.01659_fig4.csv, results/sigw_ng_validation.md and
results/fig_sigw_ng_validation.png. Step 2 needs PyMuPDF (pip install pymupdf).

Usage: python scripts/validate_sigw_ng.py [--src DIR]
"""

import argparse
import csv
import io
import tarfile
import urllib.request
from pathlib import Path

import numpy as np

from ngpbh.sigw_ng import disconnected_terms, mc_term
from ngpbh.spectra import LogNormal

from make_figures import AXIS, CATEGORICAL, GRID, SURFACE, TEXTWIDTH, plt, save
from pta_tension import SHAPES

URL = "https://arxiv.org/e-print/2105.01659v3"
FIG4 = "lognormal-vs-monochromatic-all-contributions-linear.pdf"
PANELS = [["hybrid", "reducible"], ["C", "Z"], ["planar", "nonplanar"]]
SIGMA_COLOURS = {(0.68, 0.154, 0.12): "1/40", (0.534, 0.68, 0.12): "1/20", (0.12, 0.68, 0.378): "1/10",
                 (0.12, 0.31, 0.68): "1/5"}
CONNECTED = ("C", "Z", "planar", "nonplanar")


def fetch(src):
    pdf = src / FIG4
    if not pdf.exists():
        src.mkdir(parents=True, exist_ok=True)
        req = urllib.request.Request(URL, headers={"User-Agent": "ngpbh"})
        with urllib.request.urlopen(req, timeout=120) as r, tarfile.open(fileobj=io.BytesIO(r.read()), mode="r:*") as tar:
            pdf.write_bytes(tar.extractfile(tar.getmember(FIG4)).read())
    return pdf


def digitize(pdf):
    """{(term, sigma): array of (x, value)} for the lognormal curves and the monochromatic one."""
    import pymupdf

    page = pymupdf.open(pdf)[0]
    dr = page.get_drawings()
    frames = [d["rect"] for d in dr if d["type"] == "s" and d.get("color") == (0.0, 0.0, 0.0)
              and round(d.get("width") or 0, 2) == 1.1]
    xs = sorted({round(r.x0, 4) for r in frames})
    ys = sorted({round(r.y0, 4) for r in frames} | {round(r.y1, 4) for r in frames})
    cols = [(xs[0], xs[1]), (xs[2], xs[3])]
    rows = [(ys[i], ys[i + 1]) for i in range(3)]
    # x axis: major ticks (vertical, 3.5 pt) at 10^-1 and 10^0 of the left column
    vt = sorted({round(d["rect"].x0, 3) for d in dr if d["type"] == "fs" and abs(d["rect"].height - 3.5) < 0.01
                 and d["rect"].width < 0.01})
    x_m1, x_0 = vt[0], vt[1]

    def to_x(px, c):
        return 10.0 ** ((px - (cols[c][0] - cols[0][0]) - x_0) / (x_0 - x_m1))

    # y axes: horizontal major ticks matched to the numeric tick labels (left column: labels on the
    # left; right column: on the right). A label's tick sits ~8 pt below the top of its box.
    labels = []
    for w in page.get_text("words"):
        try:
            labels.append((float(w[4].replace("−", "-")), w[:4]))
        except ValueError:
            pass
    hticks = [d["rect"] for d in dr if d["type"] == "fs" and abs(d["rect"].width - 3.5) < 0.01 and d["rect"].height < 0.01]
    to_y = {}
    for r, (y0, y1) in enumerate(rows):
        for c, (x0, x1) in enumerate(cols):
            ticks = sorted({round(t.y0, 3) for t in hticks if x0 - 1 <= t.x0 <= x1 + 1 and y0 - 0.5 <= t.y0 <= y1 + 0.5})
            pairs = set()
            for v, bb in labels:
                side = (20 < bb[0] and bb[2] < x0 + 2) if c == 0 else (x1 - 2 < bb[0] < x1 + 10)
                yc = bb[1] + 8.0
                if side and y0 - 6 <= yc <= y1 + 6:
                    t = min(ticks, key=lambda t: abs(t - yc))
                    if abs(t - yc) < 3:
                        pairs.add((t, v))
            pairs = np.array(sorted(pairs))
            fit = np.polyfit(pairs[:, 0], pairs[:, 1], 1)
            assert len(pairs) >= 3 and np.max(np.abs(np.polyval(fit, pairs[:, 0]) - pairs[:, 1])) < 1e-4, pairs
            to_y[(r, c)] = fit
    curves = {}
    for d in dr:
        if d["type"] != "s" or d.get("color") is None:
            continue
        col = tuple(round(v, 3) for v in d["color"])
        sigma = next((s for cc, s in SIGMA_COLOURS.items() if np.allclose(cc, col, atol=2e-3)), None)
        if sigma is None and not (d.get("dashes") not in (None, "[] 0") and col == (0.0, 0.0, 0.0)):
            continue
        pts = np.array([(p.x, p.y) for it in d["items"] if it[0] == "l" for p in (it[1], it[2])])
        if len(pts) < 20:  # legend swatches
            continue
        cx, cy = pts[:, 0].mean(), pts[:, 1].mean()
        c = 0 if cx < cols[1][0] else 1
        r = next(i for i, (y0, y1) in enumerate(rows) if y0 <= cy <= y1)
        x, v = to_x(pts[:, 0], c), np.polyval(to_y[(r, c)], pts[:, 1])
        _, keep = np.unique(np.round(x, 9), return_index=True)
        curves[(PANELS[r][c], sigma or "mono")] = np.column_stack([x[keep], v[keep]])
    return curves


def internal_checks():
    xs = np.array([1e-3, 0.03, 0.3, 1.0, 2.0, 5.0])
    out = []
    for name in ("BPL(4,3,1)", "LN D=0.5", "LN D=1", "LN D=2"):
        det = disconnected_terms(SHAPES[name], xs)
        for term in ("gaussian", "hybrid", "reducible"):
            for i, x in enumerate(xs):
                m, e = mc_term(SHAPES[name], term, x, n_log2=14, n_rand=16)
                out.append((name, term, x, det[term][i], m, e))
    return out


def external(curves, xs, n_log2=15, n_rand=8):
    ours = {}
    for sigma, width in (("1/5", 0.2), ("1/10", 0.1)):
        sp = LogNormal(1.0, 1.0, width)
        det = disconnected_terms(sp, xs)
        for term in ("hybrid", "reducible"):
            ours[(term, sigma)] = np.column_stack([det[term], np.zeros_like(xs)])
        for term in CONNECTED:
            ours[(term, sigma)] = np.array([mc_term(sp, term, x, n_log2=n_log2, n_rand=n_rand) for x in xs])
    return ours


def figure(curves, ours, xs, out):
    fig, axes = plt.subplots(3, 2, figsize=(0.84 * TEXTWIDTH, 5.4), sharex=True, layout="constrained")
    for r, row in enumerate(PANELS):
        for c, term in enumerate(row):
            ax = axes[r, c]
            for k, sigma in enumerate(("1/5", "1/10")):
                col = CATEGORICAL[k]
                ref = curves[(term, sigma)]
                ax.plot(ref[:, 0], ref[:, 1], color=col, linewidth=1.4, alpha=0.55,
                        label=f"Adshead et al., $\\Delta={sigma}$" if (r, c) == (0, 0) else None)
                val, err = ours[(term, sigma)][:, 0], ours[(term, sigma)][:, 1]
                ax.errorbar(xs, val, yerr=err, fmt="o" if k == 0 else "s", markersize=2.8, color=col,
                            markeredgecolor=SURFACE, markeredgewidth=0.4, elinewidth=0.7, linestyle="none",
                            label=f"ngpbh, $\\Delta={sigma}$" if (r, c) == (0, 0) else None)
            ax.axhline(0.0, color=AXIS, linewidth=0.8, zorder=0)
            ax.set_xscale("log")
            ax.grid(True, color=GRID, linewidth=0.5)
            ax.set_axisbelow(True)
            for side in ("top", "right"):
                ax.spines[side].set_visible(False)
            power = "A^3 F^2" if term in ("hybrid", "C", "Z") else "A^4 F^4"
            ax.set_title(f"{term}: $\\Omega_{{GW}}/{power}$", loc="left")
            if r == 2:
                ax.set_xlabel("$k/k_*$")
    axes[0, 0].legend(loc="upper left", fontsize=6.0)
    save(fig, out)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--src", default="data/arxiv_2105.01659")
    args = ap.parse_args()
    curves = digitize(fetch(Path(args.src)))
    with open("results/digitized_2105.01659_fig4.csv", "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["term", "sigma", "x", "value"])
        for (term, sigma), xy in sorted(curves.items()):
            w.writerows([term, sigma, f"{x:.6g}", f"{v:.6g}"] for x, v in xy)
    xs = np.geomspace(0.05, 3.0, 22)
    ours = external(curves, xs)
    figure(curves, ours, xs, "results/fig_sigw_ng_validation.png")
    md = ["# Validation of the NG corrections to the SIGW spectrum", "",
          "## 1. Monte Carlo engine against quadrature", "",
          "Terms with both: the Gaussian one and the disconnected NG ones (hybrid = 2 Omega_bil[s, s2], reducible = "
          "Omega_G[s2], s2 the one-loop spectrum). Each term per unit A and F, k* = 1; MC: 16 randomised Sobol "
          "sequences of 2^14 points per strategy; pull = (MC - quadrature)/standard error.", "",
          "| shape | term | x | quadrature | MC | MC/quad | pull |", "|---|---|---|---|---|---|---|"]
    pulls = []
    for name, term, x, d, m, e in internal_checks():
        pulls.append((m - d) / e)
        md.append(f"| {name} | {term} | {x:g} | {d:.5e} | {m:.5e} ± {e:.1e} | {m / d:.4f} | {(m - d) / e:+.1f} |")
    pulls = np.array(pulls)
    rms, worst = np.sqrt(np.mean(pulls**2)), np.abs(pulls).max()
    md += ["", f"{len(pulls)} comparisons: rms pull {rms:.2f}, max |pull| {worst:.1f}.", "",
           "## 2. Against Adshead, Lozanov & Weiner (2105.01659) Fig. 4", "",
           "Lognormal peak, A = F_NL = 1, read from the figure's vector paths. Hybrid and reducible by quadrature, "
           "the connected terms by randomised QMC (8 x 2^15 points per strategy), at 22 points in 0.05 < k/k* < 3. "
           "Columns: the largest |ALW| value of the curve, the largest |ngpbh - ALW| as a fraction of it, and the "
           "same where the curve is smooth (excluding points within 3% in k/k* of a steep edge, where the "
           "digitised curve's interpolation dominates).", "",
           "| term | sigma | max |ALW| | max |diff| / max |ALW| | smooth parts |", "|---|---|---|---|---|"]
    for term in ("hybrid", "reducible") + CONNECTED:
        for sigma in ("1/5", "1/10"):
            ref = curves[(term, sigma)]
            val = ours[(term, sigma)][:, 0]
            alw = np.interp(xs, ref[:, 0], ref[:, 1])
            lo = np.interp(xs * 0.97, ref[:, 0], ref[:, 1])
            hi = np.interp(xs * 1.03, ref[:, 0], ref[:, 1])
            scale = np.abs(ref[:, 1]).max()
            smooth = np.abs(hi - lo) < 0.1 * scale
            diff = np.abs(val - alw) / scale
            md.append(f"| {term} | {sigma} | {scale:.3f} | {diff.max():.3f} | {diff[smooth].max():.3f} |")
    md += ["", "Figure: results/fig_sigw_ng_validation.png.", ""]
    Path("results/sigw_ng_validation.md").write_text("\n".join(md))
    print("\n".join(md))


if __name__ == "__main__":
    main()
