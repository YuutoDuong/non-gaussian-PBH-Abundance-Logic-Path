"""Extract the f_PBH bands and NANOGrav15 contours of 2306.17149 Fig. 3 from the arXiv source.

The figure (Abundance.pdf) is vector graphics: the band edges are stroked polylines and the axis
ticks are line segments, so the curves are read exactly (no pixel digitising). Axis calibration:
major ticks at log10 k* = 6.0 ... 9.0 and log10 A = 0.0 ... -2.5. Writes
results/digitized_2306.17149_fig3.csv with columns curve, edge, log10_kstar, log10_A, where edge is
"f=1" / "f=1e-3" for the bands (upper / lower edge) and "1sigma" / "2sigma" / "3sigma" for the
NANOGrav15 posterior (marginalised over the BPL shape; from their Fig. 1).

Needs PyMuPDF (pip install pymupdf). Usage: python scripts/digitize_2306_17149.py [--src DIR]
"""

import argparse
import csv
import io
import tarfile
import urllib.request
from pathlib import Path

import numpy as np

URL = "https://arxiv.org/e-print/2306.17149v3"
NAMES = {(0.4, 0.8, 0.933): "fNL=-2", (0.267, 0.467, 0.667): "curvaton r=0.9", (0.0, 0.0, 0.0): "gaussian",
         (0.933, 0.4, 0.467): "USR beta=3", (0.133, 0.533, 0.2): "NG15"}


def calibration(drawings):
    """Pixel-to-data maps from the major tick marks (length 4 pt) on the left and bottom axes."""
    xs, ys = set(), set()
    for d in drawings:
        if d["type"] != "s" or d.get("color") != (0.0, 0.0, 0.0):
            continue
        for it in d["items"]:
            if it[0] != "l":
                continue
            p, q = it[1], it[2]
            if abs(p.x - q.x) < 0.01 and abs(abs(p.y - q.y) - 4.0) < 0.01:
                xs.add(round(p.x, 3))
            if abs(p.y - q.y) < 0.01 and abs(abs(p.x - q.x) - 4.0) < 0.01:
                ys.add(round(p.y, 3))
    xs, ys = np.array(sorted(xs)), np.array(sorted(ys))
    assert len(xs) == 7 and len(ys) == 6, (xs, ys)  # 6.0..9.0 by 0.5; 0.0..-2.5 by -0.5
    fx = np.polyfit(xs, np.linspace(6.0, 9.0, 7), 1)
    fy = np.polyfit(ys, np.linspace(0.0, -2.5, 6), 1)
    return lambda x: np.polyval(fx, x), lambda y: np.polyval(fy, y)


def main():
    import pymupdf

    ap = argparse.ArgumentParser()
    ap.add_argument("--src", default="data/arxiv_2306.17149")
    ap.add_argument("--out", default="results/digitized_2306.17149_fig3.csv")
    args = ap.parse_args()
    src = Path(args.src)
    pdf = src / "Abundance.pdf"
    if not pdf.exists():
        src.mkdir(parents=True, exist_ok=True)
        req = urllib.request.Request(URL, headers={"User-Agent": "ngpbh"})
        with urllib.request.urlopen(req) as r, tarfile.open(fileobj=io.BytesIO(r.read()), mode="r:gz") as tar:
            member = tar.getmember("Abundance.pdf")
            pdf.write_bytes(tar.extractfile(member).read())
    drawings = pymupdf.open(pdf)[0].get_drawings()
    to_x, to_y = calibration(drawings)
    curves = {}
    for d in drawings:
        col = tuple(round(c, 3) for c in d["color"]) if d.get("color") else None
        width = round(d.get("width") or 0, 2)
        if d["type"] != "s" or col not in NAMES or width not in (1.15, 2.16):
            continue
        pts = [q for it in d["items"] for q in ((it[1], it[2]) if it[0] == "l" else (it[1], it[4]) if it[0] == "c" else ())]
        xy = np.array([(to_x(p.x), to_y(p.y)) for p in pts])
        if len(xy) < 50:  # legend swatches and the dashed <m> = Msun line
            continue
        curves.setdefault(NAMES[col], []).append(xy)
    rows = []
    for name, polys in curves.items():
        polys = sorted(polys, key=lambda xy: -np.mean(xy[:, 1]) if name != "NG15" else np.ptp(xy[:, 0]))
        labels = ["f=1", "f=1e-3"] if name != "NG15" else ["1sigma", "2sigma", "3sigma"]
        for label, xy in zip(labels, polys):
            if name != "NG15":
                xy = xy[np.argsort(xy[:, 0])]
            rows += [dict(curve=name, edge=label, log10_kstar=round(x, 4), log10_A=round(y, 4)) for x, y in xy]
    with open(args.out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=["curve", "edge", "log10_kstar", "log10_A"])
        w.writeheader()
        w.writerows(rows)
    print(f"wrote {len(rows)} points for {sorted(curves)} to {args.out}")


if __name__ == "__main__":
    main()
