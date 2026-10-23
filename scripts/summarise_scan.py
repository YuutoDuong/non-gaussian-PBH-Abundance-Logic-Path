"""Markdown tables from results/phase2_amplitudes.csv: R for Phase 3a and 3b peaks, the double
ratio R / R_Gauss, the 3b/3a shift, and the NG shift of the threshold-statistics amplitude.

Usage: python scripts/summarise_scan.py [csv]
"""

import csv
import sys

WIDTHS = (0.1, 0.3, 0.5, 1.0, 2.0)


def main():
    path = sys.argv[1] if len(sys.argv) > 1 else "results/phase2_amplitudes.csv"
    rows = list(csv.DictReader(open(path)))
    templates = list(dict.fromkeys(r["template"] for r in rows))
    a = {(r["template"], r["statistic"], float(r["width"]), r["x_cut"]): float(r["A_star"]) for r in rows}

    def get(t, s, w, cut="4.49"):
        return a.get((t, s, w, cut), float("nan"))

    def table(title, fn):
        print(f"\n**{title}**\n")
        print("| template | " + " | ".join(f"Δ={w:g}" for w in WIDTHS) + " |")
        print("|---|" + "---|" * len(WIDTHS))
        for t in templates:
            print(f"| {t} | " + " | ".join(f"{fn(t, w):.3f}" for w in WIDTHS) + " |")

    r = lambda t, w, s: get(t, s, w) / get(t, "threshold", w)
    table("R = A_peaks / A_threshold, Phase 3b (maxima of C_l)", lambda t, w: r(t, w, "peaks_3b"))
    table("R, Phase 3a (maxima of C_g)", lambda t, w: r(t, w, "peaks"))
    table("Double ratio R / R_Gauss, Phase 3b", lambda t, w: r(t, w, "peaks_3b") / r("gaussian", w, "peaks_3b"))
    table("A_3b / A_3a", lambda t, w: get(t, "peaks_3b", w) / get(t, "peaks", w))
    table("A_threshold(NG) / A_threshold(Gauss)", lambda t, w: get(t, "threshold", w) / get("gaussian", "threshold", w))


if __name__ == "__main__":
    main()
