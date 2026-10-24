# Phase 6 audit: robustness of the Phase 5 significances

sigma_P (posterior mass of the PBH-allowed region, one-sided Gaussian equivalent), threshold statistics -> peaks 3b, base variant, f_PBH = 1. One change at a time from the reference. Last two columns: range over the cases of the statistics shift Delta sigma_P = 3b - threshold, and the largest change of any sigma_P relative to the reference.

| variation | BPL(4,3,1), gaussian | BPL(4,3,1), fNL=-2 | BPL(4,3,1), curvaton r=0.9 | LN D=0.5, gaussian | LN D=0.5, fNL=-2 | LN D=0.5, curvaton r=0.9 | LN D=1, gaussian | LN D=1, fNL=-2 | LN D=1, curvaton r=0.9 | Δσ_P range | max |change| |
|---|---|---|---|---|---|---|---|---|---|---|---|
| reference (14 bins, smooth 0.05, k* in [5,10], A in [-3,0]) | 2.98 → 3.62 | 1.87 → 2.18 | 2.01 → 2.44 | 3.02 → 3.52 | 1.88 → 2.31 | 2.00 → 2.43 | 3.54 → 4.00 | 3.43 → 3.84 | 3.09 → 3.49 | +0.31 to +0.68 | 0.00 |
| prior log10 k* in [6, 9] | 2.98 → 3.62 | 1.87 → 2.18 | 2.01 → 2.44 | 3.02 → 3.52 | 1.88 → 2.31 | 2.00 → 2.43 | 3.54 → 4.00 | 3.43 → 3.84 | 3.09 → 3.49 | +0.31 to +0.68 | 0.00 |
| prior log10 A in [-4, 0] | 2.98 → 3.62 | 1.87 → 2.18 | 2.01 → 2.44 | 3.02 → 3.52 | 1.88 → 2.31 | 2.00 → 2.43 | 3.54 → 4.00 | 3.43 → 3.84 | 3.09 → 3.49 | +0.31 to +0.68 | 0.00 |
| 5 frequency bins | 1.89 → 2.26 | 1.04 → 1.39 | 1.22 → 1.60 | 1.75 → 2.07 | 1.12 → 1.43 | 1.22 → 1.49 | 2.20 → 2.73 | 2.07 → 2.55 | 1.76 → 2.15 | +0.27 to +0.72 | 1.45 |
| 30 frequency bins | 3.80 → 4.36 | 2.93 → 3.16 | 3.04 → 3.36 | 3.82 → 4.26 | 2.92 → 3.24 | 3.02 → 3.33 | 4.05 → 4.47 | 3.96 → 4.32 | 3.68 → 4.01 | +0.23 to +0.55 | 1.06 |
| KDE tails as distributed (no smoothing) | 3.05 → 3.70 | 1.93 → 2.23 | 2.07 → 2.49 | 3.10 → 3.58 | 1.92 → 2.35 | 2.04 → 2.47 | 3.62 → 4.08 | 3.51 → 3.91 | 3.15 → 3.58 | +0.30 to +0.72 | 0.09 |
| KDE smoothing 0.03 dex | 3.01 → 3.65 | 1.90 → 2.21 | 2.04 → 2.46 | 3.06 → 3.54 | 1.90 → 2.33 | 2.02 → 2.45 | 3.57 → 4.03 | 3.46 → 3.86 | 3.12 → 3.52 | +0.31 to +0.69 | 0.03 |
| KDE smoothing 0.10 dex | 2.87 → 3.51 | 1.73 → 2.06 | 1.88 → 2.34 | 2.91 → 3.45 | 1.76 → 2.22 | 1.88 → 2.34 | 3.46 → 3.93 | 3.34 → 3.77 | 2.98 → 3.41 | +0.34 to +0.69 | 0.14 |
| posterior grid 2x finer | 2.98 → 3.62 | 1.87 → 2.18 | 2.01 → 2.44 | 3.02 → 3.52 | 1.88 → 2.31 | 2.00 → 2.43 | 3.54 → 4.00 | 3.43 → 3.84 | 3.09 → 3.49 | +0.31 to +0.68 | 0.00 |
| peaks normalisation: peak volume (2 pi)^3/2 r^3 | 2.98 → 3.72 | 1.87 → 2.26 | 2.01 → 2.56 | 3.02 → 3.61 | 1.88 → 2.42 | 2.00 → 2.56 | 3.54 → 4.11 | 3.43 → 3.94 | 3.09 → 3.57 | +0.40 to +0.81 | 0.14 |
| peaks normalisation: PS factor 2 | 2.98 → 3.56 | 1.87 → 2.13 | 2.01 → 2.38 | 3.02 → 3.47 | 1.88 → 2.25 | 2.00 → 2.36 | 3.54 → 3.94 | 3.43 → 3.78 | 3.09 → 3.44 | +0.27 to +0.59 | 0.09 |
| peaks normalisation: pbhstat v0.0.6 code | 2.98 → 3.37 | 1.87 → 2.03 | 2.01 → 2.22 | 3.02 → 3.32 | 1.88 → 2.10 | 2.00 → 2.19 | 3.54 → 3.77 | 3.43 → 3.63 | 3.09 → 3.30 | +0.16 to +0.39 | 0.33 |
| peaks normalisation: high-peak formula (2406.20089), indicative | 2.98 → 4.07 | 1.87 → 2.69 | 2.01 → 3.09 | 3.02 → 3.93 | 1.88 → 2.93 | 2.00 → 3.05 | 3.54 → 4.59 | 3.43 → 4.39 | 3.09 → 3.91 | +0.82 to +1.29 | 0.65 |

Cases in the ranges: shapes BPL(4,3,1), LN D=0.5, LN D=1; templates gaussian, fNL=-2, fNL=-1, fNL=-0.5, fNL=+0.5, curvaton r=0.9, USR mu=5/2; only cases with sigma_P(3b) <= 4 in the reference enter the last two columns.

