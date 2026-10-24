# Phase 5 validation

Likelihood: NANOGrav 15-yr HD free spectrum, 14 bins, KDE refit (ceffyl data, Zenodo 10495907).
Intervals are median and 2-sigma (2.3%-97.7%) unless stated.

## 1. Power law h_c = A (f/yr^-1)^((3-gamma)/2)

| fit | this refit (median, 90%) | NANOGrav 2306.16213 (median, 90%) |
|---|---|---|
| gamma = 13/3: A | 2.52 (+0.59/-0.61) x 1e-15 | 2.4 (+0.7/-0.6) x 1e-15 |
| gamma free (uniform 0-7): gamma | 3.07 (+0.63/-0.56) | 3.2 (+0.6/-0.6) |

## 2. Lognormal SIGW, Delta marginalised (uniform 0.1-4)

| parameter | this work | 2306.17149 Fig. 2, NANOGrav15 |
|---|---|---|
| log10 k* | 8.16 (+1.23/-0.83) | 8.14 (+0.77/-0.76) |
| log10 A | -0.57 (+0.54/-0.84) | -0.56 (+0.53/-0.74) |
| Delta | 0.90 (+1.80/-0.80) | 0.98 (+1.13/-0.74) |

## 3. Fixed shapes (log10 k* in [5, 10], log10 A in [-3, 0])

| shape | log10 k* | log10 A | best fit (log10 k*, log10 A) | max ln L |
|---|---|---|---|---|
| LN Delta=0.5 | 7.84 (+0.65/-0.71) | -0.82 (+0.76/-0.68) | (7.73, -0.96) | -103.61 |
| LN Delta=1 | 8.09 (+0.49/-0.76) | -0.57 (+0.54/-0.71) | (8.14, -0.54) | -103.72 |
| LN Delta=2 | 8.82 (+0.34/-0.86) | -0.29 (+0.27/-0.66) | (9.13, -0.02) | -104.14 |
| BPL (4,3,1) | 7.78 (+0.67/-0.63) | -0.84 (+0.78/-0.61) | (7.65, -1.01) | -103.55 |

2306.17149 Fig. 1 (BPL, alpha, beta, gamma marginalised), NANOGrav15: log10 k* = 7.92 (+0.89/-0.62), log10 A = -0.86 (+0.79/-0.65).
