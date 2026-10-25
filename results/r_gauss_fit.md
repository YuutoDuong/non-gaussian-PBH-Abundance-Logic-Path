# R_Gauss: table, fit and the factorisation for the BPL

R = A*_peaks/A*_threshold, peaks of C_l, f_PBH = 1, k* = 1e7/Mpc, Musco threshold with the transfer function, k < 4.49/r, K = 3.3, gamma = 0.36, BBKS normalisation with V = 4 pi r^3/3.

| spectrum | R_Gauss | fit | A*_threshold |
|---|---|---|---|
| lognormal, Delta = 0.1 | 0.7698 | 0.7682 | 1.9159e-02 |
| lognormal, Delta = 0.15 | 0.7708 | 0.7709 | 1.9801e-02 |
| lognormal, Delta = 0.2 | 0.7721 | 0.7727 | 2.0584e-02 |
| lognormal, Delta = 0.3 | 0.7746 | 0.7754 | 2.2347e-02 |
| lognormal, Delta = 0.4 | 0.7767 | 0.7773 | 2.4248e-02 |
| lognormal, Delta = 0.5 | 0.7784 | 0.7788 | 2.6293e-02 |
| lognormal, Delta = 0.7 | 0.7809 | 0.7810 | 3.0869e-02 |
| lognormal, Delta = 1 | 0.7834 | 0.7833 | 3.8696e-02 |
| lognormal, Delta = 1.4 | 0.7859 | 0.7855 | 5.0090e-02 |
| lognormal, Delta = 2 | 0.7886 | 0.7879 | 6.7850e-02 |
| BPL(4,2,1) | 0.7779 | | 2.3834e-02 |
| BPL(4,3,1) | 0.7761 | | 2.7174e-02 |
| BPL(4,4,1) | 0.7750 | | 3.0035e-02 |

Fit for the lognormal, 0.1 <= Delta <= 2: R_Gauss = 0.7833 + 0.0066 ln Delta, largest residual 0.0016.

Double ratio R/R_Gauss for BPL(4,3,1):

| template | R | R/R_Gauss |
|---|---|---|
| gaussian | 0.7761 | 1.0000 |
| fNL=-2 | 0.7536 | 0.9710 |
| fNL=-1 | 0.7539 | 0.9714 |
| fNL=-0.5 | 0.7619 | 0.9817 |
| fNL=+0.5 | 0.7875 | 1.0147 |
| fNL=+1 | 0.7945 | 1.0236 |
| fNL=+2 | 0.8018 | 1.0331 |
| USR mu=5/2 | 0.7869 | 1.0138 |
| curvaton r=0.9 | 0.8007 | 1.0316 |

Quadratic templates, |f_NL| <= 2: R/R_Gauss = 0.971 to 1.033.
