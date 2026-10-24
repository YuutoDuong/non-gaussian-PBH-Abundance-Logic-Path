# Validation of the NG corrections to the SIGW spectrum

## 1. Monte Carlo engine against quadrature

Terms with both: the Gaussian one and the disconnected NG ones (hybrid = 2 Omega_bil[s, s2], reducible = Omega_G[s2], s2 the one-loop spectrum). Each term per unit A and F, k* = 1; MC: 16 randomised Sobol sequences of 2^14 points per strategy; pull = (MC - quadrature)/standard error.

| shape | term | x | quadrature | MC | MC/quad | pull |
|---|---|---|---|---|---|---|
| BPL(4,3,1) | gaussian | 0.001 | 8.67448e-08 | 8.67444e-08 ± 5.1e-13 | 1.0000 | -0.8 |
| BPL(4,3,1) | gaussian | 0.03 | 5.44644e-04 | 5.44642e-04 ± 3.2e-09 | 1.0000 | -0.8 |
| BPL(4,3,1) | gaussian | 0.3 | 6.54789e-02 | 6.54799e-02 ± 1.8e-06 | 1.0000 | +0.6 |
| BPL(4,3,1) | gaussian | 1 | 3.22221e-01 | 3.23206e-01 ± 6.3e-04 | 1.0031 | +1.6 |
| BPL(4,3,1) | gaussian | 2 | 6.82676e-02 | 6.85178e-02 ± 2.5e-04 | 1.0037 | +1.0 |
| BPL(4,3,1) | gaussian | 5 | 3.79977e-04 | 3.78733e-04 ± 4.6e-06 | 0.9967 | -0.3 |
| BPL(4,3,1) | hybrid | 0.001 | 1.58073e-07 | 1.58045e-07 ± 4.8e-11 | 0.9998 | -0.6 |
| BPL(4,3,1) | hybrid | 0.03 | 1.00966e-03 | 1.00945e-03 ± 3.0e-07 | 0.9998 | -0.7 |
| BPL(4,3,1) | hybrid | 0.3 | 1.34038e-01 | 1.34004e-01 ± 4.7e-05 | 0.9997 | -0.7 |
| BPL(4,3,1) | hybrid | 1 | 5.75068e-01 | 5.73793e-01 ± 1.9e-03 | 0.9978 | -0.7 |
| BPL(4,3,1) | hybrid | 2 | 5.60557e-01 | 5.58809e-01 ± 3.9e-03 | 0.9969 | -0.4 |
| BPL(4,3,1) | hybrid | 5 | 4.10236e-03 | 4.14153e-03 ± 4.7e-05 | 1.0095 | +0.8 |
| BPL(4,3,1) | reducible | 0.001 | 9.61724e-08 | 9.62629e-08 ± 8.6e-11 | 1.0009 | +1.1 |
| BPL(4,3,1) | reducible | 0.03 | 6.39893e-04 | 6.40375e-04 ± 5.2e-07 | 1.0008 | +0.9 |
| BPL(4,3,1) | reducible | 0.3 | 9.91727e-02 | 9.91674e-02 ± 7.1e-05 | 0.9999 | -0.1 |
| BPL(4,3,1) | reducible | 1 | 3.58164e-01 | 3.57749e-01 ± 1.7e-03 | 0.9988 | -0.2 |
| BPL(4,3,1) | reducible | 2 | 9.25145e-01 | 9.23128e-01 ± 7.7e-03 | 0.9978 | -0.3 |
| BPL(4,3,1) | reducible | 5 | 1.09693e-02 | 1.10501e-02 ± 9.8e-05 | 1.0074 | +0.8 |
| LN D=0.5 | gaussian | 0.001 | 1.11231e-07 | 1.11231e-07 ± 5.0e-13 | 1.0000 | +0.5 |
| LN D=0.5 | gaussian | 0.03 | 6.64315e-04 | 6.64316e-04 ± 2.8e-09 | 1.0000 | +0.5 |
| LN D=0.5 | gaussian | 0.3 | 6.51754e-02 | 6.51762e-02 ± 2.5e-06 | 1.0000 | +0.3 |
| LN D=0.5 | gaussian | 1 | 2.99639e-01 | 2.99470e-01 ± 9.0e-04 | 0.9994 | -0.2 |
| LN D=0.5 | gaussian | 2 | 9.82726e-02 | 9.84249e-02 ± 3.0e-04 | 1.0015 | +0.5 |
| LN D=0.5 | gaussian | 5 | 1.10084e-04 | 1.10536e-04 ± 9.8e-07 | 1.0041 | +0.5 |
| LN D=0.5 | hybrid | 0.001 | 2.26818e-07 | 2.26817e-07 ± 7.5e-11 | 1.0000 | -0.0 |
| LN D=0.5 | hybrid | 0.03 | 1.40923e-03 | 1.40918e-03 ± 4.2e-07 | 1.0000 | -0.1 |
| LN D=0.5 | hybrid | 0.3 | 1.69761e-01 | 1.69764e-01 ± 6.2e-05 | 1.0000 | +0.0 |
| LN D=0.5 | hybrid | 1 | 7.57266e-01 | 7.60428e-01 ± 2.1e-03 | 1.0042 | +1.5 |
| LN D=0.5 | hybrid | 2 | 7.47198e-01 | 7.38554e-01 ± 3.2e-03 | 0.9884 | -2.7 |
| LN D=0.5 | hybrid | 5 | 1.77890e-03 | 1.82902e-03 ± 3.0e-05 | 1.0282 | +1.7 |
| LN D=0.5 | reducible | 0.001 | 1.49992e-07 | 1.50094e-07 ± 1.1e-10 | 1.0007 | +0.9 |
| LN D=0.5 | reducible | 0.03 | 9.84231e-04 | 9.84744e-04 ± 6.6e-07 | 1.0005 | +0.8 |
| LN D=0.5 | reducible | 0.3 | 1.46335e-01 | 1.46303e-01 ± 8.8e-05 | 0.9998 | -0.4 |
| LN D=0.5 | reducible | 1 | 5.59474e-01 | 5.60695e-01 ± 2.9e-03 | 1.0022 | +0.4 |
| LN D=0.5 | reducible | 2 | 1.24545e+00 | 1.23873e+00 ± 4.6e-03 | 0.9946 | -1.5 |
| LN D=0.5 | reducible | 5 | 7.11269e-03 | 7.19381e-03 ± 1.4e-04 | 1.0114 | +0.6 |
| LN D=1 | gaussian | 0.001 | 2.04781e-07 | 2.04784e-07 ± 3.1e-12 | 1.0000 | +0.9 |
| LN D=1 | gaussian | 0.03 | 8.57559e-04 | 8.57562e-04 ± 9.6e-09 | 1.0000 | +0.3 |
| LN D=1 | gaussian | 0.3 | 4.03050e-02 | 4.03478e-02 ± 9.6e-05 | 1.0011 | +0.4 |
| LN D=1 | gaussian | 1 | 1.03175e-01 | 1.03261e-01 ± 2.3e-04 | 1.0008 | +0.4 |
| LN D=1 | gaussian | 2 | 7.20964e-02 | 7.22031e-02 ± 2.2e-04 | 1.0015 | +0.5 |
| LN D=1 | gaussian | 5 | 1.15168e-02 | 1.15192e-02 ± 8.9e-05 | 1.0002 | +0.0 |
| LN D=1 | hybrid | 0.001 | 2.41720e-07 | 2.41663e-07 ± 1.6e-10 | 0.9998 | -0.3 |
| LN D=1 | hybrid | 0.03 | 1.23231e-03 | 1.23237e-03 ± 5.3e-07 | 1.0000 | +0.1 |
| LN D=1 | hybrid | 0.3 | 9.39596e-02 | 9.40222e-02 ± 1.8e-04 | 1.0007 | +0.4 |
| LN D=1 | hybrid | 1 | 3.68851e-01 | 3.67367e-01 ± 1.4e-03 | 0.9960 | -1.1 |
| LN D=1 | hybrid | 2 | 4.03507e-01 | 4.00592e-01 ± 1.5e-03 | 0.9928 | -2.0 |
| LN D=1 | hybrid | 5 | 9.32455e-02 | 9.30035e-02 ± 3.5e-04 | 0.9974 | -0.7 |
| LN D=1 | reducible | 0.001 | 1.16116e-07 | 1.16377e-07 ± 1.1e-10 | 1.0022 | +2.3 |
| LN D=1 | reducible | 0.03 | 6.84940e-04 | 6.85922e-04 ± 6.1e-07 | 1.0014 | +1.6 |
| LN D=1 | reducible | 0.3 | 7.67124e-02 | 7.67107e-02 ± 9.6e-05 | 1.0000 | -0.0 |
| LN D=1 | reducible | 1 | 3.48013e-01 | 3.50215e-01 ± 2.1e-03 | 1.0063 | +1.1 |
| LN D=1 | reducible | 2 | 5.57341e-01 | 5.60594e-01 ± 2.7e-03 | 1.0058 | +1.2 |
| LN D=1 | reducible | 5 | 1.88091e-01 | 1.88855e-01 ± 1.6e-03 | 1.0041 | +0.5 |
| LN D=2 | gaussian | 0.001 | 7.62505e-06 | 7.62067e-06 ± 3.4e-09 | 0.9994 | -1.3 |
| LN D=2 | gaussian | 0.03 | 2.84180e-03 | 2.83254e-03 ± 6.7e-06 | 0.9967 | -1.4 |
| LN D=2 | gaussian | 0.3 | 2.27943e-02 | 2.27965e-02 ± 9.6e-05 | 1.0001 | +0.0 |
| LN D=2 | gaussian | 1 | 3.00734e-02 | 3.02729e-02 ± 1.4e-04 | 1.0066 | +1.4 |
| LN D=2 | gaussian | 2 | 2.64440e-02 | 2.64162e-02 ± 1.4e-04 | 0.9989 | -0.2 |
| LN D=2 | gaussian | 5 | 1.58810e-02 | 1.58161e-02 ± 6.1e-05 | 0.9959 | -1.1 |
| LN D=2 | hybrid | 0.001 | 2.00701e-06 | 2.00404e-06 ± 1.2e-08 | 0.9985 | -0.2 |
| LN D=2 | hybrid | 0.03 | 2.71351e-03 | 2.70385e-03 ± 1.2e-05 | 0.9964 | -0.8 |
| LN D=2 | hybrid | 0.3 | 5.72175e-02 | 5.76277e-02 ± 4.4e-04 | 1.0072 | +0.9 |
| LN D=2 | hybrid | 1 | 1.21015e-01 | 1.20402e-01 ± 6.7e-04 | 0.9949 | -0.9 |
| LN D=2 | hybrid | 2 | 1.32924e-01 | 1.32722e-01 ± 6.6e-04 | 0.9985 | -0.3 |
| LN D=2 | hybrid | 5 | 9.95790e-02 | 1.00212e-01 ± 6.3e-04 | 1.0064 | +1.0 |
| LN D=2 | reducible | 0.001 | 3.63050e-07 | 3.72172e-07 ± 4.8e-09 | 1.0251 | +1.9 |
| LN D=2 | reducible | 0.03 | 9.99661e-04 | 9.98783e-04 ± 2.7e-06 | 0.9991 | -0.3 |
| LN D=2 | reducible | 0.3 | 4.04003e-02 | 4.07934e-02 ± 2.8e-04 | 1.0097 | +1.4 |
| LN D=2 | reducible | 1 | 1.24868e-01 | 1.23657e-01 ± 6.3e-04 | 0.9903 | -1.9 |
| LN D=2 | reducible | 2 | 1.68217e-01 | 1.67541e-01 ± 6.8e-04 | 0.9960 | -1.0 |
| LN D=2 | reducible | 5 | 1.56122e-01 | 1.55457e-01 ± 1.1e-03 | 0.9957 | -0.6 |

72 comparisons: rms pull 0.97, max |pull| 2.7.

## 2. Against Adshead, Lozanov & Weiner (2105.01659) Fig. 4

Lognormal peak, A = F_NL = 1, read from the figure's vector paths. Hybrid and reducible by quadrature, the connected terms by randomised QMC (8 x 2^15 points per strategy), at 22 points in 0.05 < k/k* < 3. Columns: the largest |ALW| value of the curve, the largest |ngpbh - ALW| as a fraction of it, and the same where the curve is smooth (excluding points within 3% in k/k* of a steep edge, where the digitised curve's interpolation dominates).

| term | sigma | max |ALW| | max |diff| / max |ALW| | smooth parts |
|---|---|---|---|---|
| hybrid | 1/5 | 2.165 | 0.007 | 0.002 |
| hybrid | 1/10 | 2.478 | 0.008 | 0.002 |
| reducible | 1/5 | 2.080 | 0.003 | 0.003 |
| reducible | 1/10 | 2.304 | 0.002 | 0.002 |
| C | 1/5 | 0.679 | 0.012 | 0.003 |
| C | 1/10 | 1.001 | 0.011 | 0.008 |
| Z | 1/5 | 0.137 | 0.027 | 0.027 |
| Z | 1/10 | 0.183 | 0.024 | 0.001 |
| planar | 1/5 | 1.185 | 0.004 | 0.004 |
| planar | 1/10 | 1.563 | 0.024 | 0.014 |
| nonplanar | 1/5 | 0.151 | 0.016 | 0.011 |
| nonplanar | 1/10 | 0.204 | 0.022 | 0.022 |

Figure: results/fig_sigw_ng_validation.png.
