# The threshold of 2306.17149: their Fig. S3 against three readings of the prescription

delta_th from the Musco et al. prescription (their eqs. S4-S6), read from the vector paths of their Fig. S3 (lognormal: delta_th vs Delta; BPL(4, beta, 1): delta_th vs beta), against ngpbh.threshold.MuscoThreshold without the transfer function, with T at eta = r_m/10 (Musco et al.'s superhorizon smoothing, 2011.03014 Sec. IV) and with T at horizon crossing (eta = r_m, the transfer function of their variances, which their notation T(k, r_m) suggests). Digitised curves: results/digitized_2306.17149_figS3.csv; figure: results/fig_threshold_check.png.

## lognormal

| Delta | 2306.17149 | no transfer function | T at eta = r_m/10 (Musco et al.) | T at horizon crossing, eta = r_m |
|---|---|---|---|---|
| 0.1 | 0.5925 | 0.5926 | 0.5926 | 0.5926 |
| 0.15 | 0.5897 | 0.5897 | 0.5897 | 0.5900 |
| 0.2 | 0.5854 | 0.5854 | 0.5854 | 0.5865 |
| 0.25 | 0.5795 | 0.5795 | 0.5795 | 0.5822 |
| 0.3 | 0.5717 | 0.5717 | 0.5718 | 0.5775 |
| 0.35 | 0.5624 | 0.5624 | 0.5625 | 0.5726 |
| 0.4 | 0.5519 | 0.5520 | 0.5521 | 0.5679 |
| 0.45 | 0.5410 | 0.5410 | 0.5412 | 0.5634 |
| 0.5 | 0.5301 | 0.5301 | 0.5304 | 0.5593 |
| 0.55 | 0.5196 | 0.5196 | 0.5199 | 0.5556 |
| 0.6 | 0.5098 | 0.5098 | 0.5102 | 0.5524 |
| 0.65 | 0.5007 | 0.5007 | 0.5012 | 0.5495 |
| 0.7 | 0.4925 | 0.4925 | 0.4930 | 0.5469 |
| 0.75 | 0.4851 | 0.4851 | 0.4856 | 0.5447 |
| 0.8 | 0.4784 | 0.4784 | 0.4789 | 0.5427 |
| 0.85 | 0.4724 | 0.4724 | 0.4728 | 0.5410 |
| 0.9 | 0.4669 | 0.4669 | 0.4674 | 0.5396 |
| 0.95 | 0.4621 | 0.4621 | 0.4626 | 0.5383 |
| 1 | 0.4579 | 0.4577 | 0.4581 | 0.5371 |

## BPL(4, beta, 1)

| beta | 2306.17149 | no transfer function | T at eta = r_m/10 (Musco et al.) | T at horizon crossing, eta = r_m |
|---|---|---|---|---|
| 1.5 | 0.4986 | 0.5054 | 0.5059 | 0.5592 |
| 1.6 | 0.5073 | 0.5090 | 0.5095 | 0.5602 |
| 1.7 | 0.5140 | 0.5124 | 0.5129 | 0.5611 |
| 1.8 | 0.5187 | 0.5156 | 0.5161 | 0.5619 |
| 1.9 | 0.5203 | 0.5187 | 0.5191 | 0.5628 |
| 2 | 0.5220 | 0.5216 | 0.5220 | 0.5636 |
| 2.1 | 0.5241 | 0.5243 | 0.5247 | 0.5643 |
| 2.2 | 0.5264 | 0.5269 | 0.5273 | 0.5650 |
| 2.3 | 0.5294 | 0.5294 | 0.5298 | 0.5658 |
| 2.4 | 0.5320 | 0.5318 | 0.5321 | 0.5664 |
| 2.5 | 0.5341 | 0.5340 | 0.5344 | 0.5671 |
| 2.6 | 0.5361 | 0.5362 | 0.5365 | 0.5677 |
| 2.7 | 0.5382 | 0.5382 | 0.5385 | 0.5683 |
| 2.8 | 0.5402 | 0.5402 | 0.5405 | 0.5689 |
| 2.9 | 0.5420 | 0.5420 | 0.5423 | 0.5694 |
| 3 | 0.5438 | 0.5438 | 0.5441 | 0.5700 |

## Summary: largest |ngpbh - 2306.17149| in delta_th

| reading of the prescription | lognormal, 0.1 <= Delta <= 1 | BPL, 1.9 <= beta <= 3 | BPL, 1.5 <= beta <= 3 |
|---|---|---|---|
| no transfer function | 0.0002 | 0.0016 | 0.0068 |
| T at eta = r_m/10 (Musco et al.) | 0.0005 | 0.0012 | 0.0073 |
| T at horizon crossing, eta = r_m | 0.0792 | 0.0425 | 0.0606 |

Reading: their thresholds are those of the prescription without the transfer function, or with Musco et al.'s superhorizon smoothing; the two differ by <= 0.0005 here. The horizon-crossing transfer function would raise delta_th by up to 0.08 over this range (0.11 at Delta = 2). For beta < 2 the integrals without T converge only conditionally (UV tail k^-beta); their curve has a kink at beta ~ 1.8 and lies up to 0.007 below ours there, presumably a numerical cut-off on their side.

Consistency of their two curves through eq. S6 (delta_th = delta_c(alpha_s)):

* Delta = 1: alpha_s = 0.512 gives delta_th = 0.4579, plotted 0.4579
* Delta = 0.5: alpha_s = 2.045 gives delta_th = 0.5301, plotted 0.5301
* beta = 2: alpha_s = 1.304 gives delta_th = 0.5037, plotted 0.5220
* beta = 3: alpha_s = 2.717 gives delta_th = 0.5472, plotted 0.5438

For the lognormal the plotted alpha_s and delta_th are consistent (and agree with ours); for the BPL the plotted alpha_s is not the one behind the plotted delta_th, which is what enters the abundance.
