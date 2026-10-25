# Phase 5: does the choice of statistics change the PTA tension?

**Question.** 2306.17149 (Franciolini, Iovino, Vaskonen, Veermäe) find that PBH overproduction
disfavours a Gaussian scalar-induced GW (SIGW) explanation of NANOGrav 15 "at over 2σ" and the
quasi-inflection-point (USR) case at 3σ, using threshold statistics. Their footnote 1 says that
peaks theory gives "slightly smaller amplitudes", so their conclusions are conservative. Phase 5
quantifies that statement and compares it with the effect of NG. It uses the true peaks of the
non-Gaussian compaction (Phase 3b) and the same pipeline for every statistic.

## Pipeline

1. **SIGW** (`ngpbh.sigw`). This is Ω_GW for Gaussian ζ in radiation domination, with the kernel of
   2306.17149 eqs. (6)–(7). The present-day spectrum uses the Standard-Model g*(T) and g*s(T)
   (`ngpbh.thermal`, the Saikawa–Shirai table that PTArcade ships). Checks:
   * a flat spectrum gives Ω_GW = 0.8222436 A² (Kohri & Terada);
   * the Dirac closed form is reproduced exactly;
   * the quadrature agrees with scipy `dblquad` to 10⁻⁸–10⁻¹²;
   * the infrared k³ tail matches its asymptotic form.
2. **Likelihood** (`ngpbh.pta`). This is the NANOGrav 15-yr Hellings–Downs free spectrum (14 bins)
   as KDE log-densities, the ceffyl refit used by PTArcade. The data come from Zenodo 10495907
   via `scripts/fetch_ng15.py`, with MD5 checked and no pickles. Checks
   (`results/phase5_validation.md`):
   * **γ = 13/3 amplitude:** 2.52 (+0.59/−0.61) × 10⁻¹⁵, against 2.4 (+0.7/−0.6) × 10⁻¹⁵
     published by NANOGrav.
   * **Free γ:** 3.07 (+0.63/−0.56), against 3.2 ± 0.6.
   * **Lognormal SIGW fit, Δ marginalised:** log10 A = −0.57 (+0.54/−0.84), against
     −0.56 (+0.53/−0.74) in 2306.17149 Fig. 2; log10 k* = 8.16 against 8.14; Δ = 0.90 against
     0.98. Our upper tails in k* and Δ are longer, consistent with a k* prior that ends at 10⁹ in
     their fit.
   * **BPL fit, shape marginalised:** log10 A = −0.84 (+0.79/−0.74) against −0.86 (+0.79/−0.65)
     (their Fig. 1), and log10 k* = 7.95 against 7.92. Their digitised 1σ/2σ/3σ contours lie at
     our 69/94/99% HPD levels (`results/phase5_compare.md`).
3. **Overproduction lines** (`AbundanceModel.amplitude_curve`). The variances depend only on k r,
   so moving the spectrum to another k* only re-weights the f_PBH integral over M_H. One set of
   β grids then gives A_*(k*) at all k*. It agrees with a direct solve at each k* to 10⁻⁷.
4. **Tension.** The NG15 posterior in (log10 k*, log10 A) uses uniform priors on [5, 10] × [−3, 0]
   (A ≤ 1, as in 2306.17149) and the same spectral shape as the line. For each line we compute:
   * **σ_HPD:** the highest-posterior-density level at which the allowed region A < A_*(k*) is
     first reached, quoted as a two-sided 1D Gaussian equivalent. This is 2306.17149's "excluded
     at n σ": their 1σ/2σ/3σ contours lie at our 69/94/99% HPD levels, i.e. they use 68.3/95.4/99.7%.
   * **σ_P:** the posterior mass P in that region, as a one-sided Gaussian equivalent. This is the
     posterior probability that PBHs are not overproduced.
   * the profile Δχ², which does not depend on the prior.

   σ_P is the primary measure. σ_HPD saturates where the posterior ridge is lumpy: for f_NL < 0,
   several lines first touch the same local bump at k* ≈ 10^7.05.

   **The KDE tails need regularising.** The free-spectrum KDEs use a compact kernel with
   bandwidths down to 0.004 dex. Outside the sampled range the density is floored at ln p = −36,
   and the sparse tails are combs of single-sample spikes. The two most constraining bins (3.9 and
   5.9 nHz) have almost no low-ρ tail. `FreeSpectrum(smooth=0.05)` replaces each density by
   max(KDE, KDE ⊗ 0.05-dex Gaussian):
   * the bulk is unchanged (the γ = 13/3 fit moves by 0.004 dex);
   * σ changes by ≤ 0.1 relative to the KDEs as distributed (`sigma_ceffyl` in
     `phase5_tension.csv`).

   Values above 4σ lie where the KDEs have no samples and are quoted only as ">4".
5. **NG corrections to the SIGW** (`ngpbh.sigw_ng`, stage `ngsigw`). With
   ζ = ζ_g + F(ζ_g² − ⟨ζ_g²⟩), F = 3f_NL/5, the spectrum is
   Ω = A² T2 + F²A³ (hybrid + C + Z) + F⁴A⁴ (reducible + planar + nonplanar). These are the seven
   terms of Adshead, Lozanov & Weiner (2105.01659); A is the amplitude of the Gaussian component,
   as in the PBH lines. Stages 1–4 use the Gaussian SIGW, as 2306.17149 do. Stage `ngsigw` gives
   every NG template its own posterior. The series is exact at A⁴ for the quadratic template.
   For the curvaton and USR maps it is truncated at A³, where only F and the cubic coefficient G
   enter (G through Ω += 12Gσ² A² T2). Method and conventions: `docs/conventions.md`.
   Checks (`results/sigw_ng_validation.md`):
   * the Monte Carlo engine reproduces the quadrature of the Gaussian, hybrid and reducible terms
     (72 comparisons, rms pull 0.97);
   * all six NG terms reproduce Adshead et al. Fig. 4 for lognormal peaks of width 1/5 and 1/10,
     read from the figure's vector paths, to ≤ 3% of each curve's maximum
     (`results/fig_sigw_ng_validation.png`).

## Exit check: 2306.17149 Fig. 3

The figure is vector graphics. `scripts/digitize_2306_17149.py` reads the band edges and the NG15
contours exactly, from the stroked paths and the tick marks. With their setup (see below), our
threshold-statistics bands for BPL(4,3,1) at k* ≥ 10^7.5 match theirs as follows (in log10 A):

| template | ngpbh − 2306.17149, f_PBH = 1 | f_PBH = 10⁻³ |
|---|---|---|
| Gaussian (non-linearities only) | −0.011 to −0.026 | −0.014 to −0.027 |
| f_NL = −2 | −0.005 to −0.024 | −0.008 to −0.024 |
| USR β = 3 (their quasi-inflection point) | −0.005 to −0.010 | −0.006 to −0.010 |
| curvaton r_dec = 0.9 | −0.035 to −0.048 | −0.025 to −0.042 |

This is 1–6% in A, or 8–11% for the curvaton. Near k* ≈ 10^6.5 their bands dip by up to 0.05 dex:
that is the QCD softening of the formation parameters, which we do not model. Reaching this
agreement required two conventions that 2306.17149 use but do not state where one would look:

* **The threshold: Musco et al.'s superhorizon smoothing.** Their eqs. (S4)–(S5) write
  P^T = T²(k, r_m)P, the notation of their horizon-crossing variances. Their thresholds, however,
  follow Musco et al. (2011.03014), who apply T on superhorizon scales (r̂_m ≃ 10 r_H), where
  T ≈ 1. We read their Fig. S3 (δ_th against the lognormal width and the BPL slope) from the
  vector paths (`scripts/check_threshold_2306_17149.py`, `results/threshold_2306.17149_check.md`,
  `results/fig_threshold_check.png`):
  * T at r̂_m/10, or no T at all, reproduces it to ≤ 0.0005 for 0.1 ≤ Δ ≤ 1 and ≤ 0.0016 for
    1.9 ≤ β ≤ 3;
  * T at horizon crossing is off by up to 0.08.

  The horizon-crossing reading would give BPL(4,3,1) δ_th = 0.570 instead of 0.544 and an A_* 18%
  higher (0.071 dex), for every statistic, because the linear threshold
  C_l,th = 2Φ[1 − √(1 − C_th/Φ)] moves from 0.761 to 0.826. Contacting the authors is therefore
  not needed to settle what they computed.
* **f_PBH normalisation.** Their eq. (14) with 7.9×10⁻¹⁰ gives 12.7 times the f_PBH that entropy
  conservation gives for the same β. The first-principles version, which agrees with Carr et al.
  2010 to 2%, has 1.0×10⁻⁸ in place of 7.9×10⁻¹⁰. Correcting it raises A_* by about 12%. The
  M_eq formula of 2211.01728 and pbhstat overestimates f_PBH 3–5 times, i.e. it puts A_* about
  8% too low.

The Phase 5 baseline uses their conventions (variant `base`), so that our numbers can be set
against their claims. Variants `T` (transfer-function threshold, the Phase 1–4 default) and
`entropy` (corrected f_PBH) show how large these systematics are.

## Results

All tables: σ_HPD / σ_P, or σ_P alone, for the region where PBHs are not overproduced. Full
tables are in `results/phase5_summary.md` (every shape, template, statistic and variant) and
`results/phase5_compare.md`. Figures:
* `results/fig_phase5_tension.png`: posteriors and lines;
* `results/fig_phase5_sigma.png`: σ_P per template and statistic;
* `results/fig_phase5_validation.png`: bands and contours against 2306.17149;
* `results/fig_sigw_ng_terms.png` and `results/fig_phase5_ngsigw.png`: the NG corrections to the
  SIGW spectrum and their effect (section 5; tables in `results/phase5_ngsigw_summary.md`).

### 1. In 2306.17149's own setup

This uses the BPL-marginalised posterior with the BPL(4,3,1) lines, as overlaid in their Fig. 3.

| template | f_PBH | their line (digitised) | ngpbh, threshold | ngpbh, peaks 3b |
|---|---|---|---|---|
| Gaussian (NL only) | 1 | 1.47 / 2.35 | 1.43 / 2.32 | 1.93 / 2.77 |
| Gaussian (NL only) | 10⁻³ | 2.08 / 2.92 | 2.05 / 2.88 | 2.63 / 3.24 |
| f_NL = −2 | 1 | 0.88 / 1.57 | 0.88 / 1.57 | 1.18 / 1.86 |
| curvaton r_dec = 0.9 | 1 | 1.01 / 1.66 | 1.07 / 1.72 | 1.28 / 2.03 |
| USR β = 3 | 1 | 3.52 / 3.97 | 3.45 / 3.92 | 3.78 / >4 |

* **QCD omission:** their lines (with QCD) and ours (without) give the same significance to
  ≤ 0.07 on this posterior and to ≤ 0.2 on the fixed-shape one (`phase5_compare.md`). Leaving out
  QCD does not change the conclusions, as they say.
* **Their claims:** "USR at 3σ" is confirmed. "Gaussian at over 2σ" holds for f_PBH ≤ 10⁻³ (2.1σ_HPD)
  and by posterior mass (2.3σ_P). At the f_PBH = 1 edge our marginalised posterior gives
  1.5σ_HPD, because its lower-left tip at k* ≈ 10^7 reaches 0.1 dex lower than theirs: their 2σ
  contour sits at our 94% level, and their priors on α, β, γ are not stated.
* **Peaks theory:** it raises every entry by 0.2–0.6σ. Their footnote is right that threshold
  statistics are conservative, and the size of the effect is now known.

### 2. Consistent analysis: each shape against its own posterior

Base variant, f_PBH = 1, σ_P from threshold statistics → peaks 3b:

| template | BPL(4,3,1) | LN Δ=0.5 | LN Δ=1 | LN Δ=2 |
|---|---|---|---|---|
| Gaussian | 2.98 → 3.62 | 3.02 → 3.52 | 3.54 → >4 | 3.84 → >4 |
| f_NL = −2 | 1.87 → 2.18 | 1.88 → 2.31 | 3.43 → 3.84 | >4 → >4 |
| f_NL = −1 | 1.98 → 2.40 | 1.90 → 2.36 | 3.05 → 3.46 | 3.89 → >4 |
| f_NL = −0.5 | 2.25 → 2.93 | 2.22 → 2.89 | 3.09 → 3.51 | 3.60 → >4 |
| f_NL = +0.5 | 3.67 → >4 | 3.60 → >4 | >4 | >4 |
| f_NL = +1, +2 | >4 | >4 | >4 | >4 |
| curvaton r_dec = 0.9 | 2.01 → 2.44 | 2.00 → 2.43 | 3.09 → 3.49 | 3.91 → >4 |
| USR μ_* = 5/2 | 3.53 → >4 | 3.49 → 3.91 | >4 | >4 |
| USR μ_* = 5/6, β = 3 | >4 | >4 | >4 | >4 |

* **Statistics shift:** over the 81 cases with σ_P ≤ 4, threshold → peaks 3b adds
  Δσ_P = +0.24 to +0.78 (median +0.42). With σ_HPD the shift is +0.0 to +1.1 (median +0.40).
  The shift has the same sign and nearly the same size for every NG template, spectral shape and
  convention. This is Phase 4's double ratio R(NG)/R(Gaussian) ≈ 1 carried over to the PTA.
* **3a vs 3b:** the Gaussian-peak approximation (3a) and the true NG peaks (3b) differ by
  |Δσ_P| ≤ 0.04. Phase 3b's correction to the peak count is real (up to 18%) but irrelevant for
  this question.
* **NG shift:** relative to the Gaussian case, threshold statistics, Δσ_P runs from −1.15
  (f_NL = −2, narrow spectra) to more than +2 (positive f_NL and USR templates). For broad
  spectra the sign of the f_NL < 0 shift flips: at Δ = 2, f_NL = −2 raises the tension by
  about +1σ, because large negative ζ_g also produce large compaction.
* **Shape-marginalised vs fixed posterior:** fixing the shape gives higher significance. For
  BPL(4,3,1), the steep causality tail disfavours the low-k* corner, where the allowed region is
  closest.

### 3. The statistics shift against the other choices

Gaussian ζ, threshold → peaks 3b, A_* at k* = 10⁸ Mpc⁻¹:

| choice | Δ log10 A_*, BPL(4,3,1) | Δ log10 A_*, LN Δ=1 | σ_P, BPL(4,3,1) | σ_P, LN Δ=1 |
|---|---|---|---|---|
| base (2306.17149 conventions) | −0.099 (threshold → 3b) | −0.096 | 2.98 → 3.62 | 3.54 → >4 |
| transfer function in the threshold | +0.071 | +0.208 | 2.50 → 3.22 | 2.69 → 3.12 |
| entropy f_PBH normalisation | +0.056 | +0.054 | 2.56 → 3.34 | 3.30 → 3.76 |
| f_PBH = 10⁻³ instead of 1 | −0.119 | −0.115 | 3.74 → >4 | >4 |

The threshold prescription matters most for broad spectra. For the lognormal at Δ = 1, the
transfer function moves δ_th from 0.458 to 0.537, and A_* by 62%.

### 4. Robustness (Phase 6 audit, `results/phase6_robustness.md`)

Each row changes one choice; the statistics shift Δσ_P (threshold → peaks 3b) is taken over the
cases with σ_P ≤ 4:

| change | statistics shift Δσ_P | largest change of any σ_P |
|---|---|---|
| none (reference) | +0.31 to +0.68 | — |
| prior log10 k* ∈ [6, 9], or log10 A down to −4 | unchanged | 0.00 |
| posterior grid twice as fine | unchanged | 0.00 |
| KDE tail smoothing 0 / 0.03 / 0.10 dex | +0.30 to +0.72 | 0.14 |
| 5 frequency bins (2–10 nHz) | +0.27 to +0.72 | 1.45 |
| 30 frequency bins (2–59 nHz) | +0.23 to +0.55 | 1.06 |
| peak volume (2π)^{3/2} r³ (Young) | +0.40 to +0.81 | 0.14 |
| Press–Schechter factor 2 | +0.27 to +0.59 | 0.09 |
| pbhstat v0.0.6 normalisation | +0.16 to +0.39 | 0.33 |
| high-peak formula of 2406.20089 (indicative) | +0.82 to +1.29 | 0.65 |

* **Statistics shift.** It survives every change: +0.2 to +0.8σ across the peak-normalisation
  conventions in use. Only the high-peak formula of 2406.20089 doubles it, to about +1σ.
* **Absolute significance.** It is set mostly by the frequency range. For the Gaussian
  BPL(4,3,1) case, threshold statistics give 1.9σ_P with 5 bins, 3.0σ_P with the conventional
  14, and 3.8σ_P with 30. The reason is that the posterior's tail nearest the allowed region
  (k* ≈ 10^7) puts the SIGW peak at 10–20 nHz, and the bins above it bound it.
* **Not sensitive:** the prior box and the grid resolution.

### 5. NG corrections to the SIGW spectrum

Relative to the Gaussian term, the corrections are large. In the NG15 band of the posterior bulk
(k/k* ≈ 0.02–0.3), T3/T2 ≈ 2–6 and T4/T2 ≈ 0.4–5 (`results/fig_sigw_ng_terms.png`). They grow
with the width of the peak, because the one-loop spectrum collects the whole variance σ² = A
while the peak of a broad lognormal is low. For |f_NL| = 2, Ω_GW at fixed A is therefore 10–40%
larger near the PBH lines (A ≈ 0.03–0.05). In the posterior bulk (A ≈ 0.1–0.3) the corrections
exceed the Gaussian term. The posterior moves to smaller A. For |f_NL| = 2 the median of log10 A
shifts by −0.12 (BPL), −0.14 (LN Δ = 0.5), −0.18 (LN Δ = 1) and −0.08 dex (LN Δ = 2). It shifts
by −0.07 to −0.16 dex for USR β = 3, −0.04 to −0.11 for the curvaton and ≤ 0.02 for
|f_NL| = 0.5.

Smaller A lies closer to the PBH-allowed region, so the exclusion weakens. Over the cases with
σ_P ≤ 4, in all variants:

| template | change of σ_P (NG SIGW − Gaussian SIGW) |
|---|---|
| f_NL = ±2 | −0.10 to −0.26 (median −0.17) |
| curvaton r_dec = 0.9 | −0.03 to −0.20 (median −0.14) |
| f_NL = ±1 | 0.00 to −0.09 (median −0.05) |
| f_NL = ±0.5, USR μ_* = 5/2 | ≤ 0.03 |
| USR β = 3, μ_* = 5/6 | σ_P > 4 with either SIGW (shifts −0.16 to +0.06) |

Base variant, f_PBH = 1, threshold → peaks 3b: f_NL = −2 for BPL(4,3,1) goes from 1.87 → 2.18 to
1.73 → 2.05, LN Δ = 0.5 from 1.88 → 2.31 to 1.73 → 2.11, LN Δ = 1 from 3.43 → 3.84 to
3.25 → 3.67. The curvaton for BPL(4,3,1) goes from 2.01 → 2.44 to 1.92 → 2.30.

* **The statistics shift is unchanged.** Over the 66 reliable pairs, threshold → peaks 3b adds
  +0.23 to +0.67 with the NG SIGW (median +0.39), against +0.24 to +0.68 (median +0.42) with the
  Gaussian SIGW. The pair-by-pair change is −0.10 to +0.10 (median 0.00; 90% within ±0.06).
* **Truncation.** For the quadratic template, dropping the A⁴ terms changes σ_P by ≤ 0.02. For
  the curvaton and USR maps, adding the F⁴A⁴ terms changes it by ≤ 0.02 too. The neglected O(A⁴)
  terms in G and H are of similar size, so the A³ truncation is adequate here.
* **Against the NG effect on the PBH line.** Through the PBH abundance, NG shifts σ_P by −1.2 to
  beyond +2σ; through the SIGW, by at most 0.26σ. The SIGW part relieves the tension for either
  sign of f_NL, because at fixed A it raises Ω through F² (and G > 0). In the reliable cases σ_P
  never rises by more than 0.01.

## Answer to the research question

The gap between threshold statistics and peaks theory changes how strongly PBH overproduction
disfavours a SIGW explanation of NANOGrav 15, but modestly and uniformly:

* **Size and sign.** Peaks theory strengthens the exclusion by about 0.4σ (0.25–0.8σ) for every
  NG template and spectral shape tested.
* **No interaction with NG.** The shift is nearly independent of NG, so the two effects add.
  Choosing peaks theory never turns an allowed scenario into an excluded one or the reverse,
  except near a chosen significance threshold.
* **Against NG.** For |f_NL| ≥ 1 and the USR and curvaton templates, the NG effect is 2–5 times
  larger; for |f_NL| ≈ 0.5 the two are comparable. Negative f_NL and the curvaton relieve the
  tension by up to 1.2σ for narrow and moderately broad spectra; positive f_NL and the USR
  templates push it beyond 4σ. Peaks theory erodes the relief from negative NG by 0.3–0.5σ. In
  2306.17149's setup, f_NL = −2 stays below 2σ_P with either statistic; the curvaton moves from
  1.7 to 2.0σ_P.
* **Absolute significance.** The largest analysis choice is the frequency range: 5, 14 or 30
  bins move every σ by −1.5 to +1σ. The statistics shift does not depend on it.
* **Against the conventions.** Three convention choices are rarely stated. Two of them shift
  the significance about as much as the choice of statistics does:
  * the transfer function in the threshold: −0.2 to −1.0σ;
  * f_PBH = 1 versus 10⁻³: +0.4 to +0.8σ.

  The third, the corrected f_PBH normalisation, shifts it less: −0.15 to −0.4σ. Quoting a
  significance without these choices leaves about ±1σ unaccounted for.
* **Peak approximation.** Using the Gaussian-peak approximation instead of the true NG peaks
  changes nothing at this level.
* **NG in the SIGW spectrum.** It relieves the tension by up to 0.26σ for |f_NL| = 2 and the
  curvaton, and by ≤ 0.1σ for |f_NL| ≤ 1. It leaves the statistics shift unchanged (±0.1σ). It
  matters more than the choice of statistics only for |f_NL| ≳ 2.

## What is not included

* **QCD era in the formation parameters.** 2306.17149 make δ_c, K, γ and Φ depend on M_H through
  the QCD crossover. Their bands therefore dip below ours by up to 0.05 dex at k* ≈ 10^6.5, i.e.
  M_H ≈ 1 M_⊙. The NG15 posterior lies at k* ≳ 10^7, where the dip has ended; `phase5_compare.md`
  quantifies the effect by putting their own lines on our posterior.
* **QCD era in the SIGW.** The sound-speed change during the QCD era modifies the infrared tail of
  Ω_GW (Abe et al. 2020, Franciolini et al. 2023). It is left out, as in 2306.17149.
* **NG corrections to Ω_GW beyond the quadratic template.** For the curvaton and USR maps the
  SIGW series stops at A³ (section 5). The comparison with 2306.17149 in section 1 uses the
  BPL-marginalised posterior, which exists only for the Gaussian SIGW, because the NG terms are
  tabulated for fixed shapes. The fixed-shape results of section 5 show the size of the change.
  (An earlier version of this document estimated the corrections as ~A(3f_NL/5)², i.e. ≲ 6% at
  the line. The computed prefactor is 2–8 rather than 1, which is why the shift reaches 0.26σ.)
* **Fixed shape per line.** Each line is compared with the posterior of the same shape. The
  posterior marginalised over the BPL shape, which is what 2306.17149 overlay on their bands, is
  used only for the comparison with their claim (`scripts/pta_bpl_marginal.py`).
* **Free-spectrum refit instead of the full PTA likelihood.** The refit reproduces NANOGrav's own
  fits within ~0.3σ (power law) and 2306.17149's SIGW fit (Pipeline item 2).
* **Priors.** Credible levels depend on the prior box: uniform in log10 k* ∈ [5, 10] and
  log10 A ∈ [−3, 0]. The profile Δχ² in `phase5_tension.csv` does not.
* **No astrophysical foreground in the tables above.** As in 2306.17149, the whole signal is
  assumed to be SIGW. `scripts/pta_smbhb.py` adds an SMBHB power law with NANOGrav's
  astrophysical prior (`results/pta_smbhb.md`): σ_P is then capped at 2.3–2.8σ, and the
  statistics shift drops to +0.04 to +0.39σ (median +0.26, against +0.43 without).
