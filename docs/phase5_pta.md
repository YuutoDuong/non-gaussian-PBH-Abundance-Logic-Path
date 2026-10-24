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
   * **Lognormal SIGW fit, Δ marginalised:** log10 A = −0.57 (+0.55/−0.83), against
     −0.56 (+0.53/−0.74) in 2306.17149 Fig. 2, and log10 k* = 8.16 against 8.14.
3. **Overproduction lines** (`AbundanceModel.amplitude_curve`). The variances depend only on k r,
   so moving the spectrum to another k* only re-weights the f_PBH integral over M_H. One set of
   β grids then gives A_*(k*) at all k*. It agrees with a direct solve at each k* to 10⁻⁷.
4. **Tension.** The NG15 posterior in (log10 k*, log10 A) uses uniform priors on [5, 10] × [−3, 0]
   (A ≤ 1, as in 2306.17149) and the same spectral shape as the line. For each line we compute:
   * the highest-posterior-density level at which the allowed region A < A_*(k*) is first
     reached, quoted as a two-sided 1D Gaussian equivalent ("excluded at n σ");
   * the posterior mass in that region;
   * the profile Δχ².

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

* **Threshold without the transfer function.** Their eq. (S4) writes P^T, but their Fig. S3
  thresholds (δ_th = 0.458 for the lognormal at Δ = 1, 0.544 for BPL(4,3,1)) are the Musco values
  computed without it (ngpbh: 0.4577 and 0.5438). With the transfer function, BPL(4,3,1) has
  δ_th = 0.570 and A_* is 18% higher (0.071 dex), for every statistic, because the linear threshold
  C_l,th = 2Φ[1 − √(1 − C_th/Φ)] moves from 0.761 to 0.826.
* **f_PBH normalisation.** Their eq. (14) with 7.9×10⁻¹⁰ gives 12.7 times the f_PBH that entropy
  conservation gives for the same β. The first-principles version, which agrees with Carr et al.
  2010 to 2%, has 1.0×10⁻⁸ in place of 7.9×10⁻¹⁰. Correcting it raises A_* by about 12%. The
  M_eq formula of 2211.01728 and pbhstat overestimates f_PBH 3–5 times, i.e. it puts A_* about
  8% too low.

The Phase 5 baseline uses their conventions (variant `base`), so that our numbers can be set
against their claims. Variants `T` (transfer-function threshold, the Phase 1–4 default) and
`entropy` (corrected f_PBH) show how large these systematics are.

## Results

Pending: `scripts/pta_tension.py` (overproduction lines and significances) is running; this
section is filled in from `results/phase5_summary.md` when it finishes.

## What is not included

* **QCD era in the formation parameters.** 2306.17149 make δ_c, K, γ and Φ depend on M_H through
  the QCD crossover. Their bands therefore dip below ours by up to 0.05 dex at k* ≈ 10^6.5, i.e.
  M_H ≈ 1 M_⊙. The NG15 posterior lies at k* ≳ 10^7, where the dip has ended; `phase5_compare.md`
  quantifies the effect by putting their own lines on our posterior.
* **QCD era in the SIGW.** The sound-speed change during the QCD era modifies the infrared tail of
  Ω_GW (Abe et al. 2020, Franciolini et al. 2023). It is left out, as in 2306.17149.
* **NG corrections to Ω_GW.** These are of relative size ~A (3f_NL/5)². At the overproduction line
  (A ≈ 0.01–0.04, |f_NL| ≤ 2) they change Ω by ≲ 6%, i.e. A by ≲ 3%: a fraction of the statistics
  shift. They are larger in the posterior bulk (A ≳ 0.1), but the significance is set by the
  posterior's low-A tail.
* **Fixed shape per line.** Each line is compared with the posterior of the same shape. The
  posterior marginalised over the BPL shape, which is what 2306.17149 overlay on their bands, is
  used only for the comparison with their claim (`scripts/pta_bpl_marginal.py`).
* **Free-spectrum refit instead of the full PTA likelihood.** The refit reproduces NANOGrav's own
  fits within ~0.3σ (power law) and 2306.17149's SIGW fit (Pipeline item 2).
* **Priors.** Credible levels depend on the prior box: uniform in log10 k* ∈ [5, 10] and
  log10 A ∈ [−3, 0]. The profile Δχ² in `phase5_tension.csv` does not.
* **No astrophysical foreground.** As in 2306.17149, the whole signal is assumed to be SIGW. A
  supermassive-black-hole-binary component would relax every tension quoted here.
