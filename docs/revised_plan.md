# Option A, revised (23 Oct 2026)

Same question, same seven phases as the plan PDF; the changes come from
`docs/phase0_gap_assessment.md` and from what the first code runs showed. Status is marked per
phase.

**Question.** Does the gap between threshold statistics and peaks theory, once local NG is
included for broad spectra, change how strongly PBH overproduction disfavours a scalar-induced
GW explanation of the PTA signal?

**Deliverables** (unchanged): an NG module for threshold statistics and peaks; the map
R(Δ, NG) = A_*^peaks / A_*^threshold; the PTA tension per prescription; a tested repo and a note.
Add: the double ratio R(Δ, NG)/R(Δ, Gaussian). It cancels the PS-factor and peak-volume
conventions to < 0.01, while the Gaussian R itself moves by 0.07 between them
(`scripts/double_ratio_conventions.py`, table in `docs/conventions.md`).

## Phase 0: confirm the gap — done (one item left for you)

* Gap open in the narrower form stated in the assessment. Add 2406.20089, 2607.08738,
  2503.05434, 2603.23025, 2501.00295 and 2512.13648 to the reading list.
* INSPIRE scan done (23 Oct 2026, 193 citing records): nothing closes the gap; re-run it before
  submission.
* **You still need to:** email Cole & Fumagalli (mention the PS factor 2 and the peaks σ₀/2
  normalisation in pbhstat v0.0.6, and ask about NG plans); write and date your prediction
  before opening the scan results.

## Phase 1: Gaussian baselines — done

* 1a (numerics). ngpbh vs pbhstat in pbhstat's conventions: A_* agrees to **0.017%** for
  Δ = 0.1–2, PS and peaks (`results/phase1_pbhstat_comparison_pbhstat.csv`, exit < 1% ✓).
* 1b (physics conventions). Threshold statistics on the full compaction with the transfer
  function and Musco threshold reproduce the Ferrante et al. broad-spectrum benchmark
  (A = 9.09×10⁻³ vs 9.2×10⁻³ for β = 10⁻¹⁶ at M_H = 10⁻¹⁵ M_⊙).
* Decision recorded: **one threshold prescription for both statistics** = Musco et al. with the
  transfer function (as 2306.17149). The no-transfer and 2/5 variants are systematics to
  show, not alternatives to mix.

## Phase 2: threshold statistics with NG — done

* Templates: quadratic f_NL, logarithmic (USR) μ_*, and the curvaton (added because it gives a
  published NG benchmark).
* Exit check ✓: Ferrante et al. Fig. 8 curvaton amplitudes reproduced (9.41×10⁻³ vs 9.5×10⁻³
  for r_dec = 0.5; 2.81×10⁻³ vs 2.8×10⁻³ for r_dec = 0.1). f_NL = 0 and μ_* → ±∞ reproduce the
  Gaussian result to 10⁻⁶ in log β.
* Scan over Δ × 11 templates: `results/phase2_amplitudes.csv`, figure
  `results/fig_phase2_threshold_amplitudes.png`.
* Optional, before the note: the q-marginalised prescription of 2402.11033 (q as a random
  variable, threshold → 2/5 for broad spectra), as a third threshold systematic.

## Phase 3: peaks with NG for broad spectra — started

Fix the definition first: peaks = maxima of the compaction field at the formation scale
(prescription (i) in the assessment), so that peaks and threshold statistics differ only by the
peak constraint. Then:

1. **3a, Gaussian-peak approximation — implemented.** Count BBKS maxima of the Gaussian field
   C_g(x; r) with ζ_r as a correlated extra scalar (exact BBKS f(x), no high-peak limit), map
   through F′. Exact for F′ = const; it is Young's approach generalised to broad spectra.
   First R map: `results/fig_phase3a_ratio_preview.png` (preliminary).
2. **3b, true NG peaks — the new physics.** Count maxima of C_l = F′(ζ_r) C_g, a product of two
   correlated Gaussian fields. ∇C_l = F′∇C_g + F″C_g∇ζ_r, so the stationary points move; the
   joint vector is (C_g, ζ_r, ∇C_g, ∇ζ_r, ∇∇C_g, ∇∇ζ_r), 20 Gaussian variables with known
   covariances (all follow from the same spectral integrals). Start from the local-NG peak
   statistics of Germani et al. 2503.05434 (§3) and Young 2201.13345 App. B. Write the LaTeX
   derivation before coding, as planned. Expected structure: integrate δ³(∇C_l) analytically
   (Gaussian in ∇ζ_r after eliminating ∇C_g), keep |det ∇∇C_l| with a high-peak expansion
   around the 3a result, and quantify the correction.
3. Monte Carlo check, part 1: sample the 20-variable Gaussian vector with importance sampling
   near threshold; must reproduce 3a when F″ = 0 and give the 3b result otherwise.
4. Monte Carlo check, part 2: 3D Gaussian random fields at large σ, apply F, smooth, find maxima
   of C_l, compare counts with 3a and 3b for Δ ≤ 1 (box ≥ 512³ for Δ = 1).
5. Robustness: the Germani–Sheth / Fumagalli non-linear statistics (pbhstat `NonLinearModel`)
   in the Gaussian limit; profile dispersion (Escrivà 2607.08738) for f_NL < 0.

Fallback as before (week 6): restrict to Δ ≤ 0.5 where 3a and 3b can be checked on lattices.

## Phase 4: discrepancy map R — machinery ready

* Grid: Δ ∈ {0.1, 0.3, 0.5, 1, 2}, f_NL ∈ {−2, −1, −0.5, 0, 0.5, 1, 2}, μ_* ∈ {5/6, 5/2, 10, −5/2}
  (`scripts/scan_amplitudes.py`, 4 cores, ~5 min). Same threshold, window, transfer function,
  k < 4.49/r cut, K, γ for both statistics.
* Exit: R unchanged when every grid is doubled — built into the tests at 2×10⁻⁴ in A_*.
* Sanity check updated: Gaussian R is **0.77–0.79** in ngpbh's BBKS normalisation, not "a little
  below 1"; the "slightly smaller amplitude" of 2402.11033 refers to other conventions. R spans
  0.75–0.89 across the normalisation conventions in use (`results/phase1_convention_spread.csv`),
  so report the conventions with every number, and lead with the double ratio.
* Preliminary (3a) trend: positive NG pushes R towards 1 (R = 0.92 for μ_* = 5/6 at Δ = 2);
  negative f_NL pulls it below the Gaussian value for narrow spectra.
* Preliminary (3a) sizes, which set the framing of the note: NG moves the threshold-statistics
  A_* by factors 0.2–2 across the grid, the choice of statistics moves it by ~20% (R ≈ 0.78),
  but the interaction is small: the double ratio R/R_Gauss stays within 0.96–1.04 for
  |f_NL| ≤ 2. The exception is μ_* = 5/6 (same f_NL = 1 at second order, but the full USR tail):
  1.05–1.17, largest at Δ = 2, against 1.02–1.03 for quadratic f_NL = 1, so the
  non-perturbative tail, not f_NL, drives it. That smallness is expected, because in 3a the
  peak constraint acts on the Gaussian field C_g. So whether NG changes the threshold–peaks gap
  materially is decided by 3b; if 3b also gives a few percent, the result is still publishable
  ("the statistics systematic is NG-independent to x%"), but the PTA headline changes.

## Phase 5: PTA tension — not started

1. SIGW spectra on a grid in (A, k_*, Δ) with SIGWfast or SIGWAY, Gaussian SIGW for all NG cases
   (justified by Atal & Domènech 2103.01056 when PBHs are not overproduced).
2. **Exit check (revised):** reproduce the f_PBH = 1 band of 2306.17149 Fig. 3 (Gaussian, BPL
   α = 4, β = 3, γ = 1, NANOGrav15) to a few percent in A, with their fitting method (direct fit to
   the 14 free-spectrum bins) and their significance definition (credible-region level). Note
   that they include QCD-dependent K, γ, C_th, Φ; either add them (Musco–Jedamzik–Young
   2303.07980 tables) or quote the size of the omission.
3. Then PTArcade / ceffyl, the posterior mass in the allowed region per prescription and NG case,
   and the question: does the choice of statistics move the significance more than the NG does?
4. Network: arXiv/INSPIRE are blocked from the cloud environment used so far; the NANOGrav data
   and PTArcade need a machine with access.

## Phase 6: package and write up — scaffold in place

Repo layout follows the plan (`nongaussian.py`, threshold and peaks on a shared integrator in
`stats.py`, `tests/`). 47 tests run in ~35 s. Remaining: CI, a Zenodo DOI or a pbhstat PR, the
10–15 page note, and an arXiv endorser.

## Pitfalls, confirmed or added

* Exponential sensitivity: log β throughout ✓ (a 1% shift in A moves β by ~50% at ν ≈ 7).
* USR template singularity at ζ_g = μ_*: integration domain cut there ✓.
* Peaks need the cut; apply it to both ✓.
* Hidden threshold differences: the transfer-function choice inside the threshold prescription
  moves C_th by 0.11 at Δ = 2 — larger than any statistics effect.
* **New:** normalisation conventions (PS factor, peak volume, Jacobian) move R by ~0.14.
* **New:** pbhstat's code and paper differ; compare codes in one convention at a time.
