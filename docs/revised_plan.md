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
  transfer function. The no-transfer and 2/5 variants are systematics to show, not alternatives
  to mix. (Phase 5 found that 2306.17149 themselves used the no-transfer threshold, see
  `docs/phase5_pta.md`; Phase 5 therefore runs both.)

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

## Phase 3: peaks with NG for broad spectra — done (robustness items left)

Definition: peaks = maxima of the compaction field at the formation scale (prescription (i) in
the assessment), so that peaks and threshold statistics differ only by the peak constraint.

1. **3a, Gaussian-peak approximation — done.** BBKS maxima of the Gaussian field C_g(x; r) with
   ζ_r as a correlated extra scalar (exact BBKS f(x)), mapped through F′. Exact for F′ = const;
   Young's approach generalised to broad spectra.
2. **3b, true NG peaks — done** (`docs/phase3b_derivation.tex`, `PeaksConfig(constraint="C_l")`).
   Maxima of C_l = F′(ζ_r) C_g, a local function of two correlated Gaussian fields. Germani et
   al. 2503.05434 cover a function of one Gaussian field, where peak positions do not move; no
   paper imposes the peak condition on the non-Gaussian compaction. Result: ∇C_l = F′∇Y with
   Y = C_g + λζ_r, λ = C_l F″/F′², so the maxima are BBKS maxima of Y, times a rank-one Hessian
   factor 1 − τR₂ that is exact to first order. That factor is identically 1 for the logarithmic
   template, because 1/F′ is linear.
3. **Monte Carlo check, part 1 — done.** Independent sampling of the 20 Gaussian variables with
   the exact Hessian (`ngpbh.validation`) reproduces the semi-analytic density to 3×10⁻⁴ at 15 of
   16 test nodes. The exception is the most extreme node (f_NL = −2, Δ = 1, density already
   suppressed twentyfold), off by 0.85%, the O(τ²) remainder.
4. **Monte Carlo check, part 2 — done** (`scripts/lattice_peaks.py`,
   `results/phase3b_lattice.csv`, `results/fig_phase3b_lattice.png`). Maxima of C_l counted on
   256³ Gaussian random fields (Δ = 1, r = 8 cells, 24 realisations, ν = 1.5–3), with predictions
   from the lattice's own mode sums:
   * Gaussian calibration: lattice/BBKS = 1.010 → 0.976, a mild grid trend.
   * f_NL = −2: lattice/3b = 1.006 → 0.949 (±0.3–2.7%), against lattice/3a = 1.18 → 1.14.
   * f_NL = +2: lattice/3b = 1.011 → 0.997, against lattice/3a = 0.91 → 0.98.
   * μ_* = 5/6: lattice/3b = 1.007 → 0.998, against lattice/3a = 0.93 → 0.98.

   3b agrees to about 1% beyond the grid trend; 3a miscounts peaks by up to 18%. A 512³ box was
   not needed for Δ = 1 with r = 8 cells (L/r = 32).
5. **Still open — robustness:**
   * the Germani–Sheth / Fumagalli non-linear statistics (pbhstat `NonLinearModel`) in the
     Gaussian limit;
   * profile dispersion (Escrivà 2607.08738) for f_NL < 0;
   * a lattice run at Δ = 2 (needs L/r ≥ 60, i.e. 512³ at r = 8).

## Phase 4: discrepancy map R — done with Phase 3b peaks

* Grid: Δ ∈ {0.1, 0.3, 0.5, 1, 2}; f_NL ∈ {−2, −1, −0.5, 0, 0.5, 1, 2}; μ_* ∈ {5/6, 5/2, 10, −5/2}.
  Statistics: threshold, peaks 3a, peaks 3b, all with the same threshold, window, transfer
  function, k < 4.49/r cut, K and γ (`scripts/scan_amplitudes.py`, ~20 min on 3–4 cores).
  Tables: `python scripts/summarise_scan.py`; figure: `results/fig_phase3_ratio.png`.
* Exit: R unchanged when every grid is doubled. This is built into the tests at 2×10⁻⁴ in A_*
  (3a) and 10⁻⁴ in log β (3b).
* Gaussian R is **0.77–0.79** in ngpbh's BBKS normalisation. It spans 0.75–0.89 across the
  normalisation conventions in use, and ~0.6–0.67 with the 2406.20089 high-peak formula
  (2602.12252). So report the conventions with every number, and lead with the double ratio,
  which cancels them to < 0.01.
* **Result 1 — true NG peaks vs the Gaussian-peak approximation:**
  * A_3b/A_3a = 0.983–1.000 over the grid. The largest shift is −1.7% (f_NL = −2, Δ = 0.5);
    most entries are under 0.5%.
  * 3a miscounts peaks by up to 18% (lattice), but β depends exponentially on A, so the
    required amplitude moves only by δln A ≈ δln n/|ln β|.
  * The approximation used implicitly in the literature is therefore good to ≲ 2% in A_* for
    |f_NL| ≤ 2. That is a new, quantified statement.
* **Result 2 — the R map (3b):**
  * R = 0.738–0.817 for |f_NL| ≤ 2.
  * The double ratio R/R_Gauss stays within 0.958–1.036 for |f_NL| ≤ 2: positive NG pushes R
    up, and negative f_NL pulls it down for narrow spectra and up for broad ones.
  * The non-perturbative USR tail (μ_* = 5/6, f_NL = 1 at second order) gives 1.05–1.16, largest
    at Δ = 2, against 1.02–1.03 for quadratic f_NL = 1.
* **Framing for the note:** across the grid, NG moves the threshold-statistics A_* by factors of
  0.2–2. The choice of statistics moves it by 18–26% (1 − R for |f_NL| ≤ 2). The interaction
  between the two is ±4% for |f_NL| ≤ 2 and up to 16% for strongly non-perturbative NG. For the PTA question this
  means the statistics systematic is close to an NG-independent ~20% rescaling of the
  overproduction line in A, except for strong USR-type tails. Phase 5 found that this rescaling
  moves the NG15 exclusion by +0.4σ (0.25–0.8σ), the same for every NG template.

## Phase 5: PTA tension — done (`docs/phase5_pta.md`)

1. **SIGW** (`ngpbh.sigw`) for Gaussian ζ, with Standard-Model g*(T) (`ngpbh.thermal`). Checks:
   flat spectrum 0.8222436 A², the Dirac closed form, scipy `dblquad` to 10⁻⁸–10⁻¹², and the k³
   tail. NG corrections to Ω_GW are left out, as in 2306.17149 (≲ 3% in A at the overproduction
   line).
2. **NANOGrav 15 likelihood** (`ngpbh.pta`): the ceffyl free-spectrum refit, 14 bins, data from
   Zenodo with MD5 checked. It reproduces:
   * NANOGrav's power-law fits within ~0.3σ;
   * 2306.17149's lognormal and BPL SIGW fits, with medians within 0.05;
   * their contour levels (68.3/95.4/99.7%).
   The KDE tails needed regularising for 3–4σ statements (compact kernels).
3. **Exit check ✓.** Their Fig. 3 was read from the vector paths of the arXiv source. Our bands
   for BPL(4,3,1) agree to 1–6% in A (curvaton 8–11%) for k* ≥ 10^7.5. They dip at k* ≈ 10^6.5
   because of QCD, and their own lines give the same significances as ours on the same posterior
   (≤ 0.07), so the QCD omission is quantified as negligible. Two unstated conventions of theirs
   were found:
   * the threshold is computed without the transfer function (+18% in A_* with it for
     BPL(4,3,1), +62% for the lognormal at Δ = 1);
   * the f_PBH normalisation is 12.7× the entropy-conservation value (−12% in A_*).
4. **Result.** Threshold → peaks 3b strengthens the exclusion by Δσ_P = +0.24 to +0.78 (median
   +0.42) for every template and shape:
   * Gaussian BPL(4,3,1): 3.0 → 3.6σ_P;
   * their shape-marginalised setup: 2.3 → 2.8σ_P.

   NG moves it by −1.2σ to beyond +3σ. The transfer-function threshold choice moves it by −0.2 to
   −1.0σ, f_PBH = 1 vs 10⁻³ by +0.4 to +0.8σ, and the f_PBH normalisation by −0.15 to −0.4σ.
   Peaks 3a and 3b agree to 0.04σ.
5. **Phase 6 audit** (`results/phase6_robustness.md`):
   * The statistics shift is +0.2 to +0.8σ across the peak-normalisation conventions (+0.8 to
     +1.3σ with the 2406.20089 high-peak formula). It is insensitive to the prior box, the grid
     and the KDE smoothing.
   * The absolute significance moves by −1.5 to +1σ with the frequency range (5/14/30 bins).
6. **Still open:**
   * the SMBHB foreground;
   * the full PTA likelihood instead of the refit;
   * QCD sound-speed effects on the SIGW tail;
   * a joint NANOGrav + EPTA fit.

## Phase 6: package and write up — scaffold in place

Repo layout follows the plan (`nongaussian.py`, threshold and peaks on a shared integrator in
`stats.py`, `sigw.py`, `pta.py`, `tests/`). 74 tests run in ~40 s. Remaining: CI, a Zenodo DOI or a pbhstat PR, the
10–15 page note, and an arXiv endorser.

## Pitfalls, confirmed or added

* Exponential sensitivity: log β throughout ✓ (a 1% shift in A moves β by ~50% at ν ≈ 7).
* USR template singularity at ζ_g = μ_*: integration domain cut there ✓.
* Peaks need the cut; apply it to both ✓.
* Hidden threshold differences: the transfer-function choice inside the threshold prescription
  moves C_th by 0.11 at Δ = 2 — larger than any statistics effect.
* **New:** normalisation conventions (PS factor, peak volume, Jacobian) move R by ~0.14.
* **New:** pbhstat's code and paper differ; compare codes in one convention at a time.
* **New:** the rank-one Hessian term of Phase 3b is formally O(σ) but grows like λ², so it
  reaches tens of percent in the peak density for f_NL = −2 at Δ = 1; it cannot be dropped.
* **New (Phase 5):** free-spectrum KDEs with compact kernels are floored outside the sampled
  range; regularise their tails before quoting 3–4σ, and do not quote beyond ~4σ.
* **New (Phase 5):** the f_PBH normalisation in common use (M_eq formula; 2306.17149 eq. 14) is
  3–13× larger than entropy conservation gives; state it with every A_*.
* **New (Phase 5):** a `pgrep -f` waiter matches its own command line; use `pgrep -f "[p]attern"`.
