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
  to mix. (Phase 5 found that 2306.17149's thresholds follow Musco et al.'s superhorizon smoothing,
  which here equals no transfer function to ≤ 0.0005 in δ_th, and that "with the transfer function"
  above means T at horizon crossing, a stronger smoothing than Musco et al.'s. Phase 5 runs both;
  R changes by ≤ 0.5% for quadratic NG, ≤ 6% for strong USR tails.)

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
   tail. The primordial-NG corrections to Ω_GW are item 6 (2306.17149 leave them out).
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
   * their threshold follows Musco et al.'s superhorizon smoothing (≈ no transfer function),
     although their eqs. (S4)–(S5) write the horizon-crossing T(k, r_m). Read literally, the
     formula gives +18% in A_* for BPL(4,3,1) and +62% for the lognormal at Δ = 1. Confirmed
     from their Fig. S3 over its whole range (`scripts/check_threshold_2306_17149.py`);
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
6. **NG corrections to the SIGW spectrum — done** (`ngpbh.sigw_ng`; `docs/phase5_pta.md`
   section 5). All seven terms of Adshead, Lozanov & Weiner for ζ = ζ_g + F(ζ_g² − σ²) are
   included: the disconnected ones by quadrature of the one-loop spectrum, the connected ones by
   randomised QMC. The series is exact at A⁴ for the quadratic template and truncated at A³ (F, G)
   for the curvaton and USR maps. Validation: the MC reproduces the quadrature terms (rms pull
   0.97), and all six NG terms reproduce Adshead et al. Fig. 4 to ≤ 3% of each curve's maximum.
   * **Result.** At fixed A, Ω_GW is larger by F²A × (2–8): the posterior moves to smaller A by up
     to 0.18 dex, and σ_P drops by 0.10–0.26 for |f_NL| = 2 and 0.03–0.20 for the curvaton
     (≤ 0.1 for |f_NL| ≤ 1).
   * The statistics shift is unchanged (median +0.42 → +0.39, each pair within ±0.1).
   * The earlier estimate "≲ 3% in A" was too small by a factor 2–8: the one-loop power collects
     the whole variance.
7. **Still open:**
   * ~~the SMBHB foreground~~: done for the journal version (see "Journal computations");
   * the full PTA likelihood instead of the refit;
   * QCD sound-speed effects on the SIGW tail;
   * a joint NANOGrav + EPTA fit.

## Phase 6: package and write up — draft done (24 Oct 2026)

* **Audit** (`results/phase6_robustness.md`): the results hold.
  * The statistics shift in the PTA significance is +0.2 to +0.8σ under every analysis change
    tried.
  * The absolute significance carries ±1–1.5σ from the frequency range, the threshold
    prescription and the f_PBH conventions. These are now documented.
  * One packaging bug fixed: NumPy ≥ 2 is required (`np.trapezoid`).
* **Packaging:**
  * GitHub Actions CI (pytest on 3.10/3.12/3.13, NG15 data cached);
  * MIT `LICENSE`;
  * `CITATION.cff` and `.zenodo.json` (author: Duong Ngoc Khoa, Phenikaa Institute for Advanced
    Study, ORCID 0009-0009-9409-3298);
  * version 0.4.0.
* **Note:** `docs/note/ngpbh_note.tex` / `.pdf`, 18 pages. It covers the construction, the
  checks, the R map, the PTA tension, the NG corrections to the induced GWs, robustness, the two
  pitfalls found in the literature, and appendices.
* **Pre-publication check (25 Oct 2026):**
  * every reference checked against INSPIRE: one author list fixed (Blas et al. has seven
    authors), ten journal references added;
  * 27 standard references added: the other PTAs, the compaction function, peaks theory for
    broad spectra, Yoo, Gong & Yokoyama 1906.06790 and Atal, Garriga & Marcos-Caballero
    1905.13202 (peaks of the Gaussian component with local NG, the second for the USR log
    template), critical collapse, the NG-SIGW literature and the trispectrum no-go theorem;
  * "Draft" removed; the note cites the concept DOI and version 0.6.0; an acknowledgement
    discloses the use of an AI assistant (required by arXiv and most journals);
  * version 0.6.0 prepared (`pyproject.toml`, `CITATION.cff`, `.zenodo.json`, README).
* **Mock referee review (25 Oct 2026).** Not a substitute for an independent expert read; the
  points a referee would most likely raise, and what was done:
  * the UV cut kr < 4.49 of the peak moments: R changes by ≤ 1.2% (R/R_Gauss ≤ 0.8%) when the
    cut is removed or moved to 3.5, 6 or 7.725 (`scripts/cut_sensitivity.py`); now in Sec. 4;
  * which peaks: R is defined for BBKS maxima of the compaction at fixed scale; typical-profile
    peaks of ζ (Yoo et al.) and the Germani–Sheth statistics are not computed — now stated in
    the limitations (still open, see Phase 3);
  * the threshold is the same for every NG template (NG deforms the profile: Kehagias, Musco &
    Riotto 1906.07135; Escrivà et al. 2202.01028), and scales are counted without a
    cloud-in-cloud correction — both common to the two statistics, now stated;
  * the typeset induced-GW formula integrates to 0.8222436 A² for a flat spectrum, as Kohri &
    Terada; "≤ 2%" for the cut in threshold statistics corrected to ≤ 2.3%;
  * the abstract is plain prose (arXiv and journals need it), and "used so far" is replaced by
    what the Gaussian-peak approximation does;
  * still a referee risk: App. B (f_PBH normalisation, factor 12.7 against eq. 14 of
    2306.17149) — be ready to defend it; profile dispersion for f_NL < 0 (Escrivà 2607.08738).
* **Remaining (your decisions or actions):**
  * ~~a tagged release and the Zenodo DOI~~: done, tag `astro` (commit 692dc2c),
    DOI 10.5281/zenodo.23624507 (all versions: 10.5281/zenodo.23624506);
  * ~~a release `v0.6.0`~~: done, tag `v0.6.0` (commit 401f560), DOI 10.5281/zenodo.23643762;
    the note cites it, and the affiliation is now Phenikaa Institute for Advanced Study,
    Phenikaa University, Hanoi, Vietnam;
  * read the AI-use acknowledgement in the note and adjust the wording;
  * the note on Zenodo as a preprint: a separate record (Publication → Preprint, CC BY 4.0),
    related identifier "is supplemented by" 10.5281/zenodo.23624506;
  * an arXiv endorser (astro-ph.CO); the arXiv upload needs the figures next to the `.tex`;
  * ~~the courtesy notes to the 2306.17149 authors and to Cole & Fumagalli~~: not needed for the
    facts (their Fig. S3, eq. 14 and the code comparison settle them); you decided against them
    on 25 Oct 2026;
  * ~~re-running the INSPIRE scan~~: re-run 25 Oct 2026 (records to 20 Oct), nothing closes the
    gap; repeat just before the arXiv submission;
  * optional: a pbhstat PR with the NG threshold statistics and peaks.
  * for a journal: the improvements in `docs/journal_review.md`. Its quick fixes are done (25 Oct
    2026): the Iovino et al. framing, Appendix B (origin of the 12.7, cross-checks), the spectra
    written out, σ_P in the table captions, QCD references, and descriptive labels instead of
    "3a/3b" in the text and figures (66 references, 19 pages). Its computations are done too
    (see below);
  * a release `v0.7.0` (GitHub release, which Zenodo archives): the note now quotes numbers that
    v0.6.0 does not contain; send the new version DOI so the note can cite it.

## Journal computations — done (25 Oct 2026)

Items 1–4 of `docs/journal_review.md`, now in the note. They need code version 0.7.0 (prepared;
the release and its DOI are still to come).

* **R_Gauss and the factorisation** (`scripts/r_gauss_fit.py`):
  * lognormal, 0.1 ≤ Δ ≤ 2: R_Gauss = 0.783 + 0.0066 ln Δ, to 0.002;
  * BPL(4,β,1), β = 2–4: 0.775–0.778;
  * R/R_Gauss for BPL(4,3,1) is 0.971–1.033 over the NG templates (1.03 for the curvaton), so
    A*_peaks ≃ R_Gauss A*_threshold to about 4% (Sec. 4 and the conclusions).
* **Other prescriptions** (literature, Gaussian ζ; Sec. 4):
  * BBKS peaks need 9–15% smaller amplitudes than Press–Schechter in Gow et al. 2021 and 11–12%
    in the pbhstat paper;
  * the high-peak formula gives R ≈ 0.6–0.67;
  * the non-linear (Germani–Sheth) statistics need 4–12% larger amplitudes in pbhstat, mostly
    through their profile-dependent threshold, so the prescription sets the sign of the gap.
  * pbhstat's `NonLinearModel` was not rerun in our conventions: it has its own threshold
    (K̃ = 6, g_c(w)), and the printed formula lacks 1/Ω_CDM.
* **Threshold ±5%** (`scripts/threshold_sensitivity.py`, `scripts/curvaton_branches.py`), as in
  Iovino et al. 2406.20089:
  * A* moves by −19 to +27%, R and R/R_Gauss by ≤ 0.3% (Gaussian, f_NL = −2 to −0.5);
  * σ_P moves by up to ±0.54σ; the statistics shift stays at +0.23 to +0.75σ (reference +0.31 to
    +0.68);
  * for the curvaton r = 0.9, R moves by up to 6%. Its map is symmetric about ζ_g = −2r/3,
    where the curvaton field changes sign, so troughs of ζ_g are also peaks of ζ. The two
    branches have their own R (0.75–0.79 and 0.81–0.84), each stable to ≤ 0.6%, and the
    threshold shifts the balance: at Δ = 0.3 the troughs give 3%, 39% and 92% of f_PBH in
    threshold statistics for C_th × 0.95, 1 and 1.05. Checked: the branches add up to f_PBH to
    6×10⁻⁶, and a ζ_r window of 24σ_r instead of 16 changes nothing.
* **SMBHB background** (`scripts/pta_smbhb.py`), NANOGrav's astrophysical prior (2306.16219
  eq. A1), marginalised:
  * Bayes factors SIGW + SMBHB vs SMBHB alone: 94–105 (BPL, lognormal Δ ≤ 1) and 38 (Δ = 2),
    against NANOGrav's 23–59;
  * σ_P is capped at 2.3–2.8σ (2.0–2.5σ with log10 A ≥ −4), because the allowed region now
    contains the SMBHB-only solutions;
  * the statistics shift is +0.04 to +0.39σ (median +0.26, against +0.43 without);
  * the profile Δχ² keeps the shift (7.7 → 12.5 for the Gaussian BPL) and saturates at ≈ 18.
* **Still to do for the journal:** item 6 of the review (prose, journal class, BibTeX, vector
  figures).

## Pitfalls, confirmed or added

* Exponential sensitivity: log β throughout ✓ (a 1% shift in A moves β by ~50% at ν ≈ 7).
* USR template singularity at ζ_g = μ_*: integration domain cut there ✓.
* Peaks need the cut; apply it to both ✓.
* **New:** in the Musco et al. threshold it matters *when* the transfer function acts. At
  r̂_m/10, as Musco et al. prescribe, it is ≈ 1 (≈ no transfer function); at horizon crossing,
  the notation of the variances, it raises δ_th by up to 0.11 for broad spectra.
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
* **New (journal):** the curvaton map has a second branch, ζ_g < −2r/3, where the curvaton field
  changes sign. For f_NL < 0 it competes with the first, and quantities that differ between
  branches (such as R) then depend on the threshold through their balance.
* **New (NG SIGW):** for broad spectra the Monte Carlo of the connected terms has heavy-tailed
  weights. They come from pairs of vertices that are nearly independent: cos 2φ averages to
  zero but carries the variance. Averaging each sample over rotations about k removes them, and
  the balance heuristic over the spectral momenta keeps soft momenta under control.
