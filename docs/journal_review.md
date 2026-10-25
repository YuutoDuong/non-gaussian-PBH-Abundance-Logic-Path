# Journal-readiness review: "Peaks or thresholds?"

State reviewed: `docs/note/ngpbh_note.pdf` of 25 Oct 2026 (19 pages, then 61 references), code
`ngpbh` v0.6.0 (doi:10.5281/zenodo.23643762, 96 tests). Written the way a referee for
Phys. Rev. D or JCAP would read it.

This is a self-review: the reviewer (Claude) also wrote much of the paper and the code. Its
blind spots are the author's. An independent expert read remains the most valuable step before
journal submission.

## Verdict

* **Zenodo and arXiv: ready.**
* **Journal:** the science is sound and unusually well validated. As it stands, expect "minor to
  moderate revision". Two things decide how smoothly it goes:
  * the presentation, which reads like a technical report;
  * whether the paper pre-empts four predictable requests: an astrophysical foreground in the
    PTA fit, another peak prescription in the Gaussian limit, numbers for strong negative NG
    (profile dispersion, threshold), and explicit contact with Iovino et al. 2024.

## What the paper does

It builds peaks theory for the non-Gaussian compaction function with broad spectra and local NG.
This is the exact counterpart of threshold statistics: same threshold, window, transfer function,
cut, K and γ. Maxima of C_ℓ = F'(ζ_r) C_g are maxima of the Gaussian field C_g + λζ_r, times a
rank-one Hessian factor computed to first order. A Monte Carlo and 256³ lattices confirm it.

Results:

* The amplitude ratio R = A*_peaks / A*_threshold is 0.74–0.82 for |f_NL| ≤ 2, and NG cancels in
  R/R_Gauss to about 4%.
* For NANOGrav 15, peaks strengthen the exclusion of the SIGW interpretation by 0.2–0.8σ for every
  NG template.
* NG in the induced GWs (all seven terms) lowers the significance by up to 0.26σ and leaves the
  statistics shift unchanged.
* Three conventions are quantified: the timing of the transfer function in the threshold, the
  f_PBH normalisation, and the frequency range.

## Strengths

1. **A precise question and a like-for-like comparison.** Every analysis choice is the same for
   both statistics, so R isolates the peak condition.
2. **A real technical result.** The Kac–Rice density for maxima of the non-Gaussian compaction
   reduces to BBKS for a shifted Gaussian field plus a rank-one term. It is checked by a 20-variable
   Monte Carlo (3×10⁻⁴) and by lattice counts (1%). It also shows that the Gaussian-peak
   approximation miscounts peaks by up to 18%.
3. **Validation against independent results**, which referees rarely see in this field:
   * pbhstat, to 0.017% in A*;
   * the benchmark of Ferrante et al.;
   * the bands of 2306.17149 (Fig. 3) and its threshold (Fig. S3), read from vector paths;
   * all six NG-SIGW terms of Adshead et al.;
   * the NANOGrav amplitude and slope;
   * Kohri & Terada's 0.8222.
4. **Useful clarifications.** When the transfer function acts in the threshold moves A* by 18–62%;
   the f_PBH constant in common use is 3–13 times too large; the frequency range moves σ by
   −1.5 to +1σ.
5. **Reproducibility.** An archived release, a script per table and figure, and tests.

## Main concerns and recommended actions

Effort: S = hours, M = a day or so of work plus compute, L = more.

### 1. The headline effect is modest: lead with the factorisation (S)

* **The numbers:**
  * true NG peaks against the Gaussian-peak approximation: ≤ 1.7% in A*, ≤ 0.04σ;
  * the statistics shift: +0.2 to +0.8σ;
  * NG in the GWs: ≤ 0.26σ.
* **The risk:** a referee may conclude that the new machinery does not change the answer.
* **Action:** state the positive result plainly. The conversion from threshold statistics to
  peaks factorises: A*_peaks(Δ, NG) ≈ R_Gauss(Δ) × A*_threshold(Δ, NG), to about 4% for
  |f_NL| ≤ 2. Give R_Gauss(Δ), and R for the BPL, as a table or fitting formula, so readers can
  convert any threshold-statistics result. That is what people will use and cite.

### 2. Which peaks? (M)

* **The issue:** R is defined for BBKS maxima of C at a fixed scale. Two other prescriptions define
  the abundance differently:
  * typical-profile peaks of ζ (Yoo et al.; Kitajima et al.; Pi et al.);
  * the non-linear statistics of Germani & Sheth.

  The limitation is now stated, but without a number.
* **Action:**
  * compute R_Gauss for at least one alternative in the Gaussian limit; pbhstat's
    `NonLinearModel` is already installed from Phase 1;
  * compare our R_Gauss = 0.77–0.79 with earlier Gaussian comparisons: Gow et al. 2021,
    Ianniccari et al. Sec. V, Kushwaha & Suyama, and the high-peak formula (R ≈ 0.6–0.67).
* **Also:** add a short paragraph on how the construction differs from each closest work:
  * Yoo, Gong & Yokoyama count peaks of ζ_G with typical profiles;
  * Riccardi et al. count curvature peaks with a linear threshold;
  * Escrivà adds shape dispersion;
  * Germani et al. treat sphericity.

  The table in `docs/phase0_gap_assessment.md` already has the content.

### 3. Strong negative NG: profile dispersion and the threshold (S–M)

* **The issue:**
  * Escrivà (2607.08738) finds that profile deformations raise the abundance by orders of
    magnitude for negative NG and broad spectra;
  * NG also shifts C_th (Kehagias, Musco & Riotto 2019; Escrivà et al. 2022);
  * the paper uses the mean profile and the same threshold for every template.

  Both choices act on the two statistics alike, so R is safe. The absolute significances for
  f_NL < 0 and the curvaton (the cases that relieve the tension) are not.
* **Action:**
  * vary C_th by ±5% for those templates, as Iovino et al. 2024 do, and report the change in R
    and σ_P;
  * add a paragraph estimating what Escrivà's result would do to the f_NL < 0 rows.

### 4. No astrophysical foreground in the PTA fit (M)

* **The issue:** the fit uses the NANOGrav 15 free spectrum alone, without an SMBHB component. The
  statistics shift is the point, but even Δσ depends on the shape of the posterior near the
  allowed region, and a foreground changes that shape.
* **Action:** add an SMBHB power law with the astrophysical prior of the NANOGrav new-physics
  search (amplitude and slope around γ = 13/3), marginalise over it, and report σ_P and Δσ with
  and without it. The likelihood code already exists.
* **Data currency:** NANOGrav 15 is still the most constraining public free spectrum that I could
  find. IPTA DR3 is described but I found no public GWB free spectrum; the InPTA DR2 search
  (2608.02808) appeared in Aug 2026. Say so; an EPTA DR2 cross-check is optional.

### 5. Make contact with Iovino et al. 2024 (S) — done

* **The issue:** "Curbing PBHs with PTAs" compares threshold statistics and Gaussian peaks for
  NANOGrav 15. It finds that peaks relax the PTA constraints on f_PBH by a few orders of
  magnitude. This is the same physics as ours, R < 1: at fixed A, peaks give many more PBHs. So the f_PBH allowed
  at the PTA bound on A grows, while the f_PBH = 1 line moves down and the overproduction tension
  strengthens. Readers will meet opposite words for one effect.
* **Action:** add one sentence reconciling the two framings, and one on their high-peak formula
  against our normalisation.
* **Done (25 Oct 2026):** both sentences are in the Discussion.

### 6. Presentation (M)

* **The issue:** the paper reads like a technical report:
  * many bullet lists and terse sentences;
  * the arXiv number "2306.17149" used as a name;
  * "pitfalls" aimed at named papers;
  * project jargon in the text and figures: "A*^3b/A*^3a", "Phase 3b: maxima of C_ℓ" as a figure
    title, "peaks 3b" in legends.
* **Actions:**
  * rewrite the introduction, Secs. 4–6 and the discussion in prose; keep lists for genuine
    enumerations;
  * write "Ref. [12]" or "Franciolini et al." (still to do), and "conventions that matter" instead
    of "pitfalls" (done);
  * rename 3a/3b in the text and figures (done: A*^(C_ℓ)/A*^(C_g), "peaks of C_ℓ",
    "Gaussian-peak approx.");
  * add one conventions table: window, transfer function, cut, K and γ, threshold, f_PBH
    normalisation, bins and priors;
  * make the figures vector PDF with larger fonts, consistent styles and self-contained captions;
  * shorten the abstract to about 200 words; say "paper", not "note";
  * use the journal class (revtex4-2 for PRD, jcappub for JCAP) and a BibTeX bibliography with
    DOIs, generated from INSPIRE.

### 7. Appendix B will be checked first (S) — done

* **The issue:** it claims that eq. (14) of 2306.17149 overstates f_PBH by 12.7, and the M_eq form
  by 3–5. It criticises a prominent group's paper, so a referee will check it.
* **Action:**
  * define β precisely: the energy fraction, including M/M_H;
  * show the numerical cross-check explicitly: Carr et al. 2010, f = 4.11×10⁸ β′(M/M⊙)^{−1/2},
    against our 1.0×10⁻⁸ with the same Ω_CDM h²;
  * say in one line where the larger constants come from, if that can be pinned down.
* **Done (25 Oct 2026):**
  * β is defined as the energy fraction weighted by M/M_H;
  * the result agrees with Carr et al. 2010 to 2% and with the exact eq. II.3 of Carr & Kühnel 2020
    to 1.5%;
  * 7.9×10⁻¹⁰ = [2.4 (M_eq/M⊙)^{1/2}]⁻¹ is the normalisation of the approximation
    f_PBH ≈ 2.4 β (M_eq/M)^{1/2} of Carr & Kühnel, which is already a dark-matter fraction;
    dividing it by Ω_DM again gives 3.8, and the approximation is 3.4 times the exact result,
    together 12.7.

## Minor comments

Done on 25 Oct 2026 unless marked optional:

* Write "about 4%" instead of "to within ±4%": the double ratio spans 0.958–1.036, which is −4.2%.
* Repeat the definition of σ_P in the table captions.
* In the QCD limitation, cite the QCD-era works:
  * on the threshold: Musco, Jedamzik & Young 2023; Franciolini et al. 2022;
  * on the SIGW: Abe, Tada & Ueda 2021; Balaji, Domènech & Franciolini 2023.
* Optional:
  * a lattice check at Δ = 2 (512³) would close the last validation gap; the Monte Carlo already
    covers the formula there;
  * the full non-perturbative SIGW for the USR template (Perna et al. 2024) would replace the
    A³ truncation.
* Notation: the lognormal and BPL spectra are now written out (Sec. 2.1), with their
  normalisations, and the clash of the BPL's β, γ with the mass fraction and the critical exponent
  is flagged.

## Journal choice

* **JCAP:** the natural home. Much of the PBH-abundance and SIGW work cited here appeared there,
  and there is no page limit.
* **Phys. Rev. D:** equally suitable, with a broader readership.
* **Not PRL:** too long. A letter would need a single sharp result.

Both journals accept papers already posted on arXiv or Zenodo.

## Submission checklist

* [ ] journal class; BibTeX bibliography with DOIs
* [ ] cover letter (Claude can draft it)
* [ ] suggested referees: experts in peaks theory with NG and in PBH thresholds who are not
  authors of 2306.17149, e.g. the authors of the works cited for peaks with NG and for the
  Gaussian comparisons
* [x] data availability (Code and data) and the AI-use acknowledgement
* [ ] funding acknowledgement, if any
* [ ] INSPIRE scan re-run just before submission
* [ ] arXiv first (endorser needed), then the journal with the arXiv number
* [ ] a new release (v0.7.0) once the journal version adds computations, so the archive matches

## Suggested order of work

1. **Quick fixes (S) — done on 25 Oct 2026, before the Zenodo preprint:**
   * items 5 and 7;
   * the wording in the minor comments;
   * the 3a/3b jargon in the text and figure labels.
2. **Computations (M):**
   * item 4 (foreground);
   * item 3 (threshold ±5%);
   * item 1 (R_Gauss table or fit);
   * item 2 (an alternative prescription).
3. **Rewrite and format (M):** item 6, the BibTeX bibliography and PDF figures.
