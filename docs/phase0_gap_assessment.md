# Phase 0: is the gap open? (checked 23 Oct 2026, re-checked 25 and 26 Oct)

## Verdict

**Yes, but it has to be stated more narrowly than in the plan.** No paper up to October 2026
computes peaks theory with local non-Gaussianity (NG) for broad spectra *on the compaction
function, as the counterpart of the Ferrante / Ianniccari threshold statistics*, and compares the
two prescriptions, let alone propagates the difference to the PTA tension. Five papers from
2024–2025 say explicitly that this is unsolved, and the INSPIRE scan (below) found nothing that
closes it. But peaks-with-NG work at finite width does
exist in other forms (Riccardi et al. 2021, Escrivà 2026), and a general formalism for peaks of
local-NG fields appeared in 2025 (Germani et al.). The paper must cite these and say what it adds.

Sharpened claim (suggested wording):

> We construct peaks theory on the compaction function with local NG (quadratic, logarithmic/USR
> and curvaton templates) for broad spectra, as the exact counterpart of threshold statistics:
> same threshold, windows, transfer function and cut. Unlike previous work, which imposes the
> peak condition on a Gaussian field, we count maxima of the non-Gaussian compaction itself, a
> local function of two correlated Gaussian fields, and validate the result with an independent
> Monte Carlo and with 3D realisations. We quantify the amplitude ratio R(Δ, NG) between the two
> statistics and propagate it to the significance of the PBH-overproduction tension for the
> NANOGrav 15-year signal.

(Added 23 Oct 2026, after Phase 3b: the derivation is in `docs/phase3b_derivation.tex`; the
Gaussian-peak approximation miscounts peaks by up to 18% but shifts A_* by ≤ 1.7%.)

## How this was checked, and its limits

First pass: several semantic searches of the alphaXiv index (2.5M arXiv papers, covering up to
Oct 2026) on peaks / threshold statistics / NG / broad spectra / PTA, plus full-text reads of
every candidate.

Second pass, the INSPIRE scan of the plan (run 23 Oct 2026 through the INSPIRE API; the
`refersto:arxiv:` form is not parsed by the API, so each paper was resolved to its record id):

| query | hits |
|---|---|
| `refersto:recid:2759520` (2402.11033), all years | 38 |
| `refersto:recid:2964325` (2508.19328) | 13 |
| `refersto:recid:2971120` (2509.14326, pbhstat) | 1 |
| `refersto:recid:2022804` (2201.13345), all years | 83 |
| `refersto:recid:1915684` (2109.00791), all years | 109 |
| `refersto:recid:3178750` (2607.08738) | 4 |
| `refersto:recid:2897904` (2503.05434) | 8 |
| full text "threshold statistics" + "peak theory" + "non-Gaussianity", since Jun 2025 | 12 |
| titles with "peak theory" / "peak statistics" / peaks + non-Gaussian, since 2025 | 4 |

193 unique citing records. I read the abstracts of all 2025–2026 records that mention peaks, NG
or finite width, and searched the LaTeX sources of the nine closest. **None computes peaks with local NG
for broad spectra on the compaction function, or compares it with threshold statistics.** What
the scan adds:

* Iovino, Perna, Veermäe [2512.13648](https://arxiv.org/abs/2512.13648) (Dec 2025), a fifth
  explicit statement of the gap: "there is no clear way to include NGs in the computation of the
  abundance [in peaks theory] in the case of a generic functional form for the curvature
  perturbation field … or for very large NGs. Hence, we limit our analysis to the TS approach,
  assuming that approximately similar results hold for the PT approach." That assumption is
  what this project tests.
* Carr, Iovino, Perna, Vaskonen [2601.06024](https://arxiv.org/abs/2601.06024) (review, Jan
  2026): "In the absence of strong primordial NG, the two methods differ only by a mild rescaling
  of the power spectrum amplitude." Strong NG is left open.
* Blas, Foster, Gouttenoire, Iovino [2602.12252](https://arxiv.org/abs/2602.12252): with the
  high-peak formula of 2406.20089, Gaussian peaks need A ≈ 6×10⁻³ against 10⁻² for threshold
  statistics, i.e. R ≈ 0.6–0.67. ngpbh gives 0.75–0.89 across the conventions in
  `docs/conventions.md`, so the convention spread in the literature is wider still: one more
  reason to lead with the double ratio.
* Not overlapping, but cite if relevant: Drees & Wang 2609.07401 and Allegrini et al. 2606.28296
  (NG to all orders, threshold statistics only); Mohammadi et al. 2609.07043 (numerical
  thresholds for finite width, stiff epochs); Escrivà 2609.03051 (multiple-peak collapse).

Re-run the scan just before submission.

Third pass (25 Oct 2026, before the Zenodo release). The same queries, plus the citers of
2306.17149, 2406.20089, 2501.00295, 2603.23025, 2512.13648 and 2601.06024, for records since
16 Oct (INSPIRE had indexed up to 20 Oct): six new records, none on peaks with NG
(backreaction, braneworld PBHs, birefringence, PBH clustering, 21 cm). The reference check for
the note added two older papers that the first two passes missed, because they predate the
citation trees scanned:

* Yoo, Gong, Yokoyama [1906.06790](https://arxiv.org/abs/1906.06790) (JCAP 09 (2019) 033):
  peaks of ζ_G with typical profiles and quadratic f_NL, for a monochromatic *and* an extended
  spectrum (P ∝ k³ e^{−3k²/2k₀²}). They assume that every peak of ζ_G is a peak of ζ (valid for
  3f_NL/5 ≳ −1/4) and do not compare with threshold statistics.
* Atal, Garriga, Marcos-Caballero [1905.13202](https://arxiv.org/abs/1905.13202) (JCAP 09
  (2019) 073): the origin of the USR/constant-roll log template; high-peak median profiles of
  ζ_g, k⁴ growth cut at the peak; amplitude estimates at fixed ν, no abundance comparison.

Both impose the peak condition on the Gaussian component, which the abstract of the note already
covers ("or with the peak condition imposed on a Gaussian field"). The gap stands.

Fourth pass (26 Oct 2026, just before submission), for records since 15 Oct:

* the citers of 26 seed papers: the 13 above, plus the closest works the paper compares with
  (1906.06790, 2102.04084, 1905.13202, 2008.03289, 2211.01728, 2011.03014, 1906.07135,
  2202.01028, 2105.01659, 2306.16219, 1912.07072, 2412.07709, 2602.12252);
* six INSPIRE keyword queries: titles with peak(s) theory or statistics; the full text with
  "threshold statistics", "peak theory" and "non-Gaussian"; titles on PBHs and PTAs; "compaction
  function" in the abstract; titles on induced GWs and NG; titles on PBH abundances or mass
  functions and NG. Each was checked to return hits over all years, so that a zero is not a
  parse failure;
* every INSPIRE record with "primordial black hole(s)" in the title;
* three arXiv API searches: PBHs with peaks, NG, thresholds, PTAs or the compaction function;
  induced GWs with NG; PTAs with priors, new physics or induced GWs.

Both INSPIRE and the arXiv API had indexed submissions up to 23 Oct, the latest arXiv
announcement. Twelve records came back: eight from the citation trees and four from the title
search. The keyword queries found none, and the arXiv searches found nothing new. None is on
peaks theory, threshold statistics or NG abundances:

* Dai, Zaldarriaga, Ali-Haïmoud [2610.28739](https://arxiv.org/abs/2610.28739) (20 Oct): the
  log-uniform priors on the intrinsic red noise of the pulsars inflate the amplitude of the
  background through a prior-volume effect. Priors uniform in power, or informed by the pulsar
  population, lower the amplitude by 1–2σ, into agreement with SMBHBs. They recompute the Bayes
  factor of NANOGrav's Gaussian SIGW template over SMBHBs with a factorised free-spectrum refit,
  the method we use, and it falls from 57 to 0.50–4.6 (their table II). This does not touch R,
  but it bears on the PTA side, so the paper now cites it in three places:
  * where it quotes Bayes factors (section 5.6);
  * in the discussion of the absolute tension, since it acts in the same direction as the SMBHB
    background;
  * in the limitations, since our refit inherits NANOGrav's priors on the red noise.
* Not relevant to the gap:
  * Braglia et al. 2610.26669, mean-field backreaction in single-field inflation;
  * Lorenzoni et al., PRD 114 L081306, PBH formation with a light spectator field;
  * Carr et al. 2610.28723, Poisson clustering of PBHs;
  * Hütsi 2610.23806, clustering of little red dots;
  * Jenks 2610.24067, GW and 21 cm probes;
  * Aldecoa-Tamayo 2610.26482, braneworld PBHs;
  * Minotti et al. 2610.27082, cosmic birefringence;
  * Morales-Zapien et al. 2610.30510, microlensing bounds;
  * Munoz et al. 2610.23218, gauge choices in simulations of dust collapse;
  * Singhal & Banerjee 2610.21658, bounces and magnetic fields;
  * Branco & Rosa 2610.21951, axion capture;
  * also checked, from just before the window: Yang, Wang & Li 2610.18460, the reheating
    history from PTA data.

All 67 earlier references were compared with their current INSPIRE records (journal, volume,
page, year, errata, withdrawals): none has changed. The gap stands. If the submission slips by
more than about a week, re-run the same queries (a few minutes). The citers of 2306.16219,
2406.20089 and 2601.06024, together with the title search, caught all twelve records.

## Papers that state the gap explicitly

| paper | what it says |
|---|---|
| Ianniccari, Iovino, Kehagias, Perrone, Riotto, [2402.11033](https://arxiv.org/abs/2402.11033) (PRD 109, 123549; arXiv title "The broader, the better") | Footnote 2: uses threshold rather than peak statistics. Sec. V: peaks need slightly smaller amplitudes in the Gaussian case; "We leave a deeper analysis for the discrepancies between threshold statistics and peak theory in presence of primordial non-gaussianities for a future work." |
| Gouttenoire, Trifinopoulos, Vanvlasselaer, [2508.19328](https://arxiv.org/abs/2508.19328) (v2 Jun 2026) | NG peaks exist only for a monochromatic spectrum (it credits **Kitajima, Tada, Yokoyama, Yoo [2109.00791](https://arxiv.org/abs/2109.00791)**, which answers Phase 0 step 3), so it uses Gaussian peaks. Conclusion asks for "one that incorporates primordial non-Gaussianities consistently within peak theory". |
| Iovino, Perna, Riotto, Veermäe, [2406.20089](https://arxiv.org/abs/2406.20089) "Curbing PBHs with PTAs" (v3 Jun 2026) | "it is not clear how to include NGs in the computation of the abundance in [peaks theory] … we limit our analysis for the peak theory approach to the case of negligible primordial NGs." It also finds that Gaussian peaks relax the PTA bounds on f_PBH by a few orders of magnitude. **Not in the plan's reading list; add it.** |
| Pi, Sasaki, Takhistov, Wang, [2501.00295](https://arxiv.org/abs/2501.00295) (JCAP 09 (2025) 045) | Gaussian peaks for finite width. "previous works primarily focused on a monochromatic spectrum … Simply replacing the Gaussian profile of R(r) by the non-Gaussian local relation seems not appropriate when the power spectrum is broad … This remains an open avenue." |
| Iovino, Perna, Veermäe, [2512.13648](https://arxiv.org/abs/2512.13648) (Dec 2025) | "no clear way to include NGs" in peaks theory for a generic F or large NG; uses threshold statistics, "assuming that approximately similar results hold for the PT approach". |
| Cole, Fumagalli, [2509.14326](https://arxiv.org/abs/2509.14326) (pbhstat) | Implementation limited to nearly Gaussian initial conditions; NG listed as future work. |
| Kushwaha, Suyama, [2603.23025](https://arxiv.org/abs/2603.23025) (v2 16 Oct 2026) | Gaussian PS vs peaks for narrow and broad spectra; they agree for narrow and **differ for broad spectra**. No NG. |
| Inui et al., [2411.07647](https://arxiv.org/abs/2411.07647) | Peaks with logarithmic NG and numerical thresholds, monochromatic only; PTA overproduction "serious for most of the parameter space"; notes that broad spectra could relax it. |

## Closest existing work: cite it and say how you differ

| paper | overlap | difference from this project |
|---|---|---|
| **Escrivà, [2607.08738](https://arxiv.org/abs/2607.08738) (Jul 2026)** | BBKS peaks of ζ_G + logarithmic NG + finite-width spectra + numerical thresholds. Finds that for negative NG and broad spectra, rare *profile deformations* raise the abundance by orders of magnitude (Q_A ≈ 0.62 for β_NG = −2). | Peaks with a one-parameter shape dispersion around a reference profile. No threshold-statistics comparison, no quadratic f_NL, no PTA posterior. **This is the biggest risk to the novelty claim and also a systematic that may exceed R−1 for negative NG.** Discuss it, and ideally reproduce one of its numbers. |
| Germani, Gorji, Uwabo-Niibo, Yamaguchi, [2503.05434](https://arxiv.org/abs/2503.05434) | General statistics of peaks of local-NG fields F = F[F_G], valid at any NG strength (applied to sphericity). | Not applied to the PBH abundance, the compaction function, broad spectra or PTA. Use it as the formal starting point for the Phase 3 derivation. |
| Riccardi, Taoso, Urbano, [2102.04084](https://arxiv.org/abs/2102.04084) | Peaks with quadratic local NG for *finite-width* (lognormal, broken power-law) spectra. | Peaks of the curvature field with a linear δ–ζ threshold, not the compaction function; compared with older exponentiated-cumulant estimates, not modern threshold statistics. The plan's summary ("NG effect on peaks smaller than threshold statistics suggests") should read: earlier cumulant-expansion estimates overstated the NG effect; a ≲ factor-2 change in P_R absorbs it. |
| Young, [2201.13345](https://arxiv.org/abs/2201.13345) | Peaks of the full non-linear compaction with f_NL, g_NL. | Delta-function spectrum only. **App. B shows that for broad spectra high peaks of C do not correspond to high peaks of ζ** (γ_Cζ → 0). This is the central obstacle for Phase 3; cite it as the source of the "assumption to test". |
| Kitajima, Tada, Yokoyama, Yoo, [2109.00791](https://arxiv.org/abs/2109.00791) | Peaks of −Δζ with typical profiles, local NG (f_NL, exponential tail); formalism written for any spectrum. | Numbers only for monochromatic spectra; the broad-spectrum f_NL case was announced "in preparation" in 2021 and no follow-up appears in the index. App. A: NG enhances or suppresses PBHs **more** in peaks than in Press–Schechter (monochromatic), unlike Young and Riccardi. The literature disagrees on the direction, which motivates the dated prediction in Phase 0 step 5. |
| Atal, Cid, Escrivà, Garriga, [1908.11357](https://arxiv.org/abs/1908.11357) | Shape dispersion + NG in single-field inflation. | High-peak median profile only; no statistics comparison. |

Active groups (who could close the gap first): Riotto–Iovino–Kehagias (explicitly deferred it),
Yoo–Tada–Yokoyama–Kitajima (in preparation since 2021), Escrivà, Germani–Gorji–Yamaguchi,
Pi–Sasaki–Wang. Expect competition on a 6–12 month horizon. A short letter (Phase 2, the
Gaussian-peak R map, and the Gaussian PTA number) before the full paper protects priority.

## Corrections to the plan PDF

### References (all opened and checked)

| plan says | status | correction |
|---|---|---|
| 2402.11033 "The PBH abundance: the broader, the better" | ✓ arXiv title | Published as "Primordial black hole abundance: the importance of broadness", PRD 109, 123549 (2024); cite the published version. Footnote 2, Sec. IV variables x₀…x₃, Sec. V variances σ_a, σ_b, σ_c, K ≈ 3.3, γ ≈ 0.36 all confirmed. |
| 2508.19328 notes NG peaks exist only for monochromatic spectra | ✓ | The credited reference is Kitajima et al. 2109.00791. |
| 2509.14326 pbhstat, nearly Gaussian only | ✓ | But see the normalisation issues below: code v0.0.6 ≠ paper. |
| 2103.01056 for the peak-correspondence requirement | ✓ (as a side remark) | Atal & Domènech mention it in one sentence (true for the log template, not its quadratic truncation). Use Young 2201.13345 App. B, Kitajima et al. §2.2 and Germani et al. 2503.05434 as the primary references. |
| 2103.01056 for "NG corrections to SIGW small if PBHs not overproduced" | ✓ | They show A_R F_NL² ≪ 1 whenever PBHs are not overproduced (for the single-field template). |
| Riccardi–Taoso–Urbano: NG effect on peaks smaller than threshold statistics suggests | ~ | See table above: they compare with exponentiated-cumulant estimates, and use finite-width spectra. |
| 2306.17149 "overproduction disfavours Gaussian SIGW at about 2σ" | ✓ | The abstract: 2σ for Gaussian, 3σ for single-field with NG (NANOGrav15). The statement is that the f_PBH = 1 band lies outside the 2σ region of the (k*, A) posterior for a BPL (α = 4, β = 3, γ = 1). They fit the 14 NANOGrav free-spectrum bins directly (not PTArcade) and use QCD-dependent K, γ, C_th, Φ. |
| 2501.11320 "SIGWAY" | ✓ | It is the LISA CosWG paper "Reconstructing primordial curvature perturbations via SIGWs with LISA" (JCAP 05 (2025) 062), which introduced the SIGWAY code. |
| 2603.23025 "Revisiting power-spectrum constraints from PBH abundances" | ✓ | Exact title: "Revisiting Constraints on Primordial Curvature Power Spectrum from PBH Abundances" (Kushwaha & Suyama). |
| 2306.16219 NANOGrav new-physics search | ✓ | Afzal et al., ApJL 951 L11. |
| 2306.16377 PTArcade | not opened | Verify before citing (as the plan already says). |

### Method (the ones that change results)

1. **"Peaks theory" must be defined before Phase 3.** Three inequivalent Gaussian prescriptions
   exist for broad spectra: (i) BBKS peaks of the linear compaction C_l(x; r) at fixed r
   (pbhstat `PeaksTheoryModel`, Young 2022, 2406.20089, 2508.19328); (ii) peaks of ζ or −Δζ with
   typical profiles (Yoo et al. 2008.02425, Kitajima et al., Pi et al.); (iii) the Germani–Sheth /
   Fumagalli–Garriga–Germani–Sheth non-linear compaction statistics with a peak condition in
   position and a maximum in r (pbhstat `NonLinearModel`). R depends on the choice. This repo
   uses (i), because it differs from threshold statistics *only* by the peak constraint, which is
   what Phase 4 needs; (iii) is the natural robustness check.
2. **Normalisation conventions move R by ~10%, the same size as the effect being measured.**
   pbhstat v0.0.6 has no factor 2 in Press–Schechter (its paper's eq. 7 has one), and its peaks
   integrate n(ν) over dg without dν/dg = 1/σ₀ and with an extra 1/2 (paper eq. 17). With those
   two conventions copied, ngpbh reproduces pbhstat to 0.017% (`results/phase1_pbhstat_comparison_pbhstat.csv`);
   with standard BBKS normalisation, peaks need ~10% smaller A than pbhstat's peaks. Fix the PS
   factor, the peak volume (4πr³/3 vs (2π)^{3/2} r³) and the Jacobian, and state them. Tell
   Cole and Fumagalli when you email them (Phase 0 step 4).
3. **The 1% match with pbhstat is a numerics check only**, in pbhstat's conventions (linear
   estimator, g_c = 0.77, no transfer function, peak-normalised lognormal). The physics runs use
   the Ferrante conventions (non-linear C, transfer function, shape-dependent threshold), which
   pbhstat's PS and peaks models do not implement.
4. **The threshold depends strongly on whether the transfer function enters the mean profile.**
   Musco's prescription gives C_th = 0.53 (with T) vs 0.42 (without) at Δ = 2, against
   0.593 vs 0.593 at Δ = 0.1. And 2402.11033 argues for 2/5 for broad spectra. That is a bigger
   lever than the statistics; fix one prescription for both methods (Phase 4 already requires
   this) and show the others as a systematic.
5. **Amplitude normalisation.** The plan's lognormal is normalised to ∫P dln k = A; pbhstat's
   is normalised to the peak value. They differ by √(2π)Δ (2.5 at Δ = 1).
6. **Phase 5 exit criterion.** "Gaussian case gives 2σ" should become: reproduce the position of
   the f_PBH = 1 band of 2306.17149 Fig. 3 (Gaussian, non-linearities only, BPL α = 4, β = 3,
   γ = 1) in the (k*, A) plane to a few percent in A, then compute the significance the same way
   they did (credible-region level), before switching to posterior-mass σ.
7. **Profile dispersion (Escrivà 2026) is a competing systematic** for negative NG and broad
   spectra. Either include it as a robustness section or scope it out explicitly.

## Dated prediction (Phase 0 step 5), to be filled in by you

Write it down before looking at `results/phase2_amplitudes.csv`. For reference, the literature
disagrees: Kitajima et al. (monochromatic) find NG effects larger in peaks than in threshold
statistics; Young (monochromatic, compaction) and Riccardi et al. (finite width, curvature
peaks) find NG effects smaller than earlier estimates.
