# Phase 0: is the gap open? (checked 23 Oct 2026)

## Verdict

**Yes, but it has to be stated more narrowly than in the plan.** No paper up to October 2026
computes peaks theory with local non-Gaussianity (NG) for broad spectra *on the compaction
function, as the counterpart of the Ferrante / Ianniccari threshold statistics*, and compares the
two prescriptions, let alone propagates the difference to the PTA tension. Four papers from
2024–2026 say explicitly that this is unsolved. But peaks-with-NG work at finite width does
exist in other forms (Riccardi et al. 2021, Escrivà 2026), and a general formalism for peaks of
local-NG fields appeared in 2025 (Germani et al.). The paper must cite these and say what it adds.

Sharpened claim (suggested wording):

> We construct peaks theory on the compaction function with local NG (quadratic and
> logarithmic/USR templates) for broad spectra, as the exact counterpart of threshold statistics:
> same threshold, windows, transfer function and cut. We quantify the amplitude ratio R(Δ, NG)
> between the two, test the peak-correspondence assumption on 3D realisations, and propagate
> R to the significance of the PBH-overproduction tension for the NANOGrav 15-year signal.

## How this was checked, and its limits

arXiv and INSPIRE are blocked from this environment, so the INSPIRE "cited by" scan in the plan
could not be run. Instead: several semantic searches of the alphaXiv index (2.5M arXiv papers,
covering up to Oct 2026) on peaks / threshold statistics / NG / broad spectra / PTA, plus full-text
reads of every candidate. **Run the INSPIRE scan yourself before submitting** (queries below);
it is cheap and catches papers that a semantic search misses.

```
refersto:arxiv:2402.11033 and date>2025
refersto:arxiv:2508.19328
refersto:arxiv:2509.14326
refersto:arxiv:2201.13345 and date>2025
refersto:arxiv:2109.00791 and date>2025
refersto:arxiv:2607.08738
refersto:arxiv:2503.05434
```

## Papers that state the gap explicitly

| paper | what it says |
|---|---|
| Ianniccari, Iovino, Kehagias, Perrone, Riotto, [2402.11033](https://arxiv.org/abs/2402.11033) (PRD 109, 123549; arXiv title "The broader, the better") | Footnote 2: uses threshold rather than peak statistics. Sec. V: peaks need slightly smaller amplitudes in the Gaussian case; "We leave a deeper analysis for the discrepancies between threshold statistics and peak theory in presence of primordial non-gaussianities for a future work." |
| Gouttenoire, Trifinopoulos, Vanvlasselaer, [2508.19328](https://arxiv.org/abs/2508.19328) (v2 Jun 2026) | NG peaks exist only for a monochromatic spectrum (it credits **Kitajima, Tada, Yokoyama, Yoo [2109.00791](https://arxiv.org/abs/2109.00791)**, which answers Phase 0 step 3), so it uses Gaussian peaks. Conclusion asks for "one that incorporates primordial non-Gaussianities consistently within peak theory". |
| Iovino, Perna, Riotto, Veermäe, [2406.20089](https://arxiv.org/abs/2406.20089) "Curbing PBHs with PTAs" (v3 Jun 2026) | "it is not clear how to include NGs in the computation of the abundance in [peaks theory] … we limit our analysis for the peak theory approach to the case of negligible primordial NGs." It also finds that Gaussian peaks relax the PTA bounds on f_PBH by a few orders of magnitude. **Not in the plan's reading list; add it.** |
| Pi, Sasaki, Takhistov, Wang, [2501.00295](https://arxiv.org/abs/2501.00295) (JCAP 09 (2025) 045) | Gaussian peaks for finite width. "previous works primarily focused on a monochromatic spectrum … Simply replacing the Gaussian profile of R(r) by the non-Gaussian local relation seems not appropriate when the power spectrum is broad … This remains an open avenue." |
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
