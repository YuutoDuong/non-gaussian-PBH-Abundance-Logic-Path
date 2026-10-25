# Conventions and formulas implemented in `ngpbh`

Every choice below changes A_* at the percent level or more, which is the size of the effect
the project measures. Change them only deliberately, and identically for both statistics.

## Fields and compaction (radiation domination, Φ = 2/3)

* Metric ds² = −dt² + a² e^{2ζ} dx², ζ = F(ζ_g) with ζ_g Gaussian.
* Linear compaction C_l(r) = −2Φ r ζ′(r); compaction C = C_l − C_l²/(4Φ) (= C_l − 3C_l²/8).
* With NG: C_l = F′(ζ_g(r)) C_g, C_g = −2Φ r ζ_g′(r) (2402.11033 eq. 40).
* Type I only: C_l < 2Φ = 4/3. Threshold in C_l: C_l,th = 2Φ[1 − √(1 − C_th/Φ)].

## Spectra (`spectra.py`)

* Lognormal, integrated normalisation (default): P = A/(√(2π)Δ) exp[−ln²(k/k*)/(2Δ²)],
  so ∫P dln k = A (2306.17149, 2402.11033, 2211.01728). `normalization="peak"` gives
  pbhstat's P = A exp[…]; the two differ by √(2π)Δ.
* Broken power law: 2306.17149 eq. (1), A = peak value.
* Flat band: P = A Θ(k − k_min) Θ(k_max − k) (2211.01728 eq. 26).

## Variances (`moments.py`), x = kr, P^T = T²(x) P

| symbol | definition | meaning |
|---|---|---|
| σ_c² | (4Φ²/9) ∫dln k x⁴ W² P^T | ⟨C_g²⟩ |
| σ_cr² | (2Φ/3) ∫dln k x² W W_s P^T | ⟨C_g ζ_r⟩ |
| σ_r² | ∫dln k W_s² P^T | ⟨ζ_r²⟩ |
| s₁² | (4Φ²/9) ∫dln k x⁶ W² P^T | r²⟨\|∇C_g\|²⟩ |
| s₂² | (4Φ²/9) ∫dln k x⁸ W² P^T | r⁴⟨(∇²C_g)²⟩ |
| s_cr1² | (2Φ/3) ∫dln k x⁴ W W_s P^T | r²⟨(−∇²C_g) ζ_r⟩ |

W(x) = 3(sin x − x cos x)/x³ (top hat), W_s = sin x/x (shell), T(x) = W(x/√3) evaluated at
horizon crossing of r (τ = r). Optional hard cut W → 0 for x > 4.49 (pbhstat). The peak
moments s₂² do not converge for broad spectra without the cut (the plan's "peaks need a
cut-off" pitfall), so **runs that compare statistics apply the cut to both**; it changes the
threshold-statistics A_* by ≤ 2.3% for Δ ≤ 2. With the transfer function the moments of these
spectra also converge without it, and R hardly depends on it: removing the cut, or moving it to
3.5, 6 or 7.725 (the second zero of W), changes R by ≤ 1.2% and R/R_Gauss by ≤ 0.8% for
lognormal Δ = 0.5–2, Gaussian and f_NL = ±2 (`scripts/cut_sensitivity.py`,
`results/r_cut_sensitivity.csv`).

## Threshold (`threshold.py`)

* δ_c(q) = (4/15) e^{−1/q} q^{1−5/(2q)} / γ(5/(2q), 1/q) (Escrivà–Germani–Sheth; Musco et al.).
* Musco et al. 2011.03014 via the summary in the supplement of 2306.17149: r_m from
  ∫dln k [(k²r²−1) sin(kr)/(kr) + cos(kr)] P^T = 0; q from
  F(q)[1+F(q)] q = −½[1 + r_m ∫dk k cos(kr_m) P^T / ∫dk sin(kr_m) P^T], F = √(1 − 3δ_c/2).
* Checks: Dirac spectrum gives k*r_m = 2.7437, δ_c = 0.5947. Lognormal with transfer function:
  δ_c = 0.593, 0.578, 0.559, 0.537, 0.529 for Δ = 0.1, 0.3, 0.5, 1, 2; without it the broad end
  drops to 0.458 (Δ = 1) and 0.419 (Δ = 2).
* **When the transfer function acts.** Musco et al. smooth the spectrum with T(k, η) on
  superhorizon scales, "when r̂_m ≃ 10 r_H or larger" (2011.03014 Sec. IV), i.e. η ≈ r̂_m/10, where
  T ≈ 1 except for modes with k r̂_m ≳ 20. `MuscoThreshold(transfer=True, eta_over_rm=0.1)` does
  this. For every Phase 5 spectrum it agrees with no transfer function (`transfer=False`) to
  ≤ 0.0005 in δ_th: BPL(4,3,1) 0.5438 vs 0.5441; lognormal 0.5301/0.5304 (Δ = 0.5), 0.4577/0.4581
  (Δ = 1) and 0.4187/0.4192 (Δ = 2). `transfer=True` with the default `eta_over_rm=1` evaluates T
  at horizon crossing, as for the variances. This strong smoothing gives 0.5700, 0.5593, 0.5371
  and 0.5288 for the same spectra.
* **2306.17149 Fig. S3** (δ_th against the lognormal width and the BPL slope β) is reproduced over
  its whole range by Musco et al.'s smoothing or by no transfer function: to ≤ 0.0005 for
  0.1 ≤ Δ ≤ 1, and to ≤ 0.0016 for BPL(4, β, 1) with 1.9 ≤ β ≤ 3. The horizon-crossing version is
  off by up to 0.08 (`scripts/check_threshold_2306_17149.py`, `results/threshold_2306.17149_check.md`).
  Their eqs. (S4)–(S5) write T(k, r_m), the notation of their horizon-crossing variances, but the
  computation follows Musco et al. For β < 2 the integrals without T converge only conditionally;
  `threshold.py` handles the tail with QUADPACK's QAWF, and their curve lies up to 0.007 below ours
  there.
* Phase 5 uses Musco et al.'s threshold as its baseline (variant `base`, `transfer=False`) and the
  horizon-crossing threshold as a variant (`T`). Phases 1–4 used the horizon-crossing threshold.
  For BPL(4,3,1) it raises A_* by 0.071 dex (18%), for every statistic:
  C_l,th = 2Φ[1 − √(1 − C_th/Φ)] grows faster than C_th (0.761 → 0.826). The ratio R between the
  statistics hardly depends on the choice: it changes by −0.1% to +0.5% for the quadratic templates,
  and up to +6% for the curvaton and strong USR tails (44 cases, Phase 5 lines at k* = 10⁸ Mpc⁻¹).

## Collapse

M/M_H = K (C − C_th)^γ with K = 3.3, γ = 0.36 (lognormal, 2402.11033 and 2211.01728);
K = 4.36 for the flat band. pbhstat's PS and peaks models use the linear estimator,
M/M_H = K (C_l − g_c)^γ with g_c = 0.77, K = 4 (`estimator="linear"`).

## Mass fraction at one scale (`stats.py`)

    β = ∫ dζ_r ∫_{C_l,th}^{2Φ} dC_l |F′(ζ_r)|⁻¹ P_G(C_l/F′, ζ_r) K (C − C_th)^γ w

with P_G the bivariate Gaussian of (C_g, ζ_r). Threshold statistics: w = 1 (2211.01728 eq. 56).
Peaks (Gaussian-peak approximation): maxima of the Gaussian field C_g(x; r) in the volume
V = 4πr³/3,

    w = V √(2π) G(m, s) / ((2π)² R_*³),   G(m, s) = ∫₀^∞ f(x) N(x; m, s²) dx,

where x = −∇²C_g/s₂, m and s² are the conditional mean and variance of x given
(ν = C_g/σ_c, z = ζ_r/σ_r), R_* = √3 s₁/s₂ · r and f is BBKS (A15). For F′ < 0 the maxima of
C_l are minima of C_g and G(−m, s) is used. With Gaussian ζ this is the full BBKS number
density (no high-peak approximation); `PeaksConfig(density="high_peak")` replaces G by (γν)³.

What the approximation leaves out: with F″ ≠ 0, ∇C_l = F′∇C_g + F″C_g∇ζ_r, so maxima of C_l are
displaced from maxima of C_g. `PeaksConfig(constraint="C_l")` counts the true maxima of C_l
(Phase 3b, `docs/phase3b_derivation.tex`): at each (C_l, ζ_r) node they are BBKS maxima of the
Gaussian field Y = C_g + λζ_r, λ = C_l F″/F′², so the same w with R_*, m and s of Y, built from
the moments of C_g and ζ_r up to fourth order (the three extra ones are in the moments table
below), times the rank-one factor Q = 1 − τR₂:

    τ = sign(F′) κ v / σ₂(Y),  κ = −c F′ (1/F′)″,  v = (s₁² s₁r² − s_cr1⁴) / (3 s₁(Y)²) / r²,
    R₂ = 27 ∫f₂ N / ∫f N,       f₂(x) = E[e₂(−ζ) Θ | x]   (→ x²/3 − 1/3; 9/m for high peaks)

κ = 0 for the logarithmic template; κ = −2λ²/c for the quadratic one. Q is exact to first order
in τ (checked against an independent Monte Carlo of the 20 Gaussian variables with the exact
Hessian, `ngpbh.validation`, to 3×10⁻⁴; 0.85% at the most extreme node). Maxima of C_l counted
on 256³ Gaussian random fields follow the 3b density to ~1% (beyond the Gaussian grid trend),
while 3a misses them by up to 18% (`scripts/lattice_peaks.py`, `results/phase3b_lattice.csv`).

| extra moment (Phase 3b) | definition | meaning |
|---|---|---|
| s1r² | ∫dln k x² W_s² P^T | r²⟨\|∇ζ_r\|²⟩ |
| s2cr² | (2Φ/3) ∫dln k x⁶ W W_s P^T | r⁴⟨∇²C_g ∇²ζ_r⟩ |
| s2r² | ∫dln k x⁴ W_s² P^T | r⁴⟨(∇²ζ_r)²⟩ |

Numerics: C_l on Gauss–Legendre nodes with C_l − C_l,th ∝ t³ (absorbs the (C − C_th)^γ cusp);
ζ_r on a trapezoid grid with step ≤ √(1−γ_cr²)/12 inside a window found by a coarse pass;
everything in log space. Agreement with scipy `dblquad` to 2×10⁻⁶ in log β, with direct
Monte Carlo sampling, and stable to 10⁻⁵ when all grids are doubled.

## Abundance (`abundance.py`, `cosmology.py`)

    f_PBH = (1/Ω_DM) ∫ dln M_H (M_eq/M_H)^{1/2} β(M_H)        (2211.01728 eq. 61)

M_eq = 2.8×10¹⁷ M_⊙, Ω_DM = 0.12/0.674². Scale–mass relation (default, 2211.01728 eq. 39):
M_H = 17 M_⊙ (g*/10.75)^{−1/6} (k_H/10⁶ Mpc⁻¹)^{−2}, k_H = 1/r, g* = 106.75 (no QCD
corrections yet). pbhstat preset: R = (M_H/M_eq)^{1/2}/k_eq, Ω_CDM = 0.26.

Two more conventions use the Standard-Model g*(T), g*s(T) (Saikawa & Shirai 2005.03544, the
table PTArcade ships; `thermal.py`), with k_H = 1/r re-entering at temperature T:

* `"fivt23"` (2306.17149 eqs. 4, 14, 15): k = 1.5×10⁷ Mpc⁻¹ (g*/106.75)^{1/2} (g*s/106.75)^{−1/3}
  T/GeV, M_H = 4.8×10⁻² M_⊙ (g*/106.75)^{−1/2} (T/GeV)⁻², and
  f_PBH = (1/Ω_DM) ∫ dln M_H (M_H/M_⊙)^{−1/2} (g*/106.75)^{3/4} (g*s/106.75)^{−1} β/7.9×10⁻¹⁰.
* `"entropy"`: k = aH and M_H = 4π M_P²/H exactly, and n_PBH/s conserved after formation, so
  ρ_PBH,0 = (3/4) β (g*/g*s) T s₀. This is the same formula with 7.9×10⁻¹⁰ replaced by 1.0×10⁻⁸;
  it agrees with Carr et al. (0912.5297) to 2%.

For the same β, the M_eq formula gives an f_PBH 3–5 times larger than `"entropy"` (it drops Ω_m
and the g* heating of the radiation after formation) and `"fivt23"` one 12.7 times larger. Since
β depends exponentially on 1/A, these move A_* by only ln(factor)/|ln β|: −8% and −12% at
|ln β| ≈ 21. All conventions cancel in R and in the double ratio; they matter for the absolute
PTA tension, so Phase 5 reports them as variants.

## Scalar-induced GWs and the PTA likelihood (`sigw.py`, `pta.py`)

    Ω_GW,RD(k) = (1/6) ∫₁^∞ dt ∫₀¹ ds K(t,s) P(k(t+s)/2) P(k(t−s)/2)        (2306.17149 eqs. 6–7)
    h²Ω_GW(f) = h²Ω_r (g*/g*0)(g*s/g*s0)^{−4/3} Ω_GW,RD(k),  h²Ω_r = 4.2×10⁻⁵,  f = k/2π

Gaussian ζ, radiation domination, no QCD sound-speed effects, as in 2306.17149. Numerics: t is
split at the resonance √3 and mapped to ln|t − √3|, Gauss–Legendre panels in t and s, and the
logarithmic bracket of the kernel rewritten as y g(e) with g(e) = −e³/12 − … so that the corner
t, s → 1 has no cancellation. Checks: a flat spectrum gives 0.8222436 A² (Kohri & Terada), a
narrow spectrum approaches their Dirac closed form, scipy `dblquad` agrees to 10⁻⁸–10⁻¹², and
the k³ infrared tail matches its asymptotic form to 10⁻⁵.

PTA likelihood: the NANOGrav 15-yr HD free spectrum as KDE log-densities per frequency bin (PTArcade
data, Zenodo 10495907), the lowest 14 bins, with ρ_j² = H₁₀₀² h²Ω_GW(f_j)/(8π⁴ f_j⁵ T) and
T = 16.03 yr (the ceffyl refit). It reproduces NANOGrav's γ = 13/3 amplitude
(2.5 vs 2.4 × 10⁻¹⁵) and free-γ slope (3.07 vs 3.2), and 2306.17149's lognormal and BPL SIGW
fits (`results/phase5_validation.md`, `results/phase5_compare.md`). For tail probabilities the KDEs
are regularised (`smooth=0.05`: max of the KDE and its 0.05-dex Gaussian blur), because their
compact kernels floor the density at ln p = −36 outside the sampled range; this leaves the bulk
unchanged and moves significances by ≤ 0.1σ.

## NG corrections to the SIGW (`sigw_ng.py`)

For ζ = ζ_g + F(ζ_g² − ⟨ζ_g²⟩), F = 3f_NL/5, with Δ²_g = A s(k/k*) the spectrum of the Gaussian
component (the A of the PBH lines), the four-point function of ζ ends at F⁴ and

    Ω_GW,RD = A² T2 + F² A³ (hybrid + C + Z) + F⁴ A⁴ (reducible + planar + nonplanar)

(Adshead, Lozanov & Weiner 2105.01659 eqs. 2.33–2.39; their prefactors 2⁵, 2⁷F², 2⁸F², 2⁹F⁴ are
reproduced by our derivation). Each term is a function of x = k/k* only, so one table per spectral
shape serves every f_NL, A and k* (`results/sigw_ng_terms.csv`, `scripts/sigw_ng_terms.py`).

* **Disconnected terms** are the Gaussian functional of P_g + P_2, with the one-loop spectrum
  Δ²_2 = F²A² s2(x), s2 = x² ∫dln p (s(p)/p) ∫_{|x−p|}^{x+p} s(r) dr/r² (exact inner integral).
  Then hybrid = 2 Ω_bil[s, s2] and reducible = Ω_G[s2], both by the quadrature of `sigw.py`.
* **Connected terms** need the oscillation average of the product of two transfer functions,
  J12 = (I_A1 I_A2 + I_B1 I_B2)/2 with I = (I_A sin x + I_B cos x)/x (Kohri & Terada's late-time
  form; J11 is the Gaussian kernel), and the polarisation factor Re[z1² z̄2²], z = (q_x + i q_y)/k.
  They are computed by randomised quasi-Monte Carlo over the spectral momenta. L of the L + 1
  momenta that carry a power spectrum are drawn from the spectrum, and the last one follows from
  momentum conservation. The strategies (which momentum is fixed) are combined with the balance
  heuristic, so the weights stay bounded when a momentum becomes soft. For C and planar, each
  sample is averaged over 8 rotations about k of one spectral momentum (C: q2 at fixed q1 − q2;
  planar: q1 − q3 and q2 − q3, independently). The weights are unchanged and the estimate stays
  unbiased, but the cos 2φ oscillation of nearly independent vertex pairs cancels; this removes a
  heavy-tailed variance for broad spectra.
* **Templates beyond quadratic order**, ζ = ζ_g + Fζ_g² + Gζ_g³ + …: at O(A³) only F and G enter,
  G through P → (1 + 6Gσ²)P, i.e. +12Gσ² A² T2. The g_NL-type trispectrum vanishes under the
  azimuthal integral (Perna et al. 2403.06962 eq. 3.11). O(A⁴) would need G², F²G and FH terms,
  so the curvaton and USR maps are truncated at O(A³). The quadratic template is exact at A⁴.
  The coefficients are (F, G) = (−0.617, 0.322) for the curvaton r_dec = 0.9 (Sasaki, Valiviita &
  Wands g_NL), and (1/2μ, 1/3μ²) for the USR maps: (0.75, 0.75) for β = 3, (0.6, 0.48) for
  μ = 5/6 and (0.2, 0.053) for μ = 5/2.

Validation (`results/sigw_ng_validation.md`): the MC engine reproduces the quadrature of the
Gaussian, hybrid and reducible terms to ≲ 1%. All six NG terms reproduce Adshead et al. Fig. 4
(lognormal σ = 1/5 and 1/10, read from the figure's vector paths) within the MC error and the
digitisation, including the sign changes of Z and the negative nonplanar term.

## Normalisation conventions that differ between codes

| choice | ngpbh default | alternatives | effect on Gaussian R |
|---|---|---|---|
| Press–Schechter factor | 1 | 2 (pbhstat paper eq. 7) | R 0.788 → 0.812 (Δ = 0.1) |
| peak density | BBKS, n(ν) dν | pbhstat v0.0.6: n(ν) integrated over dg, × 1/2 | R → 0.882 |
| peak volume | 4πr³/3 | (2π)^{3/2} r³ (Young 2201.13345) | R → 0.751 |

(`results/phase1_convention_spread.csv`, pbhstat physics setup.)

The double ratio R(NG)/R(Gaussian) cancels these to < 0.01. Physics setup of the scan,
PS factor 1 or 2 and V = 4πr³/3 or (2π)^{3/2} r³ (`scripts/double_ratio_conventions.py`):

| case | Gaussian R | R/R_Gauss |
|---|---|---|
| Δ = 0.5, f_NL = +1 | 0.736–0.806 | 1.0240–1.0246 |
| Δ = 0.5, f_NL = −2 | 0.736–0.806 | 1.0034–1.0048 |
| Δ = 2, μ_* = 5/6 | 0.748–0.815 | 1.159–1.167 |

The reason: at fixed f_PBH every template sits at nearly the same |ln β|, so a constant factor
in β shifts A_* by nearly the same fraction for all of them.
