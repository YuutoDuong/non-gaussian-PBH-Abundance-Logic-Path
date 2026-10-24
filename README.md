# Non-Gaussian PBH abundance: threshold statistics vs peaks theory

Code for Option A of the project plan: does the gap between threshold statistics and peaks
theory, once local non-Gaussianity (NG) is included for broad spectra, change how strongly
PBH overproduction disfavours a scalar-induced GW explanation of the PTA signal?

* `docs/phase0_gap_assessment.md` — literature check of the gap (Oct 2026, incl. the INSPIRE
  citation scan) and corrections to the plan.
* `docs/revised_plan.md` — the plan with those corrections folded in, and the status per phase.
* `docs/phase3b_derivation.tex` (and `.pdf`) — peaks of the non-Gaussian compaction function.
* `docs/phase5_pta.md` — Phase 5: the NANOGrav 15 tension per statistic and NG template.
* `docs/conventions.md` — every formula the code implements, with sources.

## Package `ngpbh`

| module | content |
|---|---|
| `spectra` | lognormal (integrated or peak normalisation), broken power law, flat band, Dirac delta |
| `windows`, `moments` | top-hat / shell / Gaussian windows, radiation transfer function, the variances of C_g and zeta_r and the spectral moments needed for peaks |
| `threshold` | delta_c(q), Musco et al. shape-dependent threshold |
| `nongaussian` | local templates zeta = F(zeta_g): quadratic f_NL, logarithmic (ultra-slow-roll), curvaton (sudden decay) |
| `stats` | log beta at one scale on a shared integrator: threshold statistics; peaks as maxima of the Gaussian C_g (Phase 3a) or of the non-Gaussian C_l = F'(zeta_r) C_g (Phase 3b) |
| `abundance` | f_PBH = int dln M_H w(M_H) beta, the amplitude A_* for a target f_PBH, and A_*(k*) for many k* at once |
| `cosmology`, `thermal` | horizon mass and f_PBH weight in four conventions; Standard-Model g*(T), g*s(T) (Saikawa–Shirai) |
| `sigw` | scalar-induced GWs in radiation domination and their present-day spectrum h^2 Omega_GW(f) |
| `pta` | NANOGrav 15-yr free-spectrum likelihood (ceffyl refit) and HPD credible levels |
| `validation` | independent Monte Carlo of the Phase 3b peak density (20 Gaussian variables, exact Hessian) |

```bash
python -m venv .venv && .venv/bin/pip install -e ".[dev]"
.venv/bin/python scripts/fetch_ng15.py # NANOGrav 15-yr free spectrum (Zenodo, MD5-checked) into data/
.venv/bin/python -m pytest            # 73 tests, ~40 s
.venv/bin/pip install pbhstat && .venv/bin/python scripts/compare_pbhstat.py   # Phase 1 exit check
```

| script | output in `results/` |
|---|---|
| `compare_pbhstat.py` | `phase1_pbhstat_comparison_<mode>.csv`: A_* from ngpbh and pbhstat, PS and peaks |
| `convention_spread.py` | `phase1_convention_spread.csv`: Gaussian R under four normalisation conventions |
| `double_ratio_conventions.py` | (printed) R(NG)/R(Gaussian) under three conventions: they cancel to < 0.01 |
| `scan_amplitudes.py` | `phase2_amplitudes.csv`: A_* over width × NG template × {threshold, peaks 3a, peaks 3b} (~20 min on 4 cores) |
| `lattice_peaks.py` | `phase3b_lattice.csv`: maxima of C_l counted on 256³ Gaussian random fields vs the 3a and 3b predictions |
| `pta_validate.py` | `phase5_validation.md`: the NG15 likelihood against NANOGrav's power-law fits and 2306.17149's SIGW fit |
| `digitize_2306_17149.py` | `digitized_2306.17149_fig3.csv`: their f_PBH bands and NG15 contours, read from the vector figure (needs PyMuPDF) |
| `pta_tension.py` | `phase5_lines.csv`, `phase5_tension.csv`, `phase5_summary.md`, `phase5_compare.md`: overproduction lines per statistic, NG template and shape, and the NG15 exclusion significance (~1 h on 4 cores) |
| `pta_bpl_marginal.py` | `phase5_posterior_bpl_marginal.npz`: the NG15 posterior with the BPL shape marginalised (as 2306.17149 Fig. 3) |
| `make_figures.py`, `make_phase5_figures.py` | the figures in `results/` from the CSVs |

```python
from ngpbh import AbundanceModel, CollapseModel, LogNormal, MuscoThreshold, PeaksConfig, Quadratic

spec = LogNormal(A=1.0, k_star=1e7, width=0.5)            # amplitude is solved for
C_th = MuscoThreshold(transfer=True).evaluate(spec).C_th
model = AbundanceModel(spec, CollapseModel(C_th=C_th), ng=Quadratic(1.0), statistic="threshold")
A_threshold = model.solve_amplitude(f_target=1.0)
peaks = AbundanceModel(spec, CollapseModel(C_th=C_th), ng=Quadratic(1.0), statistic="peaks",
                       peaks=PeaksConfig(constraint="C_l"))          # Phase 3b
A_peaks = peaks.solve_amplitude(f_target=1.0)
```

Peaks need the k < 4.49/r cut (`MomentsConfig(x_cut=4.49)`) for broad spectra; use it for both
statistics when comparing them.
