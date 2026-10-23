# Non-Gaussian PBH abundance: threshold statistics vs peaks theory

Code for Option A of the project plan: does the gap between threshold statistics and peaks
theory, once local non-Gaussianity (NG) is included for broad spectra, change how strongly
PBH overproduction disfavours a scalar-induced GW explanation of the PTA signal?

* `docs/phase0_gap_assessment.md` — literature check of the gap (Oct 2026) and corrections to the plan.
* `docs/revised_plan.md` — the plan with those corrections folded in.
* `docs/conventions.md` — every formula the code implements, with sources.

## Package `ngpbh`

| module | content |
|---|---|
| `spectra` | lognormal (integrated or peak normalisation), broken power law, flat band, Dirac delta |
| `windows`, `moments` | top-hat / shell / Gaussian windows, radiation transfer function, the variances of C_g and zeta_r and the spectral moments needed for peaks |
| `threshold` | delta_c(q), Musco et al. shape-dependent threshold |
| `nongaussian` | local templates zeta = F(zeta_g): quadratic f_NL, logarithmic (ultra-slow-roll), curvaton (sudden decay) |
| `stats` | log beta at one scale: threshold statistics and peaks (Gaussian-peak approximation) on a shared integrator |
| `abundance` | f_PBH = int dln M_H (M_eq/M_H)^(1/2) beta / Omega_DM, and the amplitude A_* for a target f_PBH |

```bash
python -m venv .venv && .venv/bin/pip install -e ".[dev]"
.venv/bin/python -m pytest            # 47 tests, ~35 s
.venv/bin/pip install pbhstat && .venv/bin/python scripts/compare_pbhstat.py   # Phase 1 exit check
```

| script | output in `results/` |
|---|---|
| `compare_pbhstat.py` | `phase1_pbhstat_comparison_<mode>.csv`: A_* from ngpbh and pbhstat, PS and peaks |
| `convention_spread.py` | `phase1_convention_spread.csv`: Gaussian R under four normalisation conventions |
| `scan_amplitudes.py` | `phase2_amplitudes.csv`: A_* over width × NG template × statistic (~5 min on 4 cores) |
| `make_figures.py` | `fig_phase2_threshold_amplitudes.png`, `fig_phase3a_ratio_preview.png` from the CSVs |

Peaks are counted in the Gaussian-peak approximation (maxima of the Gaussian field C_g); the
R values in the preview figure are therefore preliminary until Phase 3b (see `docs/revised_plan.md`).

```python
from ngpbh import AbundanceModel, CollapseModel, LogNormal, MuscoThreshold, Quadratic

spec = LogNormal(A=1.0, k_star=1e7, width=0.5)            # amplitude is solved for
C_th = MuscoThreshold(transfer=True).evaluate(spec).C_th
model = AbundanceModel(spec, CollapseModel(C_th=C_th), ng=Quadratic(1.0), statistic="threshold")
A_star = model.solve_amplitude(f_target=1.0)
```
