# Changelog

All notable changes are listed here. The format follows
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/), the commits follow
[Conventional Commits](https://www.conventionalcommits.org/en/v1.0.0/), and the versions follow
[Semantic Versioning](https://semver.org/spec/v2.0.0.html): the public interface is the `qgem` command
line and its JSON output.

## [Unreleased]

## [1.1.0] - 2026-11-17

### Fixed

- **obs**: the S–D mixed leptonic width had the D-wave amplitude with the wrong sign. In the Condon–Shortley
  convention of the code (⟨³D₁|S₁₂|³S₁⟩ = +√8) the amplitude is A = R_S(0) − 5R″_D(0)/(2√2 m²); the code had a
  plus sign, which is correct only for the opposite D-wave phase (⟨S|S₁₂|D⟩ = −√8). The sign was derived from the
  Dirac current expansion and checked four independent ways (`analysis/test_sd_sign.py`, CTest
  `analysis_sd_sign`; `tests/test_observables.cpp`). It changes Γ_ee of J/ψ, ψ(2S), ψ(3770), … and Υ(nS) in the
  12 physical-set models with S–D mixing in H (smeared A/B, all C/D); first-order A/B and every mass are
  unchanged. Example: Γ_ee(J/ψ), C const (s free), 5.57 → 6.34 keV; Γ_ee(ψ(3770)) 0.39 → 0.35 keV.
- **analysis**: `refresh_widths.py` recomputes the stored `widths` and `decay_benchmark` blocks of those models
  from a fixed-floor, fixed-start refit and writes them only if the refit reproduces the stored χ² and θ. All 48
  system fits matched. Gate 1, the selection and the REPORT numbers were recomputed (no model passes; the best
  compromise is still B APT). Then the same refresh was run on every other configuration with S–D mixing in the
  `main` and `floor2` sweeps (26 and 64 more system fits, all matched), and the `joint` and `joint2` fits of
  the 12 models were refit from their seeds (χ² and θ reproduced to 10⁻⁶, widths updated). The released archive of
  the first v1.1.0 commit held the `main` and `floor2` JSON *before* this refresh; it is replaced.
- **analysis**: `rerun_mcmc.py` reruns only the Markov chain of a `stats` run (same seed and steps, no bootstrap
  or LOO) and writes its widths and decay benchmark into the stored run only if the chain is the stored one
  (χ², θ, R̂, n_eff, ppp of χ² and posterior summary, to 10⁻⁶). Run on the 14 `stats` runs with S–D mixing in H:
  the 12 chains of the physical set were reproduced, so their decay posterior predictive p-values are now
  recomputed. The two contrast models, D GI run and D free run, have chains that do not converge (R̂ up to 2.5,
  and ≈ 10¹² in c c̄ for D free); their reruns are different chains, which replace the stored ones
  (`--replace-chain`, marked `"rerun": true`).
- **analysis**: `validity.mixing_angle` uses the same sign (the two roots for ψ(3770) of C const (s free) are
  −0.1° and +37.8°, previously +0.1° and −37.8°).

### Changed

- **analysis**: the decay posterior predictive p-value is reported only where the `stats` run carries the
  same widths as the main run (`tables.widths_agree`, tested in `test_model_selection.py`), so a `stats` run made
  with an older width code is shown as “–” instead of being mixed with refreshed widths.
- **docs**: every figure of the report is regenerated with ROOT 6.40.04 from the refreshed JSON, and
  `tables.md` and `data/*.csv` by `analysis/analyze.py`.
- **analysis**: the mass-and-width study (`widths` sweep) was recomputed for the 12 models with S–D mixing in H
  after the sign fix; the parameters move by at most 1.2σ, and `mass_width_fit.csv` and `tables.md` are updated.
- **docs**: the description of the leptonic amplitude in `analysis/analyze.py` and `tables.md` had the old sign
  (`R_S + 5R″_D/…`); it now reads `R_S − 5R″_D/…`.
- **analysis**: the best compromise is the model whose worst failed check is the mildest, graded in standard
  deviations (`selection.py`, `model_selection.best_compromise`), not the model with the fewest failures. The
  physical set now names B APT (2.2σ) in place of C const (s free) (4.1σ). The count rule stays as
  `best_by_count`. The winner rule is unchanged: no model passes. `selection.csv` has the new columns
  `severity` and `sev`, and `tables.md` the ranking by worst failure.
- **docs**: REPORT §5 (the worked example) and the figures of §5–§7 feature B APT, the best compromise of the
  worst-failure rule; up to v1.0.0 they featured C const (s free).
- **docs**: the |R(0)|² audit of η_c over J/ψ is stated as an NLO statement (it is proportional to K_ee/K_γγ,
  and both factors are large beyond one loop); the earlier wording that the discrepancy exceeds any plausible
  shift of K was an overclaim.

### Added

- **release**: `.zenodo.json` and an updated `CITATION.cff` for the archive at Zenodo; the fit results of the
  paper (`results/main`, `floor2`, `nested`, `widths`, `joint2`, `cv`, `stats`) as `results-archive/results-v1.1.0.tar.gz`
  with its SHA-256, and the instructions to unpack them (`results-archive/README.md`).
- **tests**: `tests/data/pdg-subset.sqlite`, a 725 kB extract of the PDG file for the fitted states (CC BY 4.0),
  and `tests/pdg_path.hpp`; the fit, run and validation tests now run in CI instead of skipping. 6j symmetry
  and triangle tests and a test of the E1 angular factor for every multiplet and its sum rule.
- **analysis**: `validity.origin_ratio_rows`, the |R(0)|² ratio of η_c(1S) to J/ψ(1S) against what the PDG widths
  require (`origin_ratio.csv`): 1.40–2.03 in all 19 models against 0.55 (1.57–3.48 before the sign fix).
- **docs**: the Godfrey–Isgur parameter cross-check (VALIDATION.md).
- **tests**: `GI.ReproducesTableIOfBarnesGodfreySwansonForCharmonium`, 31 levels of Table I of Barnes, Godfrey and
  Swanson (hep-ph/0505002, GI column; the table was read from the paper): all within 3 MeV except the 1S doublet
  (7–8 MeV below; basis converged to 0.01 MeV) and the 4S/1S tolerances noted in the test.

## [1.0.0] - 2026-11-12

First stable release: heavy-quarkonium spectroscopy with the Gaussian Expansion Method, fitted to the PDG
and compared by a full statistical analysis.

### Added

- **gem**: geometric Gaussian bases (real and complex-range), closed-form matrix elements, spin algebra
  with coupled channels, the generalized eigensolver with canonical orthogonalization, level tracking,
  variational basis optimisation and nested-size convergence (n → 2n − 1).
- **models**: A/B (Cornell, screened; constant or running couplings; spin terms none, first-order or
  smeared; optional semi-relativistic kinetic energy) and C/D (Godfrey–Isgur relativized, screened), with
  Feynman–Hellmann parameter derivatives.
- **fit**: Levenberg–Marquardt and Minuit2; the theory floor (Paule–Mandel) solved with the full PDG
  covariance, with its Q-profile interval; held-out refits (leave-one-out, leave-multiplet-out) that
  re-solve the floor, with calibration of the standardised misses; decays in the χ² (widths with
  parametric and theory errors, the held-out decay benchmark, `--fit-widths`); Bayesian inference with
  the floor marginalised and posterior predictive checks (`--mcmc`); joint c c̄ + b b̄ fits with shared
  parameters (`--shared b`) and a likelihood-ratio universality test; the truncation diagnostic of the
  omitted −⟨p⁴⟩/8m³ term; bootstrap, jackknife.
- **data**: PDG SQLite reader: summary-table masses with their common systematics, partial widths or
  summary-table branching fraction × total width.
- **obs**: wave functions (components, weights, norm, radial nodes, expansion coefficients), moments,
  ⟨p⁴⟩, annihilation widths at LO and NLO with the MS-bar α_s and a scale band, S–D mixed leptonic widths.
- **cli**: `qgem`, the end-to-end run over every PDG state with text and JSON output: every number with
  its uncertainty (the B_c from parameter draws), parameter and prediction covariances,
  `--wavefunctions`, `--wf-band`.
- **root**: optional Minuit2 bridge (Migrad, Hesse, Minos), the default minimiser in a ROOT build.
- **analysis**: sweep driver over 80 configurations and the studies; a pre-registered selection rule
  (observables: Pearson χ² of masses at one declared error and of the held-out decays, Holm-corrected;
  physics: nodes, heavy-quark masses, lattice string tension, P-wave spin structure, convergence,
  ⟨v²⟩; the decays in Γ and in ln Γ); information criteria, Wilks tests with Holm and Benjamini–Hochberg, cross-validation, model
  averaging, Kendall τ, limits of the potential description, joint fits, widths, Bayesian posteriors,
  and a validity audit; ROOT figures with uncertainties.
- **docs**: analysis report, method walk-through with the pre-registered rule and references,
  verification map, roadmap, algorithm chart, notes with the v1.0.0 addendum.
- **build**: Eigen 5.0 required.

### Fixed

- The floor equation used the diagonal of the PDG covariance while the χ² used all of it (b b̄ χ² = 8.996
  for ν = 9); it now uses the full covariance.
- Held-out predictions were scored with the full-data floor, so a state set its own error bar; the floor
  is now re-solved without it and the parametric error is included.
- Two Godfrey–Isgur parametrisations that give the same fit were both counted in the physical set.
- Smeared contacts narrower than the basis (σ·r₁ > 1) and collapsing spectra are refused by the χ².
- The registered decay test could not reject an overprediction: its theory error is proportional to the
  prediction, so its pull is bounded below by −1/f (about −2). Found after the final sweeps; the selection
  now also requires the same test on ln Γ (symmetric), so the amendment only adds rejections. With it no
  model passes; the report gives the linear, log and combined verdicts.

[1.0.0]: https://github.com/LeventKaanOguz/Quarkonium-System/releases/tag/v1.0.0
