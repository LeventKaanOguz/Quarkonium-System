# Changelog

All notable changes are listed here. The format follows
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/), the commits follow
[Conventional Commits](https://www.conventionalcommits.org/en/v1.0.0/), and the versions follow
[Semantic Versioning](https://semver.org/spec/v2.0.0.html): the public interface is the `qgem` command
line and its JSON output.

## [Unreleased]

### Added

- **examples**: `medium_dissociation.py`, an application of the stored wave functions (JSON of `results/`): the screening
  mass μ* at which each S and P state of a Model A (constant coupling) fit dissolves in the Debye-screened potential
  V_T = −C_F α_s e^{−μr}/r + (b/μ)(1 − e^{−μr}), at first order in the vacuum wave function (a variational lower bound
  of μ*) and from a new Gaussian diagonalization. `test_medium_dissociation.py` checks the Coulomb levels, the Yukawa
  critical screening (1.1906, converged from below: 1.1778/1.1856/1.1888 for r_max = 120/300/800 GeV⁻¹) and the
  variational bound. Result for the paper's representative: μ* = 0.82 GeV (J/ψ), 0.37 (ψ(2S)), 0.35 (χ_c1), and 1.56
  (Υ(1S)), 0.65 (Υ(2S)), 0.59 (χ_b1(1P)), 0.43 (Υ(3S)); `results/medium`, archived in `results-v1.2.0.tar.gz`.
- **me**: complex-range Gaussians (`--complex n,rho1,rhomax,omega`, pairs of real functions Re/Im of
  r^l exp(−η(1+iω)r²)) now also work for the smeared and screened potentials of Models A and B: the
  analytic matrix elements fall back to numerical quadrature of the same radial integrals when a basis function
  has a complex range (`analytic.cpp`). The relativistic matrices of Models C/D and the one-loop logarithm are
  real-range only. Test `ComplexRange.SmearedAndScreenedModelsAgreeWithTheRealRangeBasis`.
  `qgem_study complex-levels`: Model A levels of the ³S₁ and ³P₂ channels against a variational reference for real
  and mixed bases at equal numbers of real functions. At 48 real functions the 9th and 10th ³S₁ levels are off by
  9.6 and 35 MeV with the real basis and by 0.03 and 0.27 MeV with 24 real + 12 complex. Fits with 25 real + 6
  complex Gaussians reproduce the real-basis fits (floors identical, levels ≤ 0.0095 MeV).
- **apps**: `qgem --mu-s GeV|fit` and `--log FILE`. `--mu-s fit` fits the screening mass of Models B/D (as
  `--fit-mu-s`); `--shared mu_s` shares it between c c̄ and b b̄; `--log` copies everything printed to a file, after
  a `# qgem <args>` line.
- **analysis**: `sweep.py rep oneloop cx screen studies` (v1.2 runs, seeded from the v1.1.0 minima, a `.log` per run)
  and `summarize_v12.py` (writes `results/v1.2-summary.md`). Result of `screen`: Model B prefers μ_s → 0 in both
  systems (profile rises monotonically; shared 0.000 ± 0.049 GeV); Model D gives μ_s(b b̄) = 0.144 ± 0.016 GeV and
  μ_s(c c̄) → 0 (± 0.19), compatible with a shared 0.142 ± 0.016 GeV (Δχ² = 2.17, P = 0.14). Floors of fits with
  μ_s free are not comparable with fixed-μ_s ones (ν is smaller by one).
  The four representatives rerun with this version reproduce v1.1.0 to 5·10⁻⁸ (relative, parameters) and
  1.3·10⁻⁵ MeV (masses); the 80 configurations of the v1.1.0 sweep were not rerun. 113 tests pass.
- **me/models**: closure-free evaluation of the Godfrey–Isgur f(p) W(r) f(p) sandwiches, `qgem --gi-sandwich
  exact` (`GIModelOptions::exact_sandwich`). Every basis function is Hankel-transformed once, χₙ(r) =
  √(2/π) ∫ p² j_l(pr) f(p) R̃ₙ(p) dp, and ⟨fWf⟩ = ∫ r² χₘ W χₙ dr is integrated on a Gauss–Legendre grid, so the
  matrix is exact for the basis (a rigorous variational bound, no closure in the orthonormal basis). It applies to
  the Coulomb, contact, tensor and spin–orbit sandwiches. Slow (≈ 100 s against 0.2 s for the table below).
  Result at the published BGS (m_c = 1.628) and GM (m_b = 4.977) parameters (`qgem_study gi-levels`): the four
  lowest levels of every channel up to L = 3 move by ≤ 0.03 MeV, the median level by 0.01 MeV; the closure is not
  the origin of the difference to the published S waves (Sec. 7 of the paper). Test
  `GI.ExactSandwichesAgreeWithClosureForTheLowLyingLevels`.
- **apps**: `qgem_study`, small reproducible studies that evaluate a model without a fit: `gi-levels` (GI spectrum at
  given parameters, closure or exact), `contact-r0` (|R(0)|² of the 1S and 2S states with the contact term outside
  and inside H, as a function of the smearing σ).
- **models**: the one-loop short-distance static potential for Models A/B with constant coupling,
  `qgem --one-loop [--mu-ref GeV]` (`NRModelOptions::one_loop_nf`, `one_loop_mu`):
  V = −(C_F α_s/r)[1 + (α_s/4π)(a₁ + 2β₀γ_E + 2β₀ ln μr)], a₁ = 31/3 − 10 n_f/9, β₀ = 11 − 2n_f/3, n_f = 3 for
  c c̄ (and the B_c) and 4 for b b̄. The spin–orbit and tensor terms follow from the same potential (first order,
  `--spin first`): V′/r → κ[(1 + c₀ − c₁) + c₁ ln r]/r³, V′/r − V″ → κ[3(1 + c₀) − 4c₁ + 3c₁ ln r]/r³. The new
  matrix elements ⟨rᵏ ln r⟩ are analytic (digamma). Tests: levels against the independent finite-difference
  solver, the W matrix of ³P₂ against quadrature with numerical derivatives, Feynman–Hellmann gradients. The
  one-loop spin–spin (contact) correction is not included. α_s is the coupling at μ_ref (default 1 GeV); the fits
  depend on it (see the paper, Sec. 6), 2–3 GeV is the sensible range.
- **report**: the B_c prediction carries an interpolation systematic, the largest shift of the mass under three
  alternative interpolation rules (linear in the reduced mass, in its logarithm, in its inverse) of the parameters
  between the c c̄ and b b̄ fits; the JSON has `sigma_interp_GeV`, and the pull of the B_c includes it.

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
