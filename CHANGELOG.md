# Changelog

All notable changes are listed here. The format follows
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/), the commits follow
[Conventional Commits](https://www.conventionalcommits.org/en/v1.0.0/), and the versions follow
[Semantic Versioning](https://semver.org/spec/v2.0.0.html): the public interface is the `qgem` command
line and its JSON output.

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
