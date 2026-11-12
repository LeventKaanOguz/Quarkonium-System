# Roadmap

Source: `docs/notes/quarkonium_gem_notes.pdf`, §11.13 (milestones) and §12.2 (reference numbers).
Rule: a step is done only when its test passes. Every reference number comes from the notes (computed independently of
this code) or from an independent method; see [`docs/VALIDATION.md`](docs/VALIDATION.md) for the full map.

| | Step | Done when | |
|---|---|---|---|
| M0 | Skeleton: CMake, Eigen, GoogleTest, clang-format CI | tests pass, also under `asan` | [x] |
| M1 | Closed-form S, T, ⟨rᵏ⟩, Gaussian, erf/r, e^{−μr}; canonical-orthogonalization solver; Models A/B | **Test 1 (18 numbers)**, Airy zeros, own finite-difference solver | [x] |
| M2 | Spin terms: contact, L·S, tensor, coupled S–D blocks, m₁≠m₂ mixing | **Test 2**, centre-of-gravity theorem, Clebsch–Gordan construction | [x] |
| M3 | Observables: R(0), R′(0), R″(0), moments, ⟨p²⟩, ⟨p⁴⟩, ee/γγ/E1/M1 widths, 6j | **Test 5** | [x] |
| M4 | Feynman–Hellmann derivatives (incl. the eigenvector-derivative term) | agree with finite differences for every parameter | [x] |
| M5 | PDG reader, datasets with common systematics, bounded LM, theory floor, covariance | Toy 2; **Test 6** | [x] |
| M6 | Leverage/LOO, AIC/BIC/Wilks, profile, bootstrap, MCMC + R̂/τ_int | exact on the toy; §10.10 (`qgem_validation`) | [x] |
| M7 | Basis optimization (Ky Fan), nested doubling, σ_num | HUM bound; Aitken; fitted masses move < 1e-3 MeV | [x] |
| M8 | Relativistic kinetic energy E₁+E₂ | **Test 3** | [x] |
| M9 | Non-local factors f(p) W(r) f(p): closure and two-step | **Test 4**, both methods | [x] |
| M10 | Full Models C and D (GI, smeared G̃, S̃, S̃_scr, Φ-sandwiches, coupled channels) | structural checks + published spectra (recalled — verify parameters) | [x] |
| M11 | Complex-range Gaussians | table of §7.10.5 | [x] (A/B; not the relativistic matrices) |
| M12 | Priors, error budget, GI fits via the shared fit machinery | budget and prior tests; `qgem fit --model C` | [x] (fit of {α_k, γ_k} of §10.9 not implemented) |
| M13 | Running α_s for every model: frozen and APT shapes as NNLS Gaussian sums (closed-form matrix elements) | shape accuracy ≤ 3.3 % for 0.05 < Q < 50 GeV; Richardson refused | [x] |
| M14 | Physical domain of the χ²: collapsed or mis-ordered levels are a wall; multistart; exact LOO robust to failed refits | collapse test; every sweep fit re-solved without tracking | [x] |
| M15 | CERN ROOT bridge: Minuit2 Migrad/Hesse/Minos as a second minimiser (optional CMake target) | straight line exact; §10.10 fit equal to LM (Δχ² < 1e-6) | [x] |
| M16 | Model sweep and analysis: 80 configurations × {fixed floors, Minuit2, bootstrap/MCMC/jackknife, cross-validation, basis}; PyROOT TTrees and figures | `analysis/`, [`docs/report/REPORT.md`](docs/report/REPORT.md) | [x] |
| M17 | Seeded fits: `--start` points from the other minimiser and from the nested models; final fit = best minimum any method found | seed test; every final fit ≤ min(LM, Minuit2); every Wilks Δχ² ≥ 0 | [x] |
| M18 | Basis by principle: (r₁, r_max) by the variational level sum, n by nested convergence (10, 19, 37, …), η and quadrature checked; basis-resolution wall | dataset Ky Fan and nested-size tests; n = 37 meets 0.1 MeV in 159/160 fits | [x] |
| M19 | Model comparison: AIC/AICc/BIC, Laplace evidence, Wilks with Holm/BH, LOO/LMO/K-fold, learning curves, residual-structure tests, pseudo-BMA and stacking, model averaging, Kendall τ | `analysis/model_selection.py` and its test | [x] |
| M20 | Modified GI: smearing slope s fitted instead of σ₀ (σ₀ and s are degenerate within a system) | same minima, condition 10¹³ → 10⁵ | [x] |
| M21 | Performance: GI momentum matrices memoised across the finite-difference evaluations | 2.75× faster, results identical to 10⁻¹¹ | [x] |
| M22 | State-dependent floor: shapes ⟨v²⟩, ⟨p⁴⟩/8m³, r_rms from a reference model, selected by LOO + LMO predictive density | study; paired differences mostly within 2 SE | [x] |
| M23 | Semi-relativistic kinetic energy √(p²+m²) − m for A/B (running couplings) | exact mass derivative; first-order limit ⟨p⁴⟩/4m³ | [x] |
| M24 | Limits of the potential description: out-of-sample residuals against r_rms, ⟨v²⟩ and the threshold | per-state diagnostics in every run | [x] |
| M25 | Goodness of fit with an estimated floor: full-covariance floor, χ² test inverted into floor intervals, held-out refits with re-solved floors and calibration, flavour scaling of the floor | refits at the interval ends within 1 % | [x] |
| M26 | Truncation test: the a-priori floor from the unabsorbed −⟨p⁴⟩/8m³ term; χ²/ν at that floor; residual-pattern coefficient λ | closure: a parameter direction is absorbed exactly | [x] |
| M27 | Widths: NLO QCD factors with the MS-bar α_s and a scale band; S–D mixed leptonic widths; PDG BR × Γ reference values; ratios | LO/NLO/band printed per state | [x] |
| M28 | Joint c c̄ + b b̄ fit with shared parameters (universal b): per-system floors from leverage shares, likelihood-ratio universality test, B_c from joint parameters | separate limit reproduced; Jacobian against finite differences | [x] |
| M29 | Wave functions logged (components, weights, norm, nodes, coefficients) and a validity audit: spin structure, m_b − m_c against heavy-quark symmetry, √b against the lattice, the ψ(3770) mixing angle, wave functions rebuilt from the JSON | nodes = radial quantum number; virial theorem | [x] |
| M30 | Decays in the χ²: widths with parametric and theory errors, the held-out decay benchmark, the combined masses + widths fit | closure of the combined fit | [x] |
| M31 | Bayesian inference: the theory floor marginalised, posterior predictive checks of masses and decays | posterior σ_th inside the Q-profile interval | [x] |
| M32 | Standard goodness of fit (Pearson χ² at a pooled Paule–Mandel error, Holm) and a pre-registered selection rule | rule fixed before the final sweeps | [x] |
| M33 | Decay test in ln Γ as well as in Γ (the linear test cannot reject an overprediction); both required | amendment after the sweeps, adds rejections only; three verdicts reported | [x] |

Planned for v2.0.0: a momentum-space solver (an independent method for the relativized models C/D); the
log-form decay test and its posterior predictive check in the `qgem` output (v1.0.0 computes the log form
in `analysis/selection.py` from the logged widths, with the parametric error diagonal).

Open items, deliberately left: coupled-channel (meson-loop) effects for the states near and above threshold (a separate model, not a fix); the §10.9 optimization of the running coupling (the machinery — priors, FH/finite-difference
Jacobians, profile/bootstrap/MCMC — is in place; the parametrization with label-switching constraints is not written);
a Model C/D variant driven by an arbitrary running α_s form; complex-range matrices for the relativistic kinetic energy; fitting the frozen/APT shape (M_g, Λ) — the report finds
the fitted α(0) normalisation far from the shape's own value; m/E relativistic factors on the A/B spin terms for an O(v²)-consistent charmonium.

Added after the milestones: exact analytic forms for every model matrix element; a running coupling in all four models
(NR: Gaussian-sum shape with a free saturation value; GI: fixed, free or constant); the additive constant c as a switch;
an end-to-end default run over every PDG state (`qgem`); the M13–M21 rows above.
