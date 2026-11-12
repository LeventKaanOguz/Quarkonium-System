# What is verified, and against what

Each milestone of [`ROADMAP.md`](../ROADMAP.md) closes on a test. The reference is one of: a closed-form result, an
independent numerical method that shares no code with the solver (own finite-difference solver, Clebsch–Gordan
construction, radial convolution by quadrature, Airy zeros), or the reference numbers of the lecture notes
(`docs/notes/`, §12.2), which were computed independently of this code. The sixth column says which.

| test (file) | checks | reference |
|---|---|---|
| `Special.*`, `Quadrature.*` (`test_special`) | `erfcx`, master integral, screened moments Q_k on both recursion branches, Gauss rules, spherical Bessel | long-double evaluation, own Newton-Legendre integrator |
| `Radial.*` (`test_radial`) | S, p², ⟨rᵏ⟩, Gaussian, erf/r, e^{−μr}, mixed L | independent composite quadrature; p² by integration by parts |
| `Solver.*`, `Gep.*`, `Test1.*` (`test_solver`) | canonical orthogonalization, oscillator, hydrogen, ∂M/∂c = 1, nested monotonicity, **Test 1 (18 numbers)** | analytic; notes Test 1 |
| `Independent.*` | linear potential, Cornell and screened for l = 0–2 | **Airy zeros**; **own finite-difference solver** (Sturm bisection + Richardson) |
| `Spin.*` (`test_spin`) | ⟨S₁·S₂⟩, ⟨L·S⟩, ⟨S₁₂⟩ diagonal and off-diagonal, centre-of-gravity sums | built explicitly from Clebsch–Gordan coefficients and spherical harmonics, sphere quadrature |
| `Test2.*`, `NRModel.*`, `CoulombSmeared.*` (`test_model`) | ⟨δ_σ⟩, hyperfine, virial, a, t and 3P_J shifts; FH gradients for every parameter in every spin mode (incl. unequal masses, coupled channels, screened); smeared-Coulomb series | notes Test 2; central finite differences |
| `Test5.*`, `SixJ.*`, `Observables.*` (`test_observables`) | R(0), f_V, Γ_ee, Γ_γγ, E1 (all J), M1, ⟨v²⟩, ⟨p⁴⟩; 6j symbols; R′(0), R″(0) | notes Test 5; 6j orthogonality; finite differences of the radial function |
| `LevenbergMarquardt.*`, `TheoryFloor.*` | LM on the straight line, floor equation | notes Toy 2 (exact) |
| `Test6.*` (`test_fit`) | the whole fit chain: σ_th, θ̂, errors, correlation matrix, Fisher eigenvalues, residuals; Jacobian vs finite differences; nested-doubling error | notes §10.10 |
| `Pdg.*`, `Dataset.*` | PDG reader (summary-table values, the second error from the display text), common-systematic covariance | the PDG file itself; the notes' three gotchas |
| `Stats.*`, `Rng.*` (`test_stats`) | leverages and LOO, profile, bootstrap, MCMC with R̂ and τ_int, χ² survival function, Aitken, error budget | exact on the straight line; known distributions; AR(1) |
| `qgem_validation` (label `slow`) | leverages, LOO with exact refits, profile of m_c, the m_c–c degeneracy (ρ = −0.987), AIC/BIC, bootstrap, posterior, priors | notes §10.10 |
| `NelderMead.*`, `BasisOpt.*`, `Nested*` (`test_hyper`) | basis optimization stays an upper bound; σ_num monotone; wave-function observables converge more slowly than energies | Rosenbrock; a large reference basis |
| `Momentum.*`, `Test3.*`, `Test4.*` (`test_relativistic`) | ⟨f(p)⟩ by Gauss–Laguerre, ⟨p⁴⟩ in closed form (eq. 10.21), the relativistic kinetic energy, f W f by closure **and** by the two-step method | notes Tests 3 and 4 (a sine DVR) |
| `AlphaS.*` (`test_alpha_s`) | the running forms and their limits, the static potential (Richardson slope → b = 8πΛ²/(33−2n_f)), the smeared-Coulomb transform against the analytic erf forms of the Gaussian sum | analytic; notes Ex. 5.3 |
| `Analytic.*` (`test_analytic`) | every smeared matrix element used by the models is **exact**: erf/r, g'/r, tensor, r~, r~'/r, (e^{−μr})~ and its derivative (transfer lemmas with exact Taylor jets), J_s recursion | the pointwise functions integrated by quadrature (`local`, kept only as the reference) |
| `Run.*` (`test_run`) | the default run prints every PDG state, reproduces σ_th = 8.344 MeV; model / running / c choices | notes §10.10; error handling |
| `Smeared.*` (`test_smeared`) | r̃, (e^{−μr})~ and their derivatives in the stable erfcx form | the radial convolution formula (6.9) by quadrature |
| `GI.*` (`test_gi`) | Hermiticity, closure vs two-step, D→C limit, FH, B_c, and the published spectra | **recalled from memory (notes §12.3)**: smoke test only |
| `AlphaS.GaussianSumShapes*` (`test_alpha_s`) | frozen and APT shapes as non-negative Gaussian sums, Σw = 1, ≤ 5 % for 0.05 < Q < 50 GeV; a confining form is refused | the analytic couplings |
| `Chi2Problem.CollapsedSpectrum*` (`test_fit`) | a point-like contact term (σ = 1000 GeV) collapses the 1S0 channel and the χ² refuses the point | physical requirement (H bounded below) |
| `Run.SeedFromAnotherMinimiser*` (`test_run`) | a `--start` point (the Minuit2 minimum of A-APT-smeared b b̄, 2.27 MeV) is fitted and kept over the default start, where LM stalls at σ → 0 (107 MeV) | Minuit2 result |
| `Jackknife.DeleteOneFormulaAndFailedRows` (`test_fit`) | √((n−1)/n Σ(x−x̄)²) on known values; failed refits (NaN rows) skipped | hand computation |
| `Run.FreeFixAndLeaveMultipletOut` (`test_run`) | `--fix`/`--free` move parameters in and out of the fit; `--lmo` leaves out each (n, L) multiplet once; unknown names refused | definition |
| `BasisOpt.DatasetBasisLowersTheLevelSumAndStaysAnUpperBound` (`test_hyper`) | the dataset-level Ky Fan sum is lowered by optimising (r₁, r_max) and stays above the n = 80 reference; out-of-range bases refused | HUM / Ky Fan theorem |
| `BasisSize.NestedSequenceConvergesMonotonicallyAndChoosesTheSmallestSize` (`test_hyper`) | on n = 10, 19, 37, 73 every level decreases (nested bases), the changes shrink, and the chosen n is the smallest below the target | HUM theorem |
| `TheoryFloor.StateDependentWeightsSolveTheWeightedEquation` (`test_fit`) | the floor with per-state weights solves Σ Δ²/(s² + w²σ²) = ν; unit weights reproduce the uniform floor | definition |
| `TheoryFloor.UsesTheFullExperimentalCovariance` (`test_fit`) | the floor solves Δᵀ(V_exp + σ²I)⁻¹Δ = ν with a correlated V_exp, and differs from the diagonal approximation | definition |
| `TheoryFloor.ConfidenceIntervalInvertsTheChiSquareTest` (`test_fit`) | χ²_ν quantiles invert the CDF to 10⁻¹⁰; the 68 % ⊂ 95 % intervals bracket the estimate; closed form for V = s²I | closed form |
| `Test6.CharmoniumWorkedExampleOfTheNotes` (closure part) | a shift along ∂μ/∂b (≤ 1 MeV) is absorbed by b alone: R < 0.01 MeV, σ_a ≈ 0, b moves by the shift | linear algebra |
| `Joint.SharedParametersJacobianAndSeparateLimit` (`test_fit`) | nothing shared reproduces the separate fit and floor (ν_s = N − p); the shared column of the Jacobian equals finite differences; identical systems agree on b (Δχ² = 0); a fixed parameter cannot be shared | definition |
| `Observables.WaveFunctionNormNodesAndOriginValue` (`test_observables`) | the first four S and three P levels of a Cornell model: unit norm, n radial nodes for the n-th level, R(0) from the coefficients | Sturm–Liouville node theorem |
| `Observables.VirialTheoremOfTheCornellEigenstates` (`test_observables`) | 2⟨T⟩ = (4α/3)⟨1/r⟩ + b⟨r⟩ for S, P, D, n ≤ 2: 10⁻⁶ for ground states, 10⁻⁴ up to the third D level (first order in the wave-function error) | virial theorem |
| `Widths.ClosureOfTheCombinedFitAndTheoryCovariance` (`test_fit`) | six widths of the charmonium fit; theory error = scale half-band ⊕ ⟨v²⟩Γ; the scale part correlated only within a process; the model's own widths as data are fitted with χ²_w = 0 and unchanged parameters | closure |
| `Bayes.FloorPosteriorAgreesWithTheQProfileInterval` (`test_fit`) | ln det V equals the determinant; the posterior median of σ_th lies inside the Q-profile 95 % interval; posterior predictive p-values in (0, 1) | Q-profile interval |
| `Widths.QcdCouplingFactorsAndMixedLeptonicWidth` (`test_observables`) | α_s(M_Z) = 0.1180 and two-loop values at m_b, m_c; the one-loop factors; the S–D leptonic width reduces to the pure S and pure D forms | PDG; closed form |
| `Pdg.ReaderGivesTheNotesCharmoniumTable` (derived widths) | χ_c0 → γγ = BR × Γ = 2.43 keV (derived); J/ψ → e⁺e⁻ listed; η_b → γγ "not seen" gives none | PDG summary tables |
| `Test6.CharmoniumWorkedExampleOfTheNotes` (held-out part) | a refit without χ_c0 re-solves its floor (χ²/ν = 1 on the kept data, below the full floor) and gives a parametric error; ν = 0 left gives no floor | definition |
| `SemiRelativistic.MassDerivativeFirstOrderLimitAndHerbstBound` (`test_relativistic`) | √(p²+m²) − m: ∂H/∂m equals finite differences; at small v the shift equals ⟨p⁴⟩/4m³ (5 %); point Coulomb refused | first-order perturbation theory; Herbst bound |
| `Chi2Problem.UnresolvedSmearedContactIsOutsideTheModel` (`test_fit`) | a smeared contact with σ·r₁ > 1 is refused by the χ² (basis-resolution wall) | numerical validity |
| `analysis_model_selection` (CTest, Python) | AIC/AICc/BIC, Laplace evidence (and its refusal outside the prior), predictive densities, t-coverage, the variance-ratio test, stacking, Burnham–Anderson averaging, Holm/BH, Kendall τ against hand values | hand computation |
| `Minuit2.*` (`tests/root/`, ROOT build only) | Migrad, Hesse and Minos on the straight line; the §10.10 fit equals Levenberg–Marquardt | notes Toy 2 (exact); LM |
| `ComplexRange.*` | complex-range matrix elements, nested enlargement, the table of notes §7.10.5 | quadrature of the real functions; own finite-difference solver |

## Things that were found, not assumed

* The notes' hydrogen list "−0.22, −0.05, −0.02469" rounds loosely; the exact values for m_r = 1, κ = 2/3 are −2/9, −1/18, −2/81.
* For the m₁≠m₂ channel list, `enumerate_channels` first produced an invalid ³S₀; found by the Bc test.
* The GI tensor term has the sign of G̃′/r − G̃″ already included (G̃ carries the −4/3); a first version flipped it and put
  χ_c0 on χ_c1. The published fine structure caught it.
* Test 4's 3S level depends on the basis range at the 10⁻⁶ level (4.1247389 for a well-placed basis, 4.1247402–4.1247503
  for others), and the notes do not state their n = 30 range.
* Test 6 reproduces to the digits the notes print only with the smeared contact term in **every** L (`contact_in_all_L`);
  with S waves only the floor is 5.5 MeV, not 8.3.
* The Richardson static potential carries an infrared-divergent additive constant (the notes say "up to an infinite
  constant"), so only its slope is a meaningful check.

## Not verified against an independent reference

* The Godfrey–Isgur parameters themselves (the notes recall them "from memory — verify before use"). The model passes
  Hermiticity, closure/two-step and derivative checks, and reproduces the recalled spectra to a few MeV, but the parameter
  values should be checked against the paper before any physics claim.
* Complex-range bases are supported for the non-relativistic Models A/B; the relativistic momentum matrices (Models C/D)
  need complex Gauss–Laguerre scaling and are not implemented (`momentum_diagonal` throws).
* Running α_s inside the models is the Godfrey–Isgur Gaussian sum (exact erf forms). The other forms of `alpha_s.hpp`
  (frozen, APT, Richardson) have no closed-form potential; they are available as numerical tables only
  (`static_potential`, `SmearedCoulombTable`) and are not wired into the models.
* The momentum-space kernels of Models C/D (E₁+E₂, the Φ and (1+p²/E₁E₂)^{1/2} factors) are one-dimensional
  generalized Gauss–Laguerre integrals, spectrally accurate (Test 3 to 3×10⁻⁷) but not closed forms.
* The Bc prediction interpolates the cc̄ and bb̄ parameters in ln(reduced mass); that is an assumption, not a derivation.
