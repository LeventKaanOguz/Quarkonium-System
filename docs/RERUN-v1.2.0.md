# Rerun of the analysis with v1.2.0

All fits below were rerun from the v1.2.0 binary and compared with the archived results
(`analysis/compare_runs.py OLD NEW`; the comparison of physical quantities was made on chi2, masses,
theta and widths, because the raw JSON also holds sign-arbitrary wave-function coefficients).

| sweep | files | result |
|---|---|---|
| main | 80 | median |dM| 5e-6 MeV; 78/80 within 0.01 MeV; worst 0.57 MeV (B-apt-con-none, outside the physical set); parameter sigma median rel. diff 8.5e-8 |
| floor0, floor, floor2, joint, joint2, shape, bgrid, nested | 19, 80, 160, 19, 15, 57, 50, 80 | chi2 identical to <=1.3e-10 (relative) in shape/bgrid/nested; <1e-3 everywhere in the others |
| widths, rel, bopt, oneloop, screen, rep | 19, 12, 8, 8, 16, 4 | agree to 1e-5 or better in parameters and masses |
| studies | 7 | text output identical to the archive |

Differences larger than 1e-3 in parameters occur only in unconstrained configurations
(A-const-con-none, B-apt-con-none, B-frozen-con-none, and three floor2 files), at equal chi2.

Not rerun: `lm` and `minuit` (archived minima were used as starting points, so the minimum search was not
repeated from scratch), `cv` (k-fold, about 20 h on 2 cores) and `stats` (bootstrap/MCMC, about 5 h); their numbers are those of v1.1.0.
Regenerating the 12 main-sweep tables of the paper from the rerun reproduces them, except one last digit of a coefficient.

## Complex-range sweep over all variants of A and B (after v1.2.0)

`cx` was extended from four to all 56 configurations of Models A and B and run with the v1.2.0 binary
(`python3 analysis/sweep.py cx -j 2`; the four from the v1.2.0 run were reused, they agree with the main sweep to 1e-5).
All 112 fits (56 configurations, two systems) converged. Level changes against the real-basis fits (`main`; median over all 56: 0.09 keV):
15 physical-set configurations of A and B: at most 9.7 keV, median 0.19 keV, floors change by at most 9e-5 MeV;
the other 41 (degenerate Fisher matrix): one outlier at 4.3 MeV (A-apt-con-first, b b̄, the contact term of the singlet
S waves), the next ones 0.16 and 0.15 MeV (B-apt-con-none, B-frozen-con-none, c free and no spin terms, flat direction).
Times: 1 s to 49 min per configuration on 2 cores. Archive: `results-archive/results-v1.2.0-cx56.tar.gz`.
