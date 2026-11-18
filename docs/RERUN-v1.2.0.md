# Rerun of the analysis with v1.2.0

All fits below were rerun from the v1.2.0 binary and compared with the archived results
(`analysis/compare_runs.py OLD NEW`; the comparison of physical quantities was made on chi2, masses,
theta and widths, because the raw JSON also holds sign-arbitrary wave-function coefficients).

| sweep | files | result |
|---|---|---|
| main | 80 | median |dM| 5e-6 MeV; 78/80 within 0.01 MeV; worst 0.57 MeV (B-apt-con-none, outside the physical set); parameter sigma median rel. diff 8.5e-8 |
| floor0, floor, floor2, joint, joint2, shape, bgrid, nested | 19, 80, 160, 19, 15, 57, 50, 80 | chi2 identical to <=1.3e-10 (relative) in shape/bgrid/nested; <1e-3 everywhere in the others |
| widths, rel, bopt, oneloop, cx, screen, rep | 19, 12, 8, 8, 4, 16, 4 | agree to 1e-5 or better in parameters and masses |
| studies | 7 | text output identical to the archive |

Differences larger than 1e-3 in parameters occur only in unconstrained configurations
(A-const-con-none, B-apt-con-none, B-frozen-con-none, and three floor2 files), at equal chi2.

Not rerun: `lm` and `minuit` (archived minima were used as starting points, so the minimum search was not
repeated from scratch), `cv` (k-fold, about 20 h on 2 cores) and `stats` (bootstrap/MCMC, about 5 h); their numbers are those of v1.1.0.
Regenerating the 12 main-sweep tables of the paper from the rerun reproduces them, except one last digit of a coefficient.
