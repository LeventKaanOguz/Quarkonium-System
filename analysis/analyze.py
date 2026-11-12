#!/usr/bin/env python3
"""Analysis of the qgem sweeps.

Reads results/<sweep>/*.json (written by analysis/sweep.py) and writes

  results/qgem_sweep.root      TTrees `fits` and `states`
  docs/report/figures/*.{pdf,png}
  docs/report/data/*.csv       the numbers behind every table of the report
  docs/report/tables.md        all tables

Modules: common (paths, loading, labels), tables (rows from the runs, no
ROOT), model_selection (the statistics), figures (CERN ROOT).
Run with the system Python that ships PyROOT:  python3 analysis/analyze.py
"""
import math
import sys
from pathlib import Path

import numpy as np
from scipy import stats

sys.path.insert(0, str(Path(__file__).resolve().parent))
import model_selection as ms  # noqa: E402
from common import (REPORT, describe, load, md_table, plain, system, write_csv)  # noqa: E402
from figures import (fig_alpha_s, fig_averaged, fig_bgrid, fig_error_methods,
    fig_floor_intervals, fig_learning_curves, fig_selection, fig_wavefunctions, fig_limits, fig_loo_vs_floor,
    fig_minimizers,
    fig_params, fig_pulls, fig_ranking, fig_residuals, fig_spectrum, fig_tau,
    fig_weights, fig_widths, root_file)  # noqa: E402
from sweep import SHAPES, duplicate_fit, pooled_floors  # noqa: E402
import selection  # noqa: E402
from tables import (SCALING_LAWS, averaged_predictions, bgrid_rows,
    binned_residuals, calibration_rows, check_basis, criteria_agreement,
    cv_rows, error_methods, error_summary, floor_analysis, floor_study_rows,
    bayes_rows, common_floor_rows, decay_rows, floor_study_summary, gof_rows, ic_md, identifiability_rows, joint_rows,
    limits_rows, mass_width_rows, truncation_rows, width_ratio_rows, width_rows,
    merge_systems, minuit_analysis, nested_rows, nested_tests,
    overfitting_rows, phys_tags, predictive_analysis, rel_rows,
    residual_structure_rows, scaling_rows, summarize, universality_rows,
    wilks_summary)  # noqa: E402


# lattice QCD, n_f = 2: r_c = 2.5 r_0 with r_0 = 0.5 fm (Bali et al., SESAM,
# Phys. Rev. D 71, 114513 (2005), hep-lat/0505012)
R_STRING_BREAK_FM = 1.25


def gof_section(main_runs, floor0_runs, floor2_runs, cv_runs, stats_runs,
                nested_runs, phys):
    """Model selection (pre-registered rule) and the goodness of fit: Pearson
    chi^2 of the masses at the declared error, the decay benchmark, the
    Bayesian posterior predictive checks; then the diagnostics."""
    out = ""
    pf = pooled_floors()
    rows, win, best = selection.select(main_runs, floor2_runs, nested_runs, phys, pf)
    write_csv("selection.csv", [{**r, "fails": "; ".join(r["fails"])} for r in rows])
    fig_selection(rows)
    out += ("## Model selection (pre-registered rule)\n\n"
            "A model wins only if it (1) reproduces the observables, masses and decays: "
            "no rejection at α = 0.05 after Holm's correction over its four Pearson χ² "
            "tests, and (2) is physically sound: non-degenerate, basis converged, radial "
            "nodes = radial quantum number, m_b − m_c in the heavy-quark window, √b in "
            "0.44 ± 0.04 GeV (lattice), 1P centre of gravity − h at the PDG value, "
            "⟨v²⟩ < 0.5; each within 2σ. Among those, the lowest AIC = χ²_masses + "
            "χ²_decays + 2p. The rule and its thresholds were fixed before the final "
            "sweeps (docs/EXPLANATION.md).\n\n" + md_table(
                rows, ["label", "p_mass_ccbar", "p_mass_bbar", "p_decay_ccbar",
                       "p_decay_bbar", "aic", "n_fail", "verdict"],
                ["configuration", "p masses cc̄", "p masses bb̄", "p decays cc̄",
                 "p decays bb̄", "AIC", "failed checks", "verdict"],
                ["", ".3g", ".3g", ".3g", ".3g", ".1f", "d", ""]) + "\n")
    out += ("**Winner: " + (f"{win['label']}** (AIC {win['aic']:.1f})." if win else
            f"none.** No model passes both gates. Best compromise: {best['label']} "
            f"({best['n_fail']} failed check{'s' if best['n_fail'] != 1 else ''}: "
            f"{'; '.join(best['fails'])}).") + "\n\nFailed checks per model:\n\n" +
            "".join(f"- {r['label']}: {'; '.join(r['fails']) or 'none'}\n" for r in rows) +
            "\n")
    crows0 = common_floor_rows(floor2_runs, phys)
    if crows0:
        write_csv("chi2_common_floor.csv", crows0)
        out += ("## Goodness of fit: masses (Pearson χ² at the declared error)\n\n"
                f"Every model refitted with the declared theory error of the masses, c c̄ "
                f"{pf['ccbar']:.2f} MeV and b b̄ {pf['bbar']:.2f} MeV: the Paule–Mandel "
                "random-effects estimator (Paule and Mandel 1982) pooled over the physical "
                "set, so one error budget for all models, favouring none. χ²/ν and its "
                "p-value P(χ²_ν ≥ χ²) are the conventional test. Sorted by the combined "
                "p-value.\n\n" + md_table(
                    crows0, ["label", "ccbar_chi2nu", "ccbar_p", "bbar_chi2nu", "bbar_p",
                             "chi2nu", "p"],
                    ["configuration", "χ²/ν cc̄", "p", "χ²/ν bb̄", "p", "χ²/ν both", "p"],
                    ["", ".2f", ".3g", ".2f", ".3g", ".2f", ".3g"]) +
                f"\n{sum(r['p'] >= 0.05 for r in crows0)} of {len(crows0)} models pass at "
                "5 % (combined).\n\n")
    drows = decay_rows(floor2_runs, stats_runs, phys)
    if drows:
        write_csv("decay_benchmark.csv", drows)
        out += ("## Goodness of fit: decays (held-out benchmark)\n\n"
                "The measured annihilation widths of the fitted states, predicted by the "
                "mass-only fit at the declared error (the fit graded above), which never saw "
                "them: Pearson χ² with V = V_exp + V_theory + "
                "V_par. Theory error: the truncation of the α_s series by scale variation "
                "(μ = M/4 … M, half the band, correlated within a process) in quadrature "
                "with the O(v²) relativistic correction of NRQCD power counting, ⟨v²⟩Γ. "
                "ppp: the Bayesian posterior predictive p-value of the same widths (stats "
                "sweep).\n\n" + md_table(
                    drows, ["label", "system", "n", "chi2n", "p", "ppp", "worst",
                            "worst_pull"],
                    ["configuration", "system", "widths", "χ²/N", "p", "ppp",
                     "largest pull", "pull"],
                    ["", "", "d", ".2f", ".3g", ".3g", "", ".2f"]) + "\n")
    brows = bayes_rows(stats_runs, phys)
    if brows:
        write_csv("bayes.csv", brows)
        out += ("## Bayesian inference: the theory floor marginalised\n\n"
                "Posterior of the parameters and of σ_th (flat priors in u, i.e. in ln θ "
                "for positive parameters; flat in ln σ_th, scale invariant), Gaussian "
                "likelihood with its normalisation; 4 Metropolis chains. Posterior "
                "predictive p-values (Gelman, Meng and Stern 1996) for T = χ² and T = "
                "max|z|: small means the data are more discrepant than the model's own "
                "replicas. With σ_th marginalised the χ² check is not sensitive to the "
                "overall scale; max|z| flags single states.\n\n" + md_table(
                    brows, ["label", "system", "sth_med", "sth_lo", "sth_hi", "ppp_chi2",
                            "ppp_max", "rhat_max", "neff_min"],
                    ["configuration", "system", "σ_th median (MeV)", "16 %", "84 %",
                     "ppp χ²", "ppp max|z|", "max R̂", "min n_eff"],
                    ["", "", ".2f", ".2f", ".2f", ".3f", ".3f", ".3f", ".0f"]) + "\n")
    trows = truncation_rows(main_runs, phys)
    if trows:
        write_csv("truncation.csv", trows)
        out += ("## Diagnostic: the truncation estimate (A/B, non-relativistic kinetic energy)\n\n"
                "Not a grade (the grades above use standard methods); a physics diagnostic. "
                "σ_a is the part of the first omitted "
                "term, −⟨p⁴⟩(1/8m₁³ + 1/8m₂³), that the model's parameters cannot absorb "
                "(the model's own values shifted by the term and refitted; σ_a solves the "
                "floor equation for the unabsorbed residuals R). At σ_a, χ²/ν is the "
                "conventional test. λ: (PDG − model) = λR, 1 if the omitted term is the "
                "misfit, 0 if unrelated. Darwin, orbit-orbit and 1/m potential terms are "
                "of the same order and also omitted, so σ_a estimates one term of the "
                "truncation, not all of it. Per model, σ_a is a self-consistency check, not "
                "a ranking: a less flexible model absorbs less of the term and passes more "
                "easily.\n\n" + md_table(
                    trows, ["label", "system", "sigma_a", "floor", "ratio", "nu", "chi2_nu",
                            "p", "lam", "lam_err"],
                    ["configuration", "system", "σ_a (MeV)", "fitted σ_th (MeV)",
                     "σ_a/σ_th", "ν", "χ²/ν at σ_a", "p", "λ", "±"],
                    ["", "", ".2f", ".2f", ".2f", "d", ".2f", ".3f", ".2f", ".2f"]) + "\n")
        for sname in ("ccbar", "bbar"):
            rs = [r for r in trows if r["system"] == sname]
            if rs:
                out += (f"{sname}: median σ_a/σ_th = {np.median([r['ratio'] for r in rs]):.2f}, "
                        f"median χ²/ν = {np.median([r['chi2_nu'] for r in rs]):.2f}, "
                        f"p < 0.05 in {sum(r['p'] < 0.05 for r in rs)} of {len(rs)}; "
                        f"median λ = {np.median([r['lam'] for r in rs]):.2f}.\n\n")
    grows = gof_rows(main_runs, floor0_runs, phys)
    write_csv("goodness_of_fit.csv", grows)
    fig_floor_intervals(grows)
    worst_p = max(stats.chi2.logsf(r.get("chi2_pdg_refit", r["chi2_pdg"]), r["nu"])
                  for r in grows) / math.log(10)
    out += ("## Diagnostic: each model's own floor (physical set)\n\n"
           "σ_th: the Paule–Mandel estimator per model (χ² = N − p); its interval is the "
           "Q-profile interval (Viechtbauer 2007), the χ² test inverted. χ²_PDG: "
           "Cochran's Q, PDG errors only, at the floor fit and refitted (σ_th = 0). The "
           f"largest p-value over the set is 10^{worst_p:.0f}: no model reaches the "
           "precision of the data. σ_th: the χ² test of 'the model is accurate to "
           "σ_th' inverted (every σ_th it accepts at 68 % / 95 % CL). Shift: largest "
           "parameter change of the PDG-only refit, in units of the floor fit's error. "
           "t: standardised LOO miss, every term of its error from the refit without "
           "the state.\n\n" + md_table(
               grows, ["label", "system", "nu", "chi2_pdg", "chi2_pdg_refit", "floor",
                       "ci68_lo", "ci68_hi", "ci95_lo", "ci95_hi", "shift_sigma",
                       "shift_par", "loo_rms_t", "loo_rms_t_expected"],
               ["configuration", "system", "ν", "χ²_PDG", "χ²_PDG refit", "σ_th",
                "68 % lo", "68 % hi", "95 % lo", "95 % hi", "shift (σ)", "of",
                "RMS t (LOO)", "expected"],
               ["", "", "d", ".3g", ".3g", ".2f", ".2f", ".2f", ".2f", ".2f", ".1f", "",
                ".2f", ".2f"]) +
           "\n(σ_th and its intervals in MeV. The interval is the variance interval "
           "ν s²/σ² ~ χ²_ν: it assumes model errors independent across states, which "
           "the runs test below checks. The LOO RMS t mainly flags outliers: the floor "
           "was fitted to the same states.)\n\n")
    crows = calibration_rows(main_runs, cv_runs, phys)
    write_csv("calibration.csv", crows)
    out += ("### Out-of-sample calibration (pooled over the physical set)\n\n" + md_table(
        crows, ["system", "kind", "n", "rms_t", "rms_t_expected",
                "median_rms_t_per_model", "in1", "in1_expected", "in2", "in2_expected"],
        ["system", "held out", "points", "RMS t", "expected", "median RMS t per model",
         "|t| < 1", "expected", "|t| < 2", "expected"],
        ["", "", "d", ".2f", ".2f", ".2f", ".2f", ".2f", ".2f", ".2f"]) +
        "\n(expected: a t distribution with the refit's degrees of freedom; – : no finite "
        "RMS for ν ≤ 2. Points of different models share states, so the counts are "
        "descriptive. LMO is the external test of the scale: whole multiplets are "
        "predicted by a floor fitted without them.)\n\n")
    srows = scaling_rows(main_runs, phys)
    write_csv("floor_scaling.csv", srows)
    summ = []
    for law in SCALING_LAWS:
        ps = [r[f"{law}_p"] for r in srows]
        summ.append({"law": law, "models": len(ps),
                     "ratio": float(np.median([r[f"{law}_ratio"] for r in srows])),
                     "median_p": float(np.median(ps)),
                     "rejected": sum(p < 0.05 for p in ps)})
    out += ("### Does one accuracy scale describe both flavours?\n\n"
            "Floor of state i = κ·u_i with u = 1 (uniform), ⟨v²⟩_i, or "
            "⟨p⁴⟩_i/8m³ (the first kinetic term the non-relativistic Hamiltonian "
            "omits), from each model's own fit. κ solved per system; under one κ for "
            "both, F = (χ²_c/ν_c)/(χ²_b/ν_b) ~ F(ν_c, ν_b). With ν_c = 4 the test is "
            "weak: at 95 % it accepts κ_c/κ_b within about 0.33–2.2, so it rejects a "
            "wrong flavour scaling by a factor ≳ 3 but cannot tell ⟨v²⟩ from "
            "⟨p⁴⟩/8m³.\n\n" + md_table(
                summ, ["law", "models", "ratio", "median_p", "rejected"],
                ["floor ∝", "models", "median κ_c/κ_b", "median p", "rejected at 5 %"],
                ["", "d", ".2f", ".2f", "d"]) + "\n")
    return out


def main():
    main_runs = load("main")
    basis = check_basis({sw: load(sw) for sw in
                         ("main", "floor", "cv", "stats", "nested", "bopt")}
                        | {"floor2": load("floor2")})
    print("basis (n, r1, rmax):", basis)
    rows = summarize(main_runs)
    dups = [r["tag"] for r in rows if duplicate_fit(r["tag"], [x["tag"] for x in rows])]
    print("same fit as their -s form, left out of the physical set:", dups or "none")
    write_csv("ranking.csv", rows)
    root_file(main_runs, rows)
    tables = "# Tables generated by analysis/analyze.py\n\n"

    tables += "## Ranking by theory floor (main sweep)\n\n" + md_table(
        rows,
        ["label", "ccbar_floor", "bbar_floor", "ccbar_loo_rms", "bbar_loo_rms",
         "ccbar_pred_rms", "bbar_pred_rms", "bc_diff", "degenerate"],
        ["configuration", "σ_th cc̄", "σ_th bb̄", "LOO cc̄", "LOO bb̄",
         "pred RMS cc̄", "pred RMS bb̄", "B_c − PDG", "degenerate"],
        ["", ".1f", ".1f", ".1f", ".1f", ".1f", ".1f", ".0f",
         lambda v: "yes" if v else ""]) + "\n(all in MeV)\n\n"

    tables += gof_section(main_runs, load("floor0"), merge_systems(load("floor2")),
                          load("cv"), load("stats"), load("nested"), phys_tags(rows))
    fig_ranking(rows)
    fig_loo_vs_floor(rows)
    fig_alpha_s()
    # physical comparison set: c = 0 (c is nonphysical, degenerate with m_Q),
    # spin terms on, Fisher matrix not degenerate, one per distinct fit
    good = phys_tags(rows)
    top = good[:4]
    if top:
        fig_residuals(main_runs, top, "ccbar", "residuals_ccbar")
        fig_residuals(main_runs, top, "bbar", "residuals_bbar")
        fig_spectrum(main_runs[top[0]], "ccbar", "spectrum_ccbar",
                     f"charmonium, {describe(top[0])}")
        fig_spectrum(main_runs[top[0]], "bbar", "spectrum_bbar",
                     f"bottomonium, {describe(top[0])}")
        fig_pulls(main_runs[top[0]], "pulls_best")
        fig_params(main_runs, good[:10])

    floor_runs = load("floor")
    if floor_runs:
        frows = floor_analysis(floor_runs)
        write_csv("information_criteria.csv", frows)
        tables += "## Common 10 MeV floor: information criteria (both systems summed)\n\n" + ic_md(frows[:25])
        tests = nested_tests(frows)
        write_csv("nested_tests.csv", tests)
        tables += "## Nested-model tests (Wilks; Chernoff ½ on a boundary)\n\n" + md_table(
            tests, ["null", "alt", "system", "dchi2", "dp", "p_value", "p_holm", "p_bh",
                    "verdict"],
            ["null", "alternative", "system", "Δχ²", "Δp", "p-value", "Holm", "BH",
             "verdict"],
            ["", "", "", ".2f", "d", ".3g", ".3g", ".3g", ""]) + "\n"
        tables += wilks_summary(tests, "10 MeV")

    fl = pooled_floors()  # the declared mass error the floor2 sweep used
    FL = f"{fl['ccbar']:.1f}/{fl['bbar']:.1f}"
    floor2_runs = merge_systems(load("floor2"))
    f2rows = []
    if floor2_runs:
        f2rows = floor_analysis(floor2_runs)
        write_csv("information_criteria_floor2.csv", f2rows)
        tables += (f"## At the declared mass error (c c̄ {fl['ccbar']:.2f} MeV, b b̄ "
                   f"{fl['bbar']:.2f} MeV, Paule–Mandel pooled over the physical set): "
                   "information criteria"
                   "\n\n" + ic_md(f2rows[:25]))
        tests2 = nested_tests(f2rows)
        write_csv("nested_tests_floor2.csv", tests2)
        tables += wilks_summary(tests2, f"{FL} MeV")

    mn_runs = load("minuit")
    lm_runs = load("lm")
    if mn_runs and lm_runs:
        mrows = minuit_analysis(main_runs, lm_runs, mn_runs)
        write_csv("minimizers.csv", mrows)
        fig_minimizers(mrows)
        tables += "## Minimisers: σ_th (MeV) reached by each method\n\n" + md_table(
            mrows, ["label", "system", "lm_single_floor", "lm_multistart_floor",
                    "minuit_floor", "final_floor", "minuit_valid",
                    "max_minos_asym"],
            ["configuration", "system", "LM, 1 start", "LM, 9 starts",
             "Minuit2, 1 start", "final (seeded)", "Migrad valid",
             "max Minos asymmetry"],
            ["", "", ".1f", ".1f", ".1f", ".1f", "d", ".2f"]) + "\n"
        n = len(mrows)
        stuck_lm = sum(r["lm_ms_excess"] > 0.05 for r in mrows)
        stuck_mn = sum(r["minuit_excess"] > 0.05 for r in mrows)
        tables += (f"Of {n} system fits, LM with 9 starts ends > 5 % above the final "
                   f"floor in {stuck_lm}; Minuit2 from one start in {stuck_mn}.\n\n")

    stats_runs = load("stats")
    if stats_runs:
        erows = error_methods(stats_runs, mn_runs)
        write_csv("error_methods.csv", erows)
        fig_error_methods(erows)
        tables += "## Four error estimates\n\n" + md_table(
            erows, ["label", "system", "param", "value", "linear", "minos_lo",
                    "minos_hi", "bootstrap", "mcmc", "jackknife", "mcmc_rhat"],
            ["configuration", "system", "param", "value", "linear", "Minos −",
             "Minos +", "bootstrap", "MCMC", "jackknife", "R̂"],
            ["", "", "", ".4f", ".4f", ".4f", ".4f", ".4f", ".4f", ".4f", ".3f"]) + "\n"
        tables += error_summary(erows)

    # ---- residual structure (underfitting), K-fold, learning curves, basis plane
    srows = residual_structure_rows(main_runs, [r["tag"] for r in rows if r["spin"] != "none"])
    write_csv("residual_structure.csv", srows)
    tables += ("## Residual structure (underfitting checks)\n\n" + md_table(
        [r for r in srows if r["tag"] in phys_tags(rows)],
        ["label", "ccbar_runs", "ccbar_runs_p", "bbar_runs", "bbar_runs_p",
         "ccbar_runs_all_p", "bbar_runs_all_p", "pull_rms", "pull_normal_p"],
        ["configuration", "runs cc̄", "p cc̄", "runs bb̄", "p bb̄", "p cc̄ (all states)",
         "p bb̄ (all states)", "RMS LOO pull", "Shapiro–Wilk p"],
        ["", "d", ".2f", "d", ".2f", ".2f", ".2f", ".2f", ".2f"]) +
        "\n(runs: Wald–Wolfowitz test on the signs of model − PDG ordered by mass, exact; "
        "few runs = residuals drifting together, the signature of a missing term)\n\n")
    cvr = load("cv")
    crows = cv_rows(cvr, phys_tags(rows))
    if crows:
        write_csv("kfold.csv", crows)
        tables += ("## Repeated 4-fold cross-validation (physical set; MeV)\n\n" + md_table(
            crows, ["label", "ccbar_kf_train", "ccbar_kf_test", "ccbar_kf_median",
                    "bbar_kf_train", "bbar_kf_test", "bbar_kf_median"],
            ["configuration", "train cc̄", "test cc̄", "test median cc̄", "train bb̄",
             "test bb̄", "test median bb̄"], ["", *[".1f"] * 6]) + "\n")
        fig_learning_curves(cvr, phys_tags(rows)[:4])
    brows = bgrid_rows(load("bgrid"))
    if brows:
        write_csv("basis_grid.csv", brows)
        fig_bgrid(brows)
        tables += ("## Basis plane (n = 37): level sum above its minimum and the fit\n\n" +
                   md_table(sorted(brows, key=lambda r: (r["label"], r["system"], r["r1"], r["rmax"])),
                            ["label", "system", "r1", "rmax", "trace_excess_MeV", "floor",
                             "num_all"],
                            ["configuration", "system", "r₁", "r_max", "level sum − min (MeV)",
                             "σ_th (MeV)", "max σ_num (MeV)"],
                            ["", "", "g", "g", ".3g", ".3f", ".3g"]) + "\n")
    bo = load("bopt")
    if bo:
        borows = []
        for t, r in bo.items():
            for sname in ("ccbar", "bbar"):
                an = system(r, sname)["analysis"]
                borows.append({"label": plain(t), "system": sname, "r1": an["basis_r1"],
                               "rmax": an["basis_rmax"],
                               "floor": 1e3 * system(r, sname)["sigma_th_GeV"]})
        write_csv("basis_optimum.csv", borows)
        tables += ("## Variational basis optimum per model and system (n = 37)\n\n" + md_table(
            borows, ["label", "system", "r1", "rmax", "floor"],
            ["configuration", "system", "r₁ (GeV⁻¹)", "r_max (GeV⁻¹)", "σ_th (MeV)"],
            ["", "", ".3f", ".1f", ".3f"]) + "\n")

    # ---- predictive weights, model averaging, agreement of the criteria
    phys = phys_tags(rows)
    irows = identifiability_rows(main_runs, rows)
    write_csv("identifiability.csv", irows)
    shown = [r for r in irows if r["combination"] != "covariance not finite"]
    nonfinite = sorted({r["label"] for r in irows} - {r["label"] for r in shown})
    tables += ("## Identifiability: what a degenerate fit cannot fix\n\n" + md_table(
        shown, ["label", "system", "combination", "rel_sigma"],
        ["configuration", "system", "unconstrained direction (relative units)",
         "its relative σ"], ["", "", "", ".2g"]) +
        (f"\nNo finite covariance at all (a parameter with no effect, e.g. σ "
         f"without spin terms): {', '.join(nonfinite)}.\n" if nonfinite else "") + "\n")
    urows = universality_rows(main_runs, phys)
    write_csv("universality.csv", urows)
    ucols = [k for k in ("alpha_s_z", "b_z", "sigma_z", "s_z", "sigma0_z", "alpha1_z")
             if any(k in r for r in urows)]
    tables += ("## Flavour independence: c c̄ vs b b̄ parameters (z-scores)\n\n" + md_table(
        urows, ["label", *ucols], ["configuration", *[c[:-2] for c in ucols]],
        ["", *[".1f"] * len(ucols)]) + "\n")
    cv_runs = load("cv")
    orows = overfitting_rows(main_runs, cv_runs, [r["tag"] for r in rows])
    write_csv("overfitting.csv", orows)
    po = [r for r in orows if r["tag"] in phys]
    tables += ("## Overfitting checks (physical set; MeV)\n\n" + md_table(
        po, ["label", "ccbar_fit_rms", "ccbar_loo_rms", "ccbar_lmo_rms",
             "ccbar_inflation_ratio", "bbar_fit_rms", "bbar_loo_rms", "bbar_lmo_rms",
             "bbar_inflation_ratio"],
        ["configuration", "fit RMS cc̄", "LOO cc̄", "LMO cc̄", "LOO infl. obs/exp cc̄",
         "fit RMS bb̄", "LOO bb̄", "LMO bb̄", "LOO infl. obs/exp bb̄"],
        ["", ".1f", ".1f", ".1f", ".2f", ".1f", ".1f", ".1f", ".2f"]) +
        "\n(LMO: leave-multiplet-out, every state of one (n, L) multiplet removed "
        "together; obs/exp: observed LOO/fit RMS ratio over the ratio "
        "√(Σ 1/(1−h_i)/(N−p)) expected for a correctly specified linearised model)\n\n")
    prows = predictive_analysis(main_runs, phys)
    if prows:
        write_csv("predictive.csv", prows)
        tables += ("## Leave-one-out predictive density (physical set)\n\n" + md_table(
            prows, ["label", "elpd_ccbar", "elpd_bbar", "elpd", "pseudo_bma_weight",
                    "stacking_weight"],
            ["configuration", "elpd c c̄", "elpd b b̄", "elpd", "pseudo-BMA weight",
             "stacking weight"],
            ["", ".2f", ".2f", ".2f", ".3f", ".3f"]) +
            "\n(elpd = Σ log N(y − μ⁽⁻ⁱ⁾; 0, σ_th⁽⁻ⁱ⁾² w² + σ_exp² + σ_par⁽⁻ⁱ⁾²) over the "
            "left-out states; floor and parametric error from the refit without the state)\n\n")
    if prows and f2rows and floor_runs:
        def sub_w(frs, key, f):
            fr = {r["tag"]: r for r in frs if r["tag"] in phys}
            w = f([fr[t][key] for t in phys])
            return dict(zip(phys, map(float, w)))
        pr = {r["tag"]: r for r in prows}
        wsets = {
            "AIC (10 MeV)": sub_w(frows, "aic", ms.ic_weights),
            f"AICc ({FL})": sub_w(f2rows, "aicc", ms.ic_weights),
            f"BIC ({FL})": sub_w(f2rows, "bic", ms.ic_weights),
            f"evidence ({FL})": sub_w(f2rows, "lnz", ms.evidence_weights),
            "pseudo-BMA": {t: pr[t]["pseudo_bma_weight"] for t in phys},
            "stacking": {t: pr[t]["stacking_weight"] for t in phys},
        }
        wrows = [{"tag": t, "label": plain(t), **{k: w[t] for k, w in wsets.items()}}
                 for t in phys]
        wrows.sort(key=lambda r: -sum(r[k] for k in wsets))
        write_csv("model_weights.csv", wrows)
        fig_weights(wrows, [(k, k) for k in wsets])
        tables += ("## Model weights over the physical set, by criterion\n\n" +
                   md_table(wrows, ["label", *wsets], ["configuration", *wsets],
                            ["", *[".3f"] * len(wsets)]) + "\n")
        avg_sets = {"AICc": wsets[f"AICc ({FL})"], "pBMA": wsets["pseudo-BMA"],
                    "stack": wsets["stacking"],
                    "equal": {t: 1.0 / len(phys) for t in phys}}
        arows = averaged_predictions(main_runs, avg_sets)
        write_csv("averaged_predictions.csv", arows)
        fig_averaged(arows, list(avg_sets))
        tables += ("## Model-averaged predictions (physical set; MeV)\n\n" + md_table(
            arows, ["state", "pdg", *[f"{k}_diff" for k in avg_sets],
                    *[f"{k}_unc" for k in avg_sets], "AICc_spread"],
            ["state", "PDG", *[f"{k}: avg − PDG" for k in avg_sets],
             *[f"{k}: ±" for k in avg_sets], "model spread (AICc)"],
            ["", ".1f", *[".1f"] * (2 * len(avg_sets) + 1)]) +
            "\n(± = Burnham–Anderson unconditional error: Σ w √(σ² + (x − x̄)²), "
            "σ² = σ_th² + σ_par² + σ_num²)\n\n")
        fr10 = {r["tag"]: r for r in frows}
        fr2 = {r["tag"]: r for r in f2rows}
        rr = {r["tag"]: r for r in rows}
        oo = {r["tag"]: r for r in orows}
        scores = {
            "σ_th (iterated)": {t: rr[t]["score"] for t in phys},
            "AIC 10": {t: fr10[t]["aic"] for t in phys},
            "AICc 10": {t: fr10[t]["aicc"] for t in phys},
            "BIC 10": {t: fr10[t]["bic"] for t in phys},
            "−ln Z 10": {t: -fr10[t]["lnz"] for t in phys},
            f"AIC {FL}": {t: fr2[t]["aic"] for t in phys},
            f"AICc {FL}": {t: fr2[t]["aicc"] for t in phys},
            f"BIC {FL}": {t: fr2[t]["bic"] for t in phys},
            f"−ln Z {FL}": {t: -fr2[t]["lnz"] for t in phys},
            "LOO RMS": {t: math.hypot(rr[t]["ccbar_loo_rms"], rr[t]["bbar_loo_rms"])
                        for t in phys},
            "LMO RMS": {t: math.hypot(oo[t]["ccbar_lmo_rms"], oo[t]["bbar_lmo_rms"])
                        for t in phys},
            "−elpd LOO": {t: -pr[t]["elpd"] for t in phys},
            "above-thr. RMS": {t: math.hypot(rr[t]["ccbar_pred_rms"],
                                             rr[t]["bbar_pred_rms"]) for t in phys},
            "|B_c − PDG|": {t: abs(rr[t]["bc_diff"]) for t in phys},
        }
        names, M, used = criteria_agreement(scores)
        fig_tau(names, M)
        tau_rows = [{"criterion": a, **{b: float(M[i, j]) for j, b in enumerate(names)}}
                    for i, a in enumerate(names)]
        write_csv("criteria_tau.csv", tau_rows)
        winners = []
        for nm, sc in scores.items():
            order = sorted(sc, key=lambda t: sc[t] if math.isfinite(sc[t])
                           else math.inf)
            winners.append({"criterion": nm, "first": plain(order[0]),
                            "second": plain(order[1]), "third": plain(order[2])})
        write_csv("criteria_winners.csv", winners)
        tables += (f"## Criteria compared: Kendall τ over the physical set ({len(used)} models)\n\n"
                   + md_table(tau_rows, ["criterion", *names], ["", *names],
                              ["", *[".2f"] * len(names)]) + "\n" +
                   md_table(winners, ["criterion", "first", "second", "third"],
                            ["criterion", "1st", "2nd", "3rd"], ["", "", "", ""]) + "\n")

    # ---- universal string tension: joint c c-bar + b b-bar fits
    for name, title in (("joint", "one string tension b"),
                        ("joint2", "b and the flavour-independent coupling or smearing")):
        jrows = joint_rows(main_runs, load(name), phys)
        if not jrows:
            continue
        write_csv(f"{name}.csv", jrows)
        tables += (f"## Joint c c̄ + b b̄ fit: {title}\n\n"
                   "Each system keeps its floor, solved from its own share of the degrees "
                   "of freedom (N_s − Σ leverages). Universality: χ²(joint) − χ²(separate) "
                   "at the separate floors ~ χ²_k for k shared parameters. B_c: from the "
                   "separate and the joint parameters.\n\n" + md_table(
                       jrows, ["label", "shared", "b_c", "b_b", "b_joint", "b_joint_err",
                               "dchi2", "p", "floor_c_sep", "floor_c", "floor_b_sep",
                               "floor_b", "bc_sep", "bc_joint"],
                       ["configuration", "shared", "b c c̄", "b b b̄", "b joint", "±",
                        "Δχ²", "p", "σ_th cc̄ sep.", "joint", "σ_th bb̄ sep.", "joint",
                        "B_c − PDG sep.", "joint"],
                       ["", "", ".4f", ".4f", ".4f", ".4f", ".2f", ".3g", ".2f", ".2f",
                        ".2f", ".2f", ".0f", ".0f"]) +
                   f"\n(b in GeV², floors and B_c in MeV. Universality rejected at 5 % in "
                   f"{sum(r['p'] < 0.05 for r in jrows)} of {len(jrows)}; median |B_c − PDG| "
                   f"{np.nanmedian([abs(r['bc_sep']) for r in jrows]):.0f} MeV separate, "
                   f"{np.nanmedian([abs(r['bc_joint']) for r in jrows]):.0f} MeV joint.)\n\n")

    # ---- masses and widths fitted together (study)
    mwr = mass_width_rows(main_runs, load("widths"), phys)
    if mwr:
        write_csv("mass_width_fit.csv", mwr)
        tables += ("## Masses and widths fitted together (study)\n\n"
                   "The measured widths of the fitted states added to the χ² with their "
                   "experimental and theory errors (fixed); the mass floor iterated on the "
                   "mass block (leverage share). The selection uses the mass-only fits, with "
                   "the widths as a held-out benchmark, so that the masses are not traded "
                   "for the widths; this table shows what the widths would pull.\n\n" +
                   md_table(mwr, ["label", "system", "floor_mass_only", "floor_combined",
                                  "chi2w_n", "p_w", "shift_par", "shift"],
                            ["configuration", "system", "σ_th masses only (MeV)",
                             "with widths", "χ²/N widths", "p", "largest shift", "(σ)"],
                            ["", "", ".2f", ".2f", ".2f", ".3g", "", ".2f"]) + "\n")

    # ---- annihilation widths: LO, NLO, PDG; ratios
    wtags = phys[:4]
    wrows = width_rows(main_runs, wtags)
    if wrows:
        write_csv("widths.csv", width_rows(main_runs, phys))
        tables += ("## Annihilation widths (best four of the physical set; keV)\n\n"
                   "LO: Van Royen–Weisskopf (vector states: |R_S(0) + 5R″_D(0)/(2√2 m²)|², "
                   "so the S–D mixed ψ(3770) is included). NLO: × one-loop QCD factor with "
                   "the MS-bar α_s at μ = M/2 (band μ = M/4 … M); the series converges "
                   "poorly in charmonium. PDG values marked * are summary-table branching "
                   "fraction × total width.\n\n" + md_table(
                       wrows, ["label", "state", "process", "lo", "nlo", "nlo_err_par",
                               "nlo_err_th", "nlo_lo", "nlo_hi", "pdg", "pdg_err", "derived",
                               "pull"],
                       ["configuration", "state", "process", "LO", "NLO", "± par", "± th",
                        "band lo", "band hi", "PDG", "±", "*", "pull"],
                       ["", "", "", ".3f", ".3f", ".3f", ".3f", ".3f", ".3f", ".3f", ".3f",
                        lambda v: "*" if v else "", ".2f"]) + "\n")
        rrows = width_ratio_rows(main_runs, phys)
        write_csv("width_ratios.csv", rrows)
        tables += ("### Width ratios: the wave functions without normalisation\n\n" + md_table(
            [r for r in rrows if r["tag"] in wtags],
            ["label", "ratio", "model_lo", "model_nlo", "pdg", "pdg_err"],
            ["configuration", "ratio", "model LO", "model NLO", "PDG", "±"],
            ["", "", ".3f", ".3f", ".3f", ".3f"]) + "\n")
        fig_widths(main_runs, wtags)

    # ---- validity audit: physics against model-independent references,
    # and the wave functions rebuilt from their logged parameters
    import validity
    sp = validity.spin_structure_rows(main_runs, phys)
    write_csv("spin_structure.csv", sp)
    jr = load("joint")
    qrows = validity.quark_mass_rows(main_runs, jr, phys)
    srows = validity.string_tension_rows(main_runs, jr, phys)
    ref_tag = phys[0]
    wrows = validity.wavefunction_check(main_runs[ref_tag])
    write_csv("wavefunction_check.csv", wrows)
    mix, own = validity.mixing_angle(main_runs[ref_tag])
    fig_wavefunctions(main_runs[ref_tag], describe(ref_tag))
    tables += ("## Validity audit\n\n"
               "Physics against references the fits never saw, and the wave functions "
               "rebuilt from their logged parameters alone (best model of the physical "
               f"set, {plain(ref_tag)}). Failures are kept.\n\n" +
               "".join(f"- {line}\n" for line in validity.summary(sp, qrows, srows, wrows,
                                                                   mix, own)) +
               "\nSpin structure per model (MeV; ρ = (M₂ − M₁)/(M₁ − M₀); cog − h: the 1P "
               "centre of gravity minus the ¹P₁ state, zero for a contact spin-spin force "
               "at first order):\n\n" + md_table(
                   sp, ["label", "system", "spin", "hf1_model", "hf1_model_err", "hf1_pdg",
                        "hf2_model", "hf2_model_err", "hf2_pdg", "rho_model", "rho_model_err",
                        "rho_pdg", "cog_h_model", "cog_h_model_err", "cog_h_pdg"],
                   ["configuration", "system", "spin", "1S hf", "±", "PDG", "2S hf", "±", "PDG",
                    "ρ", "±", "PDG", "cog − h", "±", "PDG"],
                   ["", "", "", ".1f", ".1f", ".1f", ".1f", ".1f", ".1f", ".2f", ".2f", ".2f",
                    ".1f", ".1f", ".2f"]) +
               "\n(± = parametric error of the model value, from the prediction covariance; "
               "the theory error of a splitting is not included)\n\n")

    # ---- studies: limits of the potential description, state-dependent
    # floors, semi-relativistic kinetic energy
    phys_runs = {t: main_runs[t] for t in phys}
    lim_tags = [t for t in ("A-apt-coff-smeared", "C-const-coff-s", "D-const-coff-s")
                if t in phys_runs]
    lrows = limits_rows(phys_runs, phys)
    if lrows:
        write_csv("limits_states.csv", lrows)
        fig_limits([r for r in lrows if r["tag"] in lim_tags], lim_tags,
                   r_break=R_STRING_BREAK_FM)
        cols = ["system", "bin", "states", "fitted", "n_below", "rms_below_MeV",
                "n_above", "rms_above_MeV", "rms_floor_MeV"]
        fmts = ["", "", "d", "d", "d", ".1f", "d", ".1f", ".1f"]

        def heads(var):
            return ["system", var, "states", "fitted", "below thr.", "RMS (MeV)",
                    "above thr.", "RMS (MeV)", "floor (MeV)"]
        tables += (f"## Limits of the potential description (physical set, {len(phys)} models)\n\n"
                   "Out-of-sample model − PDG (exact LOO miss for a fitted state, the "
                   "prediction otherwise), split at the open-flavour threshold. RMS per "
                   "model, median over models; state counts are distinct states.\n\n"
                   "In bins of the rms radius:\n\n" + md_table(
                       binned_residuals(lrows, "r_rms_fm", [0, 0.4, 0.6, 0.8, 1.0, 1.25, 2.0]),
                       cols, heads("r_rms (fm)"), fmts) +
                   "\nIn bins of ⟨v²⟩:\n\n" + md_table(
                       binned_residuals(lrows, "v2", [0, 0.1, 0.2, 0.25, 0.3, 0.5]),
                       cols, heads("⟨v²⟩"), fmts) +
                   f"\n(lattice string breaking: r_c ≈ {R_STRING_BREAK_FM} fm, "
                   "Bali et al., PRD 71, 114513 (2005))\n\n")
    shape_runs = load("shape")
    if shape_runs:
        frs = floor_study_rows(main_runs, load("cv"), shape_runs, phys, SHAPES)
        write_csv("floor_shapes.csv", frs)
        tables += ("## State-dependent theory floor (physical set)\n\n"
                   "Shapes from the reference model, frozen and shared by every model; "
                   "selection by LOO and LMO predictive density of the fitted states "
                   "(fixed in advance). ↑/↓: models whose paired difference from the "
                   "uniform floor exceeds +2 / −2 standard errors. The rms pull above "
                   "threshold is reported only.\n\n"
                   + md_table(floor_study_summary(frs),
                              ["shape", "models", "elpd_loo", "loo_up", "loo_down",
                               "elpd_lmo", "lmo_up", "lmo_down", "above_pull_rms_median"],
                              ["floor shape", "models", "Σ elpd LOO", "↑", "↓",
                               "Σ elpd LMO", "↑", "↓", "median rms pull above threshold"],
                              ["", "d", ".1f", "d", "d", ".1f", "d", "d", ".2f"]) + "\n")
    rel_runs = load("rel")
    if rel_runs:
        rrows = rel_rows(main_runs, load("cv"), rel_runs)
        write_csv("semirelativistic.csv", rrows)
        tables += ("## Semi-relativistic kinetic energy √(p²+m²) − m (A/B, running couplings)\n\n"
                   + md_table(rrows, ["label", "nr_ccbar_floor", "rel_ccbar_floor",
                                      "d_ccbar", "se_ccbar", "nr_bbar_floor",
                                      "rel_bbar_floor", "d_bbar", "se_bbar",
                                      "nr_ccbar_pred", "rel_ccbar_pred", "nr_bbar_pred",
                                      "rel_bbar_pred", "nr_bc", "rel_bc", "rel_ccbar_n",
                                      "rel_bbar_n"],
                              ["configuration", "σ_th cc̄ NR", "rel", "Δelpd cc̄", "±",
                               "σ_th bb̄ NR", "rel", "Δelpd bb̄", "±", "above cc̄ NR",
                               "rel", "above bb̄ NR", "rel", "B_c NR", "rel", "n cc̄",
                               "n bb̄"],
                              ["", ".1f", ".1f", ".1f", ".1f", ".1f", ".1f", ".1f", ".1f",
                               ".1f", ".1f", ".1f", ".1f", ".0f", ".0f", "d", "d"]) +
                   "\n(MeV; Δelpd: LOO + LMO predictive density, rel − NR, summed per "
                   "state, with its standard error over states; n: nested basis size "
                   "meeting 0.1 MeV)\n\n")

    nrows = nested_rows(load("nested"))
    if nrows:
        write_csv("nested_convergence.csv", nrows)
        bad = [r for r in nrows if r["chosen"] == 0 or r["chosen"] > 37]
        tables += ("## Basis size: nested convergence (all configurations)\n\n"
                   f"{len(nrows)} system fits. Largest change of any state, fitted or "
                   "predicted, at each doubling (MeV), worst case over all fits:\n\n")
        for n0, n1 in ((10, 19), (19, 37), (37, 73), (73, 145)):
            k = max(nrows, key=lambda r: r.get(f"d{n0}", 0))
            tables += (f"- n = {n0} → {n1}: {k[f'd{n0}']:.4g} MeV "
                       f"({k['label']}, {k['system']}, {k[f'w{n0}']})\n")
        tables += (f"\nn = 37 meets the 0.1 MeV target in "
                   f"{len(nrows) - len(bad)} of {len(nrows)} fits"
                   + (": not in " + ", ".join(f"{r['label']} {r['system']}"
                                              for r in bad) if bad else "")
                   + f". Largest rise of a level along the sequence (dropped "
                   f"near-dependent functions): "
                   f"{max(r['rise_MeV'] for r in nrows):.4f} MeV.\n\n")

    (REPORT / "tables.md").write_text(tables)
    phys = [r for r in rows if r["tag"] in phys_tags(rows)]
    print(f"{len(rows)} configurations; lowest floor (physical set: c = 0, "
          f"spin on, non-degenerate, distinct): {phys[0]['label'] if phys else '-'}")


if __name__ == "__main__":
    main()
