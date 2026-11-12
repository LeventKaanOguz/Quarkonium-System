"""The pre-registered model-selection rule (fixed before the final sweeps; see
docs/EXPLANATION.md, "Model selection"), with one amendment made after the
results were seen (below). A model wins only if it

  Gate 1 - reproduces the observables: no rejection, at alpha = 0.05 after
           Holm's correction over its four tests, among
             masses c c-bar, masses b b-bar: Pearson chi^2 at the declared
               mass error (Paule-Mandel pooled over the physical set, refit);
             decays c c-bar, decays b b-bar: the decay benchmark, Pearson
               chi^2 of the measured widths of the fitted states against the
               same mass-only fit (V = V_exp + V_theory + V_par); amended:
               in both the linear form (as registered) and the log form
               (below), each with Holm over its four tests;
  Gate 2 - is physically sound, every check at its stated threshold:
             (a) non-degenerate fit, condition number < 1e10;
             (b) basis converged: nested size meeting 0.1 MeV <= n of the fits;
             (c) radial nodes = radial quantum number for every state;
             (d) m_b - m_c inside the heavy-quark window within 2 sigma;
             (e) sqrt(b) within 0.44 +- 0.04 GeV (lattice) within 2 sigma, both
                 systems;
             (f) 1P centre of gravity - h within 2 sigma of the PDG (sigma: the
                 declared mass error and the parametric error, propagated);
             (g) <v^2> < 0.5 for every fitted state;
and, among the models passing both gates, has the lowest
  AIC = chi2_masses(declared error) + chi2_decays + 2 p   (p: free parameters
  of both systems). If none passes, there is no winner; the best compromise is
the model with the fewest failed checks (ties: AIC).

Amendment (v1.0.0, after the final sweeps): the pre-registered decay test
was linear in Gamma, with a theory error proportional to the prediction,
sigma_th = f Gamma_NLO (scale band and <v^2> Gamma_NLO). Its pull
(Gamma_PDG - Gamma_NLO) / (f Gamma_NLO) is bounded below by -1/f (about -2 at
f = 0.5) however large the overprediction, so the test could not reject an
overprediction. The log form grades r = ln(Gamma_PDG / Gamma_NLO), with
V = D^-1 (V_th + V_par) D^-1 + diag(s_exp / Gamma_PDG)^2, D = diag(Gamma_NLO),
symmetric in over- and underprediction (and so milder than the linear form
on underpredictions). A model must pass gate 1 in both forms
(select_both): the amendment can only add rejections, never a pass the
registered rule refused. Ranking and tiebreak: the registered (linear) AIC."""
import math

import numpy as np

import model_selection as ms
from common import plain, system
from validity import SPIN, quark_mass_window

ALPHA = 0.05
SQRT_SIGMA = (0.40, 0.48)  # GeV: lattice 0.44 GeV +- 10 % (scale setting)
V2_MAX = 0.5
COND_MAX = 1e10
N_BASIS = 37


def _cog_h(s, pooled_mev):
    """(model - PDG) of the 1P centre of gravity minus h, and its 2-sigma
    tolerance (declared mass error and parametric covariance)."""
    names = SPIN[s["system"]]["P"]
    idx = [next(i for i, x in enumerate(s["states"]) if x["label"] == n) for n in names]
    a = np.array([1 / 9, 3 / 9, 5 / 9, -1.0])
    mod = sum(c * s["states"][i]["model_GeV"] for c, i in zip(a, idx))
    pdg = sum(c * s["states"][i]["pdg_GeV"] for c, i in zip(a, idx))
    C = np.array(s["pred_cov_GeV2"])[np.ix_(idx, idx)]
    var = a @ C @ a + (1e-3 * pooled_mev) ** 2 * float(a @ a)
    return 1e3 * (mod - pdg), 2e3 * math.sqrt(var)


def _width_key(w):
    """The one-loop-factor class of a width (the scale-band correlation)."""
    if w["process"] != "gamma gamma":
        return "ee"
    lab = w["label"]
    return "gg_P0" if lab.startswith("chi_") and "0(" in lab else (
        "gg_P2" if lab.startswith("chi_") else "gg_S")


def decay_log_chi2(s2):
    """The amended decay benchmark of a system: Pearson chi^2 of
    r = ln(Gamma_PDG / Gamma_NLO) over the benchmark widths. Relative errors
    as standard deviations of ln Gamma (first order): scale band (half of it,
    correlated within a one-loop-factor class), <v^2>, the parametric error
    (diagonal, as logged) and the experimental error. Returns
    (chi2, n, worst label, worst z) or None."""
    db = s2["analysis"].get("decay_benchmark")
    if not db:
        return None
    ws = [w for w in s2.get("widths", []) if w["label"] in db["labels"]
          and w.get("pdg_keV") and w["in_fit"]]
    if not ws:
        return None
    n = len(ws)
    r = np.array([math.log(w["pdg_keV"] / w["nlo_keV"]) for w in ws])
    sc = np.array([0.5 * abs(w["nlo_band_hi_keV"] - w["nlo_band_lo_keV"]) / w["nlo_keV"]
                   for w in ws])
    key = [_width_key(w) for w in ws]
    V = np.array([[sc[i] * sc[j] if key[i] == key[j] else 0.0 for j in range(n)]
                  for i in range(n)])
    V += np.diag([w["v2"] ** 2 + (w["nlo_err_par_keV"] / w["nlo_keV"]) ** 2 +
                  (w["pdg_err_keV"] / w["pdg_keV"]) ** 2 for w in ws])
    chi2 = float(r @ np.linalg.solve(V, r))
    z = r / np.sqrt(np.diag(V))
    k = int(np.argmax(abs(z)))
    return chi2, n, ws[k]["label"], float(z[k])


def select(main_runs, floor2_runs, nested_runs, tags, pooled_mev, decay_test="log"):
    """One row per model: every check, its value and verdict; then the
    winner (or None) and the best compromise. decay_test: "log" (amended)
    or "linear" (as pre-registered)."""
    lo, hi = quark_mass_window()
    rows = []
    for t in tags:
        r, f2 = main_runs[t], floor2_runs.get(t)
        row = {"tag": t, "label": plain(t), "fails": []}
        # ---- gate 1: four goodness-of-fit tests, Holm within the model
        pv, chi2m, n_par = {}, 0.0, 0
        for sname in ("ccbar", "bbar"):
            s2 = system(f2, sname) if f2 else None
            if s2 is not None and "chi2" in s2:
                pv[f"mass_{sname}"] = float(ms.stats.chi2.sf(s2["chi2"], s2["nu"]))
                chi2m += s2["chi2"]
                row[f"mass_{sname}_chi2nu"] = s2["chi2"] / s2["nu"]
            # the decays from the same fit as the masses (declared error)
            db = s2["analysis"].get("decay_benchmark") if s2 is not None else None
            if db:
                row[f"p_lin_decay_{sname}"] = db["p"]
                lg = decay_log_chi2(s2)
                if lg:
                    c2, n, wl, wz = lg
                    row[f"p_log_decay_{sname}"] = float(ms.stats.chi2.sf(c2, n))
                    row[f"decay_{sname}_worst"], row[f"decay_{sname}_worst_z"] = wl, wz
                if decay_test == "linear":
                    c2, n = db["chi2"], db["n"]
                elif lg:
                    c2, n = lg[0], lg[1]
                else:
                    c2 = None
                if c2 is not None:
                    pv[f"decay_{sname}"] = float(ms.stats.chi2.sf(c2, n))
                    row[f"decay_{sname}_chi2n"] = c2 / n
                    chi2m += c2
            n_par += system(r, sname)["n_free"]
        row.update({f"p_{k}": v for k, v in pv.items()})
        if len(pv) < 4:
            row["fails"].append("tests missing")
        adj = dict(zip(pv, ms.holm(list(pv.values())))) if pv else {}
        for k, v in adj.items():
            if v < ALPHA:
                row["fails"].append(f"{k} (Holm p = {v:.2g})")
        row["aic"] = chi2m + 2 * n_par
        # ---- gate 2: physics
        for sname in ("ccbar", "bbar"):
            s = system(r, sname)
            an = s["analysis"]
            cond = an.get("condition") or math.inf
            if cond >= COND_MAX:
                row["fails"].append(f"{sname} degenerate ({cond:.1e})")
            nj = (system(nested_runs[t], sname)["analysis"].get("nested")
                  if t in nested_runs else None)
            if nj and (not nj["chosen"] or nj["chosen"] > N_BASIS):
                row["fails"].append(f"{sname} basis not converged at n = {N_BASIS}")
            bad = [x["label"] for x in s["states"] if x.get("wavefunction")
                   and x["wavefunction"]["nodes"] != x["n"]]
            if bad:
                row["fails"].append(f"{sname} nodes: {', '.join(bad)}")
            v2 = max(x["v2"] for x in s["states"] if x["in_fit"])
            row[f"v2max_{sname}"] = v2
            if v2 >= V2_MAX:
                row["fails"].append(f"{sname} <v^2> = {v2:.2f}")
            if "b" in an["names"]:
                i = an["names"].index("b")
                sb = math.sqrt(an["value"][i])
                esb = an["sigma"][i] / (2 * sb)
                row[f"sqrtb_{sname}"], row[f"sqrtb_{sname}_err"] = sb, esb
                if sb + 2 * esb < SQRT_SIGMA[0] or sb - 2 * esb > SQRT_SIGMA[1]:
                    row["fails"].append(f"{sname} sqrt(b) = {sb:.3f} +- {esb:.3f} GeV")
            d, tol = _cog_h(s, pooled_mev[sname])
            row[f"cogh_{sname}"], row[f"cogh_{sname}_tol"] = d, tol
            if abs(d) > tol:
                row["fails"].append(f"{sname} cog - h off by {d:.1f} MeV (> {tol:.1f})")
        sc, sb_ = system(r, "ccbar"), system(r, "bbar")
        dm = sb_["theta"][0] - sc["theta"][0]
        edm = math.hypot(sb_["analysis"]["sigma"][0], sc["analysis"]["sigma"][0])
        row["dm"], row["dm_err"] = dm, edm
        if dm + 2 * edm < lo or dm - 2 * edm > hi:
            row["fails"].append(f"m_b - m_c = {dm:.3f} +- {edm:.3f} GeV outside "
                                f"[{lo:.3f}, {hi:.3f}]")
        row["n_fail"] = len(row["fails"])
        row["verdict"] = "pass" if not row["fails"] else "fail"
        rows.append(row)
    rows.sort(key=lambda x: (x["n_fail"], x["aic"]))
    winners = [x for x in rows if x["verdict"] == "pass"]
    return rows, (winners[0] if winners else None), rows[0] if rows else None


def select_both(main_runs, floor2_runs, nested_runs, tags, pooled_mev):
    """The amended rule: gate 1 in both decay forms. Rows carry the log-form
    p-values, the linear ones (p_lin_*), the three verdicts and the union of
    the failed checks; ranked by the number of failures, then the registered
    (linear) AIC. Returns (rows, winner or None, best compromise)."""
    rl, _, _ = select(main_runs, floor2_runs, nested_runs, tags, pooled_mev, "linear")
    rg, _, _ = select(main_runs, floor2_runs, nested_runs, tags, pooled_mev, "log")
    lin = {r["tag"]: r for r in rl}
    rows = []
    for r in rg:
        a = lin[r["tag"]]
        # the mass tests are the same in both forms (their Holm p differs
        # only through the decay p-values): keep the registered one
        fails = list(dict.fromkeys(
            [f"{f} [linear]" if f.startswith("decay_") else f for f in a["fails"]] +
            [f"{f} [log]" for f in r["fails"] if f.startswith("decay_")] +
            [f for f in r["fails"] if not f.startswith(("decay_", "mass_"))]))
        fails += [f for f in r["fails"] if f.startswith("mass_")
                  and not any(g.startswith(f.split(" (")[0]) for g in fails)]
        rows.append({**r, "fails": fails, "n_fail": len(fails),
                     "verdict_log": r["verdict"], "verdict_linear": a["verdict"],
                     "verdict": "pass" if not fails else "fail",
                     "aic_log": r["aic"], "aic": a["aic"]})
    rows.sort(key=lambda x: (x["n_fail"], x["aic"]))
    winners = [x for x in rows if x["verdict"] == "pass"]
    return rows, (winners[0] if winners else None), rows[0] if rows else None
