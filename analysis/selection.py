"""The pre-registered model-selection rule (fixed before the final sweeps; see
docs/EXPLANATION.md, "Model selection"). A model wins only if it

  Gate 1 - reproduces the observables: no rejection, at alpha = 0.05 after
           Holm's correction over its four tests, among
             masses c c-bar, masses b b-bar: Pearson chi^2 at the declared
               mass error (Paule-Mandel pooled over the physical set, refit);
             decays c c-bar, decays b b-bar: the decay benchmark, Pearson
               chi^2 of the measured widths of the fitted states against the
               same mass-only fit (V = V_exp + V_theory + V_par);
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
the model with the fewest failed checks (ties: AIC)."""
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


def select(main_runs, floor2_runs, nested_runs, tags, pooled_mev):
    """One row per model: every check, its value and verdict; then the
    winner (or None) and the best compromise."""
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
                pv[f"decay_{sname}"] = db["p"]
                row[f"decay_{sname}_chi2n"] = db["chi2"] / db["n"]
                chi2m += db["chi2"]
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
