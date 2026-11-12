"""Analysis tables: from the sweep JSONs to rows and Markdown (no ROOT)."""
import math

import numpy as np
from scipy import optimize, stats

import model_selection as ms
from common import md_table, plain, rms, system
from sweep import duplicate_fit, nulls


def summarize(runs):
    """One row per configuration: the numbers every ranking uses."""
    rows = []
    for tag, r in runs.items():
        row = {"tag": tag, "label": plain(tag), "model": r["config"]["model"],
               "running": r["config"]["running"],
               "c_free": r["config"]["c_free"], "spin": r["config"]["spin"],
               "fit_mu_s": r["config"]["fit_mu_s"]}
        ok = True
        for sname in ("ccbar", "bbar"):
            s = system(r, sname)
            if s is None or "sigma_th_GeV" not in s:
                ok = False
                break
            an = s["analysis"]
            miss = an.get("loo_exact_miss_GeV")
            pred = [1e3 * (x["model_GeV"] - x["pdg_GeV"]) for x in s["states"]
                    if not x["in_fit"]]
            fit = [1e3 * (x["model_GeV"] - x["pdg_GeV"]) for x in s["states"]
                   if x["in_fit"]]
            row.update({
                f"{sname}_floor": 1e3 * s["sigma_th_GeV"],
                f"{sname}_floor_lo": 1e3 * (an.get("sigma_th_ci68_GeV") or [math.nan])[0],
                f"{sname}_floor_hi": 1e3 * (an.get("sigma_th_ci68_GeV") or [0, math.nan])[1],
                f"{sname}_chi2": s["chi2"],
                f"{sname}_N": s["n_data"], f"{sname}_p": s["n_free"],
                f"{sname}_cond": (an.get("condition") if an.get("condition")
                                  is not None else float("inf")),
                f"{sname}_conv": an.get("converged", 0),
                f"{sname}_fit_rms": rms(fit),
                f"{sname}_loo_rms": 1e3 * rms(miss) if miss else float("nan"),
                f"{sname}_loo_fail": sum(1 for x in (miss or []) if x is None),
                f"{sname}_pred_rms": rms(pred),
                f"{sname}_start0": 1e3 * an["multistart_scores"][0]
                if "multistart_scores" in an and an["multistart_scores"][0]
                is not None else float("nan"),
                f"{sname}_best_start": an.get("multistart_best", 0),
            })
        if not ok:
            continue
        bc = system(r, "bc")
        bcm = bc["states"][0]["model_GeV"] if bc else None
        # None: the interpolated B_c Hamiltonian collapsed (see run.cpp
        # warning). With c free the B_c is not a prediction: c is degenerate
        # with the quark mass, so averaging it across systems is meaningless.
        row["bc_diff"] = (1e3 * (bcm - bc["states"][0]["pdg_GeV"])
                          if bcm is not None and not row["c_free"]
                          else float("nan"))
        row["degenerate"] = int(max(row["ccbar_cond"], row["bbar_cond"]) > 1e10)
        row["score"] = math.hypot(row["ccbar_floor"], row["bbar_floor"])
        rows.append(row)
    rows.sort(key=lambda x: x["score"])
    return rows


def overfitting_rows(main_runs, cv_runs, tags):
    """In-sample RMS, exact LOO RMS, leave-multiplet-out RMS, and the LOO
    inflation expected for a correctly specified (linearised) model:
    E[LOO^2]/E[res^2] = sum 1/(1-h_i) / (N - p) (leverage h_i)."""
    rows = []
    for t in tags:
        r = main_runs[t]
        row = {"tag": t, "label": plain(t)}
        for sname in ("ccbar", "bbar"):
            s = system(r, sname)
            an = s["analysis"]
            fit = [1e3 * (x["model_GeV"] - x["pdg_GeV"]) for x in s["states"]
                   if x["in_fit"]]
            h = [x for x in (an.get("leverage") or []) if x is not None]
            miss = [1e3 * m for m in (an.get("loo_exact_miss_GeV") or [])
                    if m is not None]
            N, p = s["n_data"], s["n_free"]
            row[f"{sname}_fit_rms"] = rms(fit)
            row[f"{sname}_loo_rms"] = rms(miss) if miss else float("nan")
            row[f"{sname}_inflation"] = (row[f"{sname}_loo_rms"] / row[f"{sname}_fit_rms"]
                                         if row[f"{sname}_fit_rms"] > 0 else float("nan"))
            row[f"{sname}_inflation_expected"] = (
                math.sqrt(sum(1 / (1 - x) for x in h if x < 1) / (N - p))
                if h and N > p else float("nan"))
            ie = row[f"{sname}_inflation_expected"]
            row[f"{sname}_inflation_ratio"] = (row[f"{sname}_inflation"] / ie
                                               if ie and math.isfinite(ie)
                                               else float("nan"))
            row[f"{sname}_eff_dof"] = sum(h) if h else float("nan")
            lm = []
            cv = cv_runs.get(t)
            if cv:
                for g in system(cv, sname)["analysis"].get("lmo", []):
                    lm += [1e3 * m for m in g["miss_GeV"] if m is not None]
            row[f"{sname}_lmo_rms"] = rms(lm) if lm else float("nan")
        rows.append(row)
    return rows


def identifiability_rows(main_runs, rows):
    """For every degenerate fit: the parameter combination the data do not fix,
    i.e. the leading eigenvector of the relative covariance
    diag(sigma/theta) rho diag(sigma/theta) (c and other zero-centred
    parameters use sigma itself)."""
    out = []
    for r in rows:
        if not r["degenerate"]:
            continue
        for sname in ("ccbar", "bbar"):
            an = system(main_runs[r["tag"]], sname)["analysis"]
            cond = an.get("condition") or float("inf")
            if cond < 1e10:
                continue
            v = np.array(an["value"], float)
            sg = np.array(an["sigma"], float)
            rel = np.array([sg[k] if abs(v[k]) < 1e-9 or an["names"][k] == "c"
                            else sg[k] / abs(v[k]) for k in range(len(v))])
            if not np.all(np.isfinite(rel)):
                out.append({"label": r["label"], "system": sname,
                            "combination": "covariance not finite",
                            "rel_sigma": float("inf")})
                continue
            corr = np.array([[np.nan if x is None else x for x in row]
                             for row in an["corr"]], float)
            C = np.outer(rel, rel) * corr
            if not np.all(np.isfinite(C)):
                out.append({"label": r["label"], "system": sname,
                            "combination": "covariance not finite",
                            "rel_sigma": float("inf")})
                continue
            w, V = np.linalg.eigh(C)
            e = V[:, -1]
            e = e / np.max(np.abs(e))
            comb = " ".join(f"{e[k]:+.2f}·{n}" for k, n in enumerate(an["names"])
                            if abs(e[k]) > 0.1)
            out.append({"label": r["label"], "system": sname, "combination": comb,
                        "rel_sigma": float(math.sqrt(max(w[-1], 0)))})
    return out


def universality_rows(main_runs, tags):
    """Flavour independence: the same potential parameter from the c c-bar and
    the b b-bar fit; z = (theta_c - theta_b) / sqrt(sigma_c^2 + sigma_b^2)."""
    out = []
    for t in tags:
        a = system(main_runs[t], "ccbar")["analysis"]
        b = system(main_runs[t], "bbar")["analysis"]
        row = {"tag": t, "label": plain(t)}
        for n in a["names"]:
            if n == "m1" or n not in b["names"]:
                continue
            i, j = a["names"].index(n), b["names"].index(n)
            d = a["value"][i] - b["value"][j]
            e = math.hypot(a["sigma"][i], b["sigma"][j])
            row[f"{n}_c"], row[f"{n}_b"] = a["value"][i], b["value"][j]
            row[f"{n}_z"] = d / e if e > 0 else float("nan")
        out.append(row)
    return out


def residual_structure_rows(main_runs, tags):
    """Underfitting shows as structure in the residuals: runs test on the signs
    of model - PDG ordered by mass (fitted states; and all states, where the
    above-threshold ones are expected to drift), Shapiro-Wilk on the pooled
    leave-one-out pulls of both systems."""
    out = []
    for t in tags:
        row = {"tag": t, "label": plain(t)}
        pulls = []
        for sname in ("ccbar", "bbar"):
            st = sorted(system(main_runs[t], sname)["states"], key=lambda x: x["pdg_GeV"])
            fit = [x["model_GeV"] - x["pdg_GeV"] for x in st if x["in_fit"]]
            allr = [x["model_GeV"] - x["pdg_GeV"] for x in st]
            row[f"{sname}_runs_p"], row[f"{sname}_runs"] = ms.runs_test(fit)
            row[f"{sname}_runs_all_p"], _ = ms.runs_test(allr)
            pulls += [x["pull"] for x in st if x["in_fit"] and x["pull"] is not None]
        row["pull_normal_p"] = ms.normality_p(pulls)
        row["pull_rms"] = rms(pulls)
        out.append(row)
    return out


def cv_rows(cv_runs, tags):
    """Repeated K-fold and leave-multiplet-out summaries from the cv sweep."""
    out = []
    for t in tags:
        r = cv_runs.get(t)
        if not r:
            continue
        row = {"tag": t, "label": plain(t)}
        for sname in ("ccbar", "bbar"):
            an = system(r, sname)["analysis"]
            kf = an.get("kfold") or {}
            row[f"{sname}_kf_train"] = 1e3 * (kf.get("train_rms_GeV") or float("nan"))
            row[f"{sname}_kf_test"] = 1e3 * (kf.get("test_rms_GeV") or float("nan"))
            row[f"{sname}_kf_median"] = 1e3 * (kf.get("test_median_GeV") or float("nan"))
        out.append(row)
    return out


def bgrid_rows(bgrid_runs):
    rows = []
    for key, r in bgrid_runs.items():
        tag, rest = key.split("-r", 1)
        r1, rm = rest.split("-R")
        for sname in ("ccbar", "bbar"):
            s = system(r, sname)
            if s is None or "sigma_th_GeV" not in s:
                continue
            an = s["analysis"]
            rows.append({"tag": tag, "label": plain(tag), "system": sname,
                         "r1": float(r1), "rmax": float(rm),
                         "floor": 1e3 * s["sigma_th_GeV"],
                         "trace": an.get("basis_trace_GeV", float("nan")),
                         "num_all": 1e3 * max(x["sigma_num_GeV"] for x in s["states"])})
    for tag in {r["tag"] for r in rows}:
        for sname in ("ccbar", "bbar"):
            sub = [r for r in rows if r["tag"] == tag and r["system"] == sname]
            tmin = min((r["trace"] for r in sub if math.isfinite(r["trace"])),
                       default=float("nan"))
            for r in sub:
                r["trace_excess_MeV"] = 1e3 * (r["trace"] - tmin)
    return rows


def phys_tags(rows):
    """The physical set: c = 0, spin terms on, Fisher matrix not degenerate,
    one configuration per distinct fit (see sweep.duplicate_fit)."""
    tags = [r["tag"] for r in rows if not r["degenerate"] and r["spin"] != "none"
            and not r["c_free"]]
    return [t for t in tags if not duplicate_fit(t, tags)]


def merge_systems(runs):
    """floor2 runs are one JSON per (configuration, system): tag-ccbar,
    tag-bbar. Merge them into one run per configuration."""
    out = {}
    for key, r in runs.items():
        tag, sname = key.rsplit("-", 1)
        out.setdefault(tag, {"systems": [], "config": r["config"]})
        out[tag]["systems"] += [x for x in r["systems"] if x["system"] == sname]
    return out


def floor_analysis(runs_floor):
    """Fixed floor: chi2 goodness of fit, AIC, AICc, BIC and the Laplace
    evidence per system and summed; weights of each criterion."""
    rows = []
    for tag, r in runs_floor.items():
        row = {"tag": tag, "label": plain(tag)}
        ok = True
        for sname in ("ccbar", "bbar"):
            s = system(r, sname)
            if s is None or "chi2" not in s:
                ok = False
                break
            N, p, chi2 = s["n_data"], s["n_free"], s["chi2"]
            an = s["analysis"]
            cond = an.get("condition") or float("inf")
            # Laplace needs a proper (non-degenerate) posterior
            lnz = (ms.laplace_log_evidence(chi2, an["names"], an["value"],
                                           an["sigma"], an["corr"])
                   if cond < 1e10 else float("nan"))
            lnz2 = (ms.laplace_log_evidence(chi2, an["names"], an["value"],
                                            an["sigma"], an["corr"], decades=2.0,
                                            c_range=2.0)
                    if cond < 1e10 else float("nan"))
            row.update({f"{sname}_chi2": chi2, f"{sname}_p": p,
                        f"{sname}_N": N,
                        f"{sname}_gof_p": chi2_sf(chi2, N - p),
                        f"{sname}_aic": ms.aic(chi2, p),
                        f"{sname}_aicc": ms.aicc(chi2, p, N),
                        f"{sname}_bic": ms.bic(chi2, p, N),
                        f"{sname}_lnz": lnz, f"{sname}_lnz2": lnz2})
        if ok:
            for k in ("aic", "aicc", "bic", "lnz", "lnz2"):
                row[k] = row[f"ccbar_{k}"] + row[f"bbar_{k}"]
            rows.append(row)
    rows.sort(key=lambda x: x["aicc"])
    for key, wkey, f in (("aic", "akaike_weight", ms.ic_weights),
                         ("aicc", "aicc_weight", ms.ic_weights),
                         ("bic", "bic_weight", ms.ic_weights),
                         ("lnz", "evidence_weight", ms.evidence_weights),
                         ("lnz2", "evidence2_weight", ms.evidence_weights)):
        if rows:
            for r, w in zip(rows, f([r[key] for r in rows])):
                r[wkey] = float(w)
    return rows


def chi2_sf(x, k):
    """P(chi2_k >= x)."""
    return float(stats.chi2.sf(max(x, 0.0), k)) if k > 0 else float("nan")


def nested_tests(frows):
    by = {r["tag"]: r for r in frows}
    pairs = []
    for t in by:
        for null in nulls(t):
            if "con" in t.split("-") and "coff" in null.split("-"):
                pairs.append((null, t, "c free", False))
            elif t.endswith("-fitmu") and null[0] in "AC":
                # mu_s -> 0 is the unscreened model, on the boundary
                pairs.append((null, t, "screening mass fitted (vs none)", True))
            elif t.endswith("-fitmu"):
                pairs.append((null, t, "screening mass fitted (vs fixed)", False))
            else:
                pairs.append((null, t, "free GI amplitudes", False))
    out = []
    for null, alt, what, boundary in pairs:
        if null not in by or alt not in by:
            continue
        for sname in ("ccbar", "bbar"):
            a, b = by[null], by[alt]
            dp = b[f"{sname}_p"] - a[f"{sname}_p"]
            d = a[f"{sname}_chi2"] - b[f"{sname}_chi2"]
            if dp <= 0:
                continue
            # Wilks: dchi2 ~ chi2_dp under the null; Chernoff: half of it when
            # the null value sits on the boundary (mu_s -> 0). dchi2 < 0 would
            # mean the extended fit missed the nested minimum.
            pv = chi2_sf(d, dp) * (0.5 if boundary else 1.0)
            out.append({"null": plain(null), "alt": plain(alt), "system": sname,
                        "test": what, "dchi2": d, "dp": dp, "p_value": pv,
                        "verdict": ("fit missed nested minimum" if d < -0.01 else
                                    "extension significant (p < 0.05)" if pv < 0.05
                                    else "not significant")})
    if out:
        pv = [t["p_value"] for t in out]
        for t, h, b in zip(out, ms.holm(pv), ms.benjamini_hochberg(pv)):
            t["p_holm"], t["p_bh"] = float(h), float(b)
    return out


def minuit_analysis(main_runs, lm_runs, mn_runs):
    """sigma_th from: LM one start, LM nine starts (`lm`), Minuit2 one start
    (`minuit`) and the final seeded fit (`main`)."""
    rows = []
    nan = float("nan")
    for tag, r in mn_runs.items():
        if tag not in main_runs or tag not in lm_runs:
            continue
        for sname in ("ccbar", "bbar"):
            s, m = system(r, sname), system(main_runs[tag], sname)
            l = system(lm_runs[tag], sname)
            if s is None or m is None or l is None or "sigma_th_GeV" not in s \
                    or "sigma_th_GeV" not in m or "sigma_th_GeV" not in l:
                continue
            an, al = s["analysis"], l["analysis"]
            cc = an.get("cross_check", {})
            lm0 = al["multistart_scores"][0] if "multistart_scores" in al else None
            lmb = 1e3 * l["sigma_th_GeV"]
            mnf = 1e3 * s["sigma_th_GeV"]
            fin = 1e3 * m["sigma_th_GeV"]
            rows.append({
                "tag": tag, "label": plain(tag), "system": sname,
                "lm_single_floor": 1e3 * lm0 if lm0 else nan,
                "lm_multistart_floor": lmb,
                "minuit_floor": mnf,
                "final_floor": fin,
                "lm_ms_excess": lmb / fin - 1,
                "minuit_excess": mnf / fin - 1,
                "minuit_valid": cc.get("valid", 0),
                "minuit_edm": cc.get("edm", float("nan")),
                "dchi2_restart": cc.get("dchi2_vs_lm", float("nan")),
                # Minos is skipped (null) when Migrad's minimum is not valid
                "max_minos_asym": max([abs(h + l) / (h - l)
                                       for l, h in zip(cc.get("minos_lo") or [],
                                                       cc.get("minos_hi") or [])
                                       if l is not None and h is not None
                                       and h - l > 0] or [float("nan")]),
            })
    return rows


def error_methods(stats_runs, mn_runs):
    """Linear, Minos, bootstrap and MCMC errors of each fitted parameter."""
    rows = []
    for tag, r in stats_runs.items():
        for sname in ("ccbar", "bbar"):
            s = system(r, sname)
            if s is None or "analysis" not in s:
                continue
            an = s["analysis"]
            names, val, sig = an["names"], an["value"], an["sigma"]
            bs = an.get("bootstrap", {}).get("summary", {})
            mc = an.get("mcmc", {}).get("summary", {})
            cc = {}
            if tag in mn_runs:
                sm = system(mn_runs[tag], sname)
                cc = sm["analysis"].get("cross_check", {}) if sm else {}
            for k, nm in enumerate(names):
                rows.append({
                    "tag": tag, "label": plain(tag), "system": sname,
                    "param": nm, "value": val[k], "linear": sig[k],
                    "minos_lo": cc["minos_lo"][k] if cc.get("valid") else float("nan"),
                    "minos_hi": cc["minos_hi"][k] if cc.get("valid") else float("nan"),
                    "bootstrap": (bs["std"][k] if bs.get("std") and bs["std"][k]
                                  else float("nan")),
                    "bootstrap_failed": an.get("bootstrap", {}).get("failed", 0),
                    "mcmc": mc.get("std", [float("nan")] * 9)[k],
                    "mcmc_rhat": an.get("mcmc", {}).get("rhat", [float("nan")] * 9)[k],
                    "jackknife": (an["jackknife_sigma"][k]
                                  if an.get("jackknife_sigma") and
                                  an["jackknife_sigma"][k] is not None
                                  else float("nan")),
                })
    return rows


def nested_rows(nested_runs):
    """Basis convergence on the nested sizes 10, 19, 37, 73, 145: the largest
    change of any state (fitted or predicted) at each doubling."""
    out = []
    for tag, r in nested_runs.items():
        for sname in ("ccbar", "bbar"):
            s = system(r, sname)
            nj = s["analysis"].get("nested") if s else None
            if not nj:
                continue
            row = {"tag": tag, "label": plain(tag), "system": sname,
                   "chosen": nj["chosen"],
                   "gap_MeV": 1e3 * nj["gap_to_largest_GeV"],
                   "rise_MeV": 1e3 * nj["max_increase_GeV"]}
            for n0, n1, d, w in zip(nj["n"], nj["n"][1:], nj["max_change_GeV"],
                                    nj["worst"]):
                row[f"d{int(n0)}"] = 1e3 * d
                row[f"w{int(n0)}"] = w
            out.append(row)
    return sorted(out, key=lambda r: (r["label"], r["system"]))


def check_basis(runs_by_sweep):
    """Every sweep that enters the comparison must use the default basis; lm
    and minuit (minimiser study) and bgrid (the basis plane) are exempt."""
    seen = {}
    for sweep, runs in runs_by_sweep.items():
        for tag, r in runs.items():
            c = r.get("config", {})
            key = (c.get("basis_n"), c.get("basis_r1"), c.get("basis_rmax"))
            seen.setdefault(key, set()).add(sweep)
    if len(seen) > 1:
        raise SystemExit("mixed bases in the analysed sweeps: " +
                         "; ".join(f"{k}: {sorted(v)}" for k, v in seen.items()))
    return next(iter(seen)) if seen else None


def held_out_points(run, sname, kind):
    """Held-out predictions of one system ('loo': exact leave-one-out; 'lmo':
    leave-multiplet-out): label, miss (GeV), predictive sd = sqrt((sigma_th
    w)^2 + s_exp^2 + sigma_par^2) with sigma_th and sigma_par from the refit
    without the held-out data, and nu_left, the degrees of freedom that fixed
    the refit's floor. sd is NaN where that floor is undetermined (nu_left = 0)
    or the refit failed."""
    s = system(run, sname)
    an, nu = s["analysis"], s["nu"]
    state = {x["label"]: x for x in s["states"]}

    def point(lab, m, th, sp, nu_left):
        x = state[lab]
        ok = None not in (m, th, sp)
        sd = (math.sqrt((th * x.get("floor_w", 1.0)) ** 2 + x["sigma_exp_GeV"] ** 2
                        + sp ** 2) if ok else math.nan)
        return {"system": sname, "kind": kind, "label": lab,
                "miss": math.nan if m is None else m, "sd": sd, "nu_left": nu_left}
    if kind == "loo":
        fitted = [x["label"] for x in s["states"] if x["in_fit"]]
        return [point(lab, m, th, sp, nu - 1) for lab, m, th, sp in zip(
            fitted, an["loo_exact_miss_GeV"], an["loo_sigma_th_GeV"],
            an["loo_sigma_par_GeV"])]
    return [point(lab, m, g["sigma_th_GeV"], sp, nu - len(g["labels"]))
            for g in an.get("lmo", [])
            for lab, m, sp in zip(g["labels"], g["miss_GeV"], g["sigma_par_GeV"])]


def scored(points):
    """The points with a finite miss and predictive sd, with t = miss / sd and
    their log predictive density."""
    out = [dict(q) for q in points if math.isfinite(q["miss"]) and math.isfinite(q["sd"])]
    for q in out:
        q["t"] = q["miss"] / q["sd"]
        q["logdens"] = float(ms.log_predictive_density([q["miss"]], [q["sd"]])[0])
    return out


def logdens_map(loo_run, lmo_run):
    """{(kind, system, label): log predictive density} of one configuration:
    LOO from `loo_run`, LMO from `lmo_run` (either may be None)."""
    out = {}
    for kind, run in (("loo", loo_run), ("lmo", lmo_run)):
        if run is None:
            continue
        for sname in ("ccbar", "bbar"):
            for q in scored(held_out_points(run, sname, kind)):
                out[(kind, sname, q["label"])] = q["logdens"]
    return out


def paired_elpd(a, b, kind=None, sname=None):
    """Sum and standard error of the pointwise difference b - a of two
    logdens_maps over the points both have: (delta, se, n)."""
    keys = [k for k in a if k in b and (kind is None or k[0] == kind)
            and (sname is None or k[1] == sname)]
    d = np.array([b[k] - a[k] for k in keys])
    if len(d) < 2:
        return math.nan, math.nan, len(d)
    return float(d.sum()), float(math.sqrt(len(d)) * d.std(ddof=1)), len(d)


def paired_elpd_by_state(a, b, sname):
    """As paired_elpd for one system, with the LOO and LMO differences of the
    same state summed first: the two are correlated, so the standard error is
    taken over states."""
    per = {}
    for k in a:
        if k in b and k[1] == sname:
            per[k[2]] = per.get(k[2], 0.0) + b[k] - a[k]
    d = np.array(list(per.values()))
    if len(d) < 2:
        return math.nan, math.nan, len(d)
    return float(d.sum()), float(math.sqrt(len(d)) * d.std(ddof=1)), len(d)


def expected_rms_t(nus):
    """sqrt(mean nu/(nu - 2)): the RMS of t_nu variables; NaN if any nu <= 2."""
    return (math.sqrt(float(np.mean([n / (n - 2) for n in nus])))
            if nus and min(nus) > 2 else math.nan)


def loo_logdens(run):
    """Pointwise LOO log predictive densities, c c-bar then b b-bar; None if
    any point is missing."""
    out = []
    for sname in ("ccbar", "bbar"):
        pts = held_out_points(run, sname, "loo")
        sc = scored(pts)
        if len(sc) != len(pts):
            return None
        out += [q["logdens"] for q in sc]
    return out


def predictive_analysis(main_runs, tags):
    """elpd_loo, pseudo-BMA and stacking weights over `tags`."""
    rows, L = [], []
    for t in tags:
        ld = loo_logdens(main_runs[t])
        if ld is None:
            continue
        n_c = sum(x["in_fit"] for x in system(main_runs[t], "ccbar")["states"])
        rows.append({"tag": t, "label": plain(t), "elpd_ccbar": sum(ld[:n_c]),
                     "elpd_bbar": sum(ld[n_c:]), "elpd": sum(ld)})
        L.append(ld)
    if rows:
        for r, w in zip(rows, ms.pseudo_bma_weights([r["elpd"] for r in rows])):
            r["pseudo_bma_weight"] = float(w)
        for r, w in zip(rows, ms.stacking_weights(np.array(L))):
            r["stacking_weight"] = float(w)
    rows.sort(key=lambda r: -r["elpd"])
    return rows


def averaged_predictions(main_runs, weight_sets):
    """Model-averaged predictions of every state that was not fitted (above
    threshold) and of the B_c, for each set of weights {tag: w}."""
    first = next(iter(weight_sets.values()))
    tags = list(first)
    ref = main_runs[tags[0]]
    rows = []
    for sname in ("ccbar", "bbar", "bc"):
        for k, x0 in enumerate(system(ref, sname)["states"]):
            if x0["in_fit"]:
                continue
            row = {"system": sname, "state": x0["label"], "pdg": 1e3 * x0["pdg_GeV"]}
            vals, sigs = {}, {}
            for t in tags:
                r = main_runs[t]
                x = system(r, sname)["states"][k]
                if sname == "bc":
                    th = 0.5 * (system(r, "ccbar")["sigma_th_GeV"] +
                                system(r, "bbar")["sigma_th_GeV"])
                else:
                    th = system(r, sname)["sigma_th_GeV"]
                m = x["model_GeV"]
                vals[t] = float("nan") if m is None else 1e3 * m
                sigs[t] = 1e3 * math.hypot(th, x.get("sigma_par_GeV") or 0.0,
                                           x.get("sigma_num_GeV") or 0.0)
            for wname, w in weight_sets.items():
                mean, unc, spread = ms.model_average(
                    [vals[t] for t in tags], [sigs[t] for t in tags],
                    [w.get(t, 0.0) for t in tags])
                row[f"{wname}_mean"] = mean
                row[f"{wname}_diff"] = mean - row["pdg"]
                row[f"{wname}_unc"] = unc
                row[f"{wname}_spread"] = spread
            rows.append(row)
    return rows


def criteria_agreement(scores):
    """Kendall tau between every pair of criteria; `scores`: name -> {tag:
    score}, lower = better."""
    names = list(scores)
    tags = sorted(set.intersection(*(set(v) for v in scores.values())))
    M = np.full((len(names), len(names)), np.nan)
    for i, a in enumerate(names):
        for j, b in enumerate(names):
            M[i, j] = ms.kendall_tau([scores[a][t] for t in tags],
                                     [scores[b][t] for t in tags])
    return names, M, tags


def ic_md(frows):
    return md_table(
        frows, ["label", "ccbar_chi2", "bbar_chi2", "ccbar_gof_p", "bbar_gof_p",
                "ccbar_p", "aic", "aicc", "bic", "lnz", "akaike_weight", "aicc_weight",
                "bic_weight", "evidence_weight"],
        ["configuration", "χ² cc̄", "χ² bb̄", "p(χ²) cc̄", "p(χ²) bb̄", "params", "AIC",
         "AICc", "BIC", "ln Z", "w AIC", "w AICc", "w BIC", "w Z"],
        ["", ".1f", ".1f", ".2f", ".2f", "d", ".1f", ".1f", ".1f", ".1f", ".3f", ".3f",
         ".3f", ".3f"]) + (
        "\n(sorted by AICc. At a fixed floor the variance is known, and for a linear "
        "Gaussian model AIC's 2p is then the exact optimism; AICc, derived for an "
        "estimated variance, over-penalises and is shown for reference. ln Z: Laplace "
        "evidence with flat priors over one decade per log-parameter and 1 GeV in c, up "
        "to a common constant; – = degenerate fit or posterior wider than the prior)\n\n")


def wilks_summary(tests, floor):
    n = len(tests)
    k = sum(t["p_value"] < 0.05 for t in tests)
    kh = sum(t["p_holm"] < 0.05 for t in tests)
    kb = sum(t["p_bh"] < 0.05 for t in tests)
    miss = sum(t["dchi2"] < -0.01 for t in tests)
    return (f"At the {floor} floor: {n} tests, {k} with p < 0.05 (≈ {0.05 * n:.0f} "
            f"expected by chance), {kh} after Holm, {kb} after Benjamini–Hochberg; "
            f"{miss} with Δχ² < 0.\n\n")


def error_summary(erows):
    """Ratio of each method to the linear error, per model: median and range."""
    out = "| configuration | Minos | bootstrap | MCMC | jackknife | max R̂ |\n|---|---|---|---|---|---|\n"
    by = {}
    for r in erows:
        by.setdefault(r["label"], []).append(r)
    for lab, rs in by.items():
        cells = []
        for key in ("minos", "bootstrap", "mcmc", "jackknife"):
            v = [(0.5 * (r["minos_hi"] - r["minos_lo"]) if key == "minos" else r[key])
                 / r["linear"] for r in rs if r["linear"] > 0]
            v = [x for x in v if math.isfinite(x)]
            cells.append(f"{np.median(v):.2f} [{min(v):.2f}, {max(v):.2f}]" if v else "–")
        rh = [r["mcmc_rhat"] for r in rs if math.isfinite(r["mcmc_rhat"])]
        out += f"| {lab} | " + " | ".join(cells) + f" | {max(rh):.2f} |\n" if rh else \
            f"| {lab} | " + " | ".join(cells) + " | – |\n"
    return "\nRatio to the linear error, median [min, max] over the parameters of both systems:\n\n" + out + "\n"


# ------------------------------------------------ goodness of fit --
def gof_rows(main_runs, floor0_runs, tags):
    """Per model and system: the conventional chi^2 (PDG errors only) at the
    floor fit and refitted, the floor with its confidence intervals (the chi^2
    test inverted), the largest parameter shift of the PDG-only refit in units
    of the floor fit's error, and the RMS of the standardised LOO misses."""
    out = []
    for t in tags:
        for sname in ("ccbar", "bbar"):
            s = system(main_runs[t], sname)
            an = s["analysis"]
            row = {"tag": t, "label": plain(t), "system": sname, "nu": s["nu"],
                   "chi2_pdg": an["chi2_exp_only"],
                   "floor": 1e3 * s["sigma_th_GeV"],
                   "ci68_lo": 1e3 * an["sigma_th_ci68_GeV"][0],
                   "ci68_hi": 1e3 * an["sigma_th_ci68_GeV"][1],
                   "ci95_lo": 1e3 * an["sigma_th_ci95_GeV"][0],
                   "ci95_hi": 1e3 * an["sigma_th_ci95_GeV"][1]}
            z = scored(held_out_points(main_runs[t], sname, "loo"))
            row["loo_rms_t"] = rms([q["t"] for q in z]) if z else math.nan
            row["loo_rms_t_expected"] = expected_rms_t([q["nu_left"] for q in z])
            f0 = floor0_runs.get(t)
            s0 = system(f0, sname) if f0 else None
            if s0 and "chi2" in s0:
                row["chi2_pdg_refit"] = s0["chi2"]
                a0 = s0["analysis"]
                shifts = {n: abs(a0["value"][a0["names"].index(n)] - v) / e
                          for n, v, e in zip(an["names"], an["value"], an["sigma"])
                          if n in a0["names"] and e and e > 0}
                worst = max(shifts, key=shifts.get) if shifts else ""
                row["shift_par"] = worst
                row["shift_sigma"] = shifts.get(worst, math.nan)
            out.append(row)
    return out


def calibration_rows(main_runs, cv_runs, tags):
    """Out-of-sample calibration, pooled over `tags`, per system and method:
    the standardised misses t = miss / sd (every term of sd from the refit
    without the held-out data) against the coverage expected for a scale
    estimated with nu_left degrees of freedom (t distribution). Points of
    different models share the states, so the pooled counts are descriptive."""
    out = []
    for sname in ("ccbar", "bbar"):
        for kind in ("loo", "lmo"):
            pts = []
            for t in tags:
                run = main_runs.get(t) if kind == "loo" else cv_runs.get(t)
                if run:
                    pts += scored(held_out_points(run, sname, kind))
            if not pts:
                continue
            tt = np.array([q["t"] for q in pts])
            out.append({
                "system": sname, "kind": kind, "n": len(pts),
                "rms_t": rms(tt),
                "rms_t_expected": expected_rms_t([q["nu_left"] for q in pts]),
                "in1": float(np.mean(np.abs(tt) < 1)),
                "in1_expected": float(np.mean([ms.t_coverage(1, q["nu_left"]) for q in pts])),
                "in2": float(np.mean(np.abs(tt) < 2)),
                "in2_expected": float(np.mean([ms.t_coverage(2, q["nu_left"]) for q in pts])),
                "median_rms_t_per_model": float(np.median([
                    rms([q["t"] for q in scored(held_out_points(
                        (main_runs if kind == "loo" else cv_runs)[t], sname, kind))])
                    for t in tags if (main_runs if kind == "loo" else cv_runs).get(t)]))})
    return out


SCALING_LAWS = {"uniform": None, "v2": "v2", "p4": "p4_term_GeV"}


def _solve_scale(d, se, u, target):
    """kappa with sum d^2 / (se^2 + kappa^2 u^2) = target (diagonal errors)."""
    f = lambda k: float(np.sum(d ** 2 / (se ** 2 + (k * u) ** 2))) - target  # noqa: E731
    if f(0.0) <= 0:
        return 0.0
    hi = 1.0
    while f(hi) > 0:
        hi *= 2
    return float(optimize.brentq(f, 0.0, hi, xtol=1e-14 * hi))


def scaling_rows(main_runs, tags):
    """Is the model error one universal number across flavours? For each law
    the floor of state i is kappa u_i (u = 1, <v^2>_i, or <p^4>_i/8m^3 of the
    model's own fit). kappa is solved per system (chi2 = nu) and pooled (sum
    chi2 = nu_c + nu_b); the variance-ratio statistic F = (chi2_c/nu_c) /
    (chi2_b/nu_b) at the pooled kappa is ~ F(nu_c, nu_b) under one kappa.
    Parameters held at the uniform-floor fit; experimental errors diagonal."""
    out = []
    for t in tags:
        row = {"tag": t, "label": plain(t)}
        sysd = {}
        for sname in ("ccbar", "bbar"):
            s = system(main_runs[t], sname)
            fit = [x for x in s["states"] if x["in_fit"]]
            sysd[sname] = (np.array([x["model_GeV"] - x["pdg_GeV"] for x in fit]),
                           np.array([x["sigma_exp_GeV"] for x in fit]),
                           {law: np.array([1.0 if key is None else x[key] for x in fit])
                            for law, key in SCALING_LAWS.items()}, s["nu"])
        for law in SCALING_LAWS:
            k = {sn: _solve_scale(d, se, u[law], nu) for sn, (d, se, u, nu) in sysd.items()}
            nu_c, nu_b = sysd["ccbar"][3], sysd["bbar"][3]

            def chi2(sn, kap):
                d, se, u, _ = sysd[sn]
                return float(np.sum(d ** 2 / (se ** 2 + (kap * u[law]) ** 2)))
            f = lambda kap: chi2("ccbar", kap) + chi2("bbar", kap) - nu_c - nu_b  # noqa: E731
            kp = 0.0
            if f(0.0) > 0:
                hi = 1.0
                while f(hi) > 0:
                    hi *= 2
                kp = optimize.brentq(f, 0.0, hi, xtol=1e-14 * hi)
            F = (chi2("ccbar", kp) / nu_c) / (chi2("bbar", kp) / nu_b)
            row[f"{law}_kc"], row[f"{law}_kb"] = k["ccbar"], k["bbar"]
            row[f"{law}_ratio"] = k["ccbar"] / k["bbar"] if k["bbar"] > 0 else math.nan
            row[f"{law}_p"] = ms.variance_ratio_test(F, nu_c, nu_b)
        out.append(row)
    return out


# ------------------------------------------------ limits of the potential model --
THRESHOLD_GEV = {"ccbar": 3.7297, "bbar": 10.5590}  # D0 D0bar, B+ B-


def limits_rows(runs, tags):
    """Every state of each model: size, velocity, neglected kinetic term,
    distance to the open-flavour threshold, and the residual model - PDG.
    `resid_oos_MeV` is out of sample for every state: the exact leave-one-out
    miss for a fitted state, the prediction for the others. A fitted state
    without its LOO miss is left out."""
    out = []
    for t in tags:
        r = runs.get(t)
        if not r:
            continue
        for s in r["systems"]:
            sname = s["system"]
            if sname not in THRESHOLD_GEV or "sigma_th_GeV" not in s:
                continue
            pts = {q["label"]: q for q in held_out_points(r, sname, "loo")}
            miss = {k: q["miss"] for k, q in pts.items()}
            for x in s["states"]:
                resid = 1e3 * (x["model_GeV"] - x["pdg_GeV"])
                if x["in_fit"]:
                    if not math.isfinite(miss.get(x["label"], math.nan)):
                        continue
                    oos = 1e3 * miss[x["label"]]
                    err = 1e3 * pts[x["label"]]["sd"]
                else:
                    oos = resid
                    err = 1e3 * math.sqrt((s["sigma_th_GeV"] * x.get("floor_w", 1.0)) ** 2
                                          + (x.get("sigma_par_GeV") or 0) ** 2
                                          + (x.get("sigma_num_GeV") or 0) ** 2
                                          + (x.get("sigma_exp_GeV") or 0) ** 2)
                out.append({
                    "tag": t, "label": plain(t), "system": sname,
                    "state": x["label"], "in_fit": x["in_fit"],
                    "r_rms_fm": x["r_rms_fm"], "v2": x["v2"],
                    "p4_term_MeV": 1e3 * x["p4_term_GeV"],
                    "below_threshold_MeV": 1e3 * (THRESHOLD_GEV[sname] - x["pdg_GeV"]),
                    "resid_MeV": resid, "resid_oos_MeV": oos, "err_MeV": err,
                    "floor_MeV": 1e3 * s["sigma_th_GeV"] * x.get("floor_w", 1.0)})
    return out


def binned_residuals(rows, key, edges):
    """Out-of-sample residuals in bins of `key`, per system, split at the
    open-flavour threshold (size and threshold distance are nearly collinear).
    Each model is one sample: the RMS is taken per model, then the median over
    models; `states` counts distinct states, not model x state pairs."""
    out = []
    for sname in ("ccbar", "bbar"):
        for lo, hi in zip(edges, edges[1:]):
            b = [r for r in rows if r["system"] == sname and lo <= r[key] < hi]
            if not b:
                continue
            row = {"system": sname, "bin": f"{lo:g}–{hi:g}",
                   "states": len({r["state"] for r in b}),
                   "fitted": len({r["state"] for r in b if r["in_fit"]})}
            for side, keep in (("below", lambda r: r["below_threshold_MeV"] > 0),
                               ("above", lambda r: r["below_threshold_MeV"] <= 0)):
                per_model = {}
                for r in b:
                    if keep(r):
                        per_model.setdefault(r["tag"], []).append(r["resid_oos_MeV"])
                row[f"n_{side}"] = len({r["state"] for r in b if keep(r)})
                row[f"rms_{side}_MeV"] = (float(np.median([rms(v) for v in per_model.values()]))
                                          if per_model else math.nan)
            per_model = {}
            for r in b:
                per_model.setdefault(r["tag"], []).append(r["floor_MeV"])
            row["rms_floor_MeV"] = float(np.median([rms(v) for v in per_model.values()]))
            out.append(row)
    return out


# ------------------------------------------------ state-dependent floor study --
def floor_study_rows(base_runs, base_cv, shape_runs, tags, shapes):
    """Per model and floor shape: the rms floor per system, the LOO and LMO
    predictive density of the fitted states, its paired difference from the
    uniform floor with its standard error, and the rms pull of the states
    above threshold (reported only: not a selection criterion)."""
    out = []
    for t in tags:
        if t not in base_runs:
            continue
        ref = logdens_map(base_runs[t], base_cv.get(t))
        variants = [("uniform", base_runs[t], base_cv.get(t))] + [
            (sh, shape_runs.get(f"{t}-{sh}"), shape_runs.get(f"{t}-{sh}"))
            for sh in shapes]
        for sh, run, cvrun in variants:
            if not run:
                continue
            ld = logdens_map(run, cvrun)
            row = {"tag": t, "label": plain(t), "shape": sh}
            for sname in ("ccbar", "bbar"):
                row[f"{sname}_floor_MeV"] = 1e3 * system(run, sname)["sigma_th_GeV"]
            for kind in ("loo", "lmo"):
                row[f"elpd_{kind}"] = sum(v for k, v in ld.items() if k[0] == kind)
                row[f"d_{kind}"], row[f"se_{kind}"], _ = paired_elpd(ref, ld, kind)
            row["above_pull_rms"] = rms([x["pull"] for sname in ("ccbar", "bbar")
                                         for x in system(run, sname)["states"]
                                         if not x["in_fit"] and x["pull"] is not None])
            out.append(row)
    return out


def floor_study_summary(frows):
    """Per shape: totals, and how many models gain or lose by more than two
    standard errors of the paired difference from the uniform floor."""
    out = {}
    for r in frows:
        t = out.setdefault(r["shape"], {"shape": r["shape"], "models": 0, "elpd_loo": 0.0,
                                        "elpd_lmo": 0.0, "loo_up": 0, "loo_down": 0,
                                        "lmo_up": 0, "lmo_down": 0, "pulls": []})
        t["models"] += 1
        for kind in ("loo", "lmo"):
            t[f"elpd_{kind}"] += r[f"elpd_{kind}"]
            d, se = r[f"d_{kind}"], r[f"se_{kind}"]
            if math.isfinite(d) and math.isfinite(se):
                t[f"{kind}_up"] += d > 2 * se
                t[f"{kind}_down"] += d < -2 * se
        t["pulls"].append(r["above_pull_rms"])
    for t in out.values():
        t["above_pull_rms_median"] = float(np.median(t.pop("pulls")))
    return list(out.values())


# --------------------------------------------- semi-relativistic kinetic study --
def rel_rows(base_runs, base_cv, rel_runs):
    """Non-relativistic vs semi-relativistic kinetic energy, same model: the
    floor, the paired LOO + LMO predictive-density difference (rel - NR) with
    its standard error, the RMS miss above threshold and the B_c."""
    out = []
    for key, rr in sorted(rel_runs.items()):
        t = key[: -len("-rel")]
        nr = base_runs.get(t)
        if not nr:
            continue
        row = {"tag": t, "label": plain(t)}
        a, b = logdens_map(nr, base_cv.get(t)), logdens_map(rr, rr)
        for sname in ("ccbar", "bbar"):
            for kind, run in (("nr", nr), ("rel", rr)):
                s = system(run, sname)
                row[f"{kind}_{sname}_floor"] = 1e3 * s["sigma_th_GeV"]
                row[f"{kind}_{sname}_pred"] = rms([
                    1e3 * (x["model_GeV"] - x["pdg_GeV"]) for x in s["states"]
                    if not x["in_fit"]])
            row[f"d_{sname}"], row[f"se_{sname}"], _ = paired_elpd_by_state(a, b, sname)
            nj = system(rr, sname)["analysis"].get("nested")
            row[f"rel_{sname}_n"] = nj["chosen"] if nj else 0
        for kind, run in (("nr", nr), ("rel", rr)):
            bc = system(run, "bc")
            m = bc["states"][0]["model_GeV"] if bc else None
            row[f"{kind}_bc"] = (1e3 * (m - bc["states"][0]["pdg_GeV"])
                                 if m is not None else math.nan)
        out.append(row)
    return out


# ------------------------------------------------ truncation test (a-priori floor) --
def truncation_rows(main_runs, tags):
    """A/B with the non-relativistic kinetic energy: the a-priori floor sigma_a
    (the part of the omitted -<p^4>/8m^3 term the parameters cannot absorb),
    the fitted floor, chi2/nu at sigma_a with its p-value (the conventional
    goodness-of-fit test, with a floor not tuned to the data), and the
    pattern coefficient lambda in (PDG - model) = lambda R."""
    out = []
    for t in tags:
        for sname in ("ccbar", "bbar"):
            s = system(main_runs[t], sname)
            tr = s["analysis"].get("truncation") if s else None
            if not tr:
                continue
            out.append({"tag": t, "label": plain(t), "system": sname,
                        "sigma_a": 1e3 * tr["sigma_a_GeV"],
                        "floor": 1e3 * s["sigma_th_GeV"],
                        "ratio": tr["sigma_a_GeV"] / s["sigma_th_GeV"],
                        "nu": tr["nu"], "chi2_nu": tr["chi2"] / tr["nu"],
                        "p": tr["p"], "lam": tr["lambda"], "lam_err": tr["lambda_err"]})
    return out


# ------------------------------------------------ joint fit: universal b --
def joint_rows(main_runs, joint_runs, tags):
    """Separate and joint string tensions, the universality test, the floors
    and the B_c from the separate and the joint parameters."""
    out = []
    for t in tags:
        j = joint_runs.get(t)
        if not j or "joint" not in j:
            continue
        jj = j["joint"]
        row = {"tag": t, "label": plain(t), "shared": ",".join(jj["shared"]),
               "dchi2": jj["dchi2"], "k": jj["k"], "p": jj["p"],
               "floor_c": 1e3 * jj["sigma_th_ccbar_GeV"],
               "floor_b": 1e3 * jj["sigma_th_bbar_GeV"]}
        for sname, key in (("ccbar", "c"), ("bbar", "b")):
            a = system(main_runs[t], sname)["analysis"]
            if "b" in a["names"]:
                i = a["names"].index("b")
                row[f"b_{key}"], row[f"b_{key}_err"] = a["value"][i], a["sigma"][i]
            row[f"floor_{key}_sep"] = 1e3 * system(main_runs[t], sname)["sigma_th_GeV"]
        if "b" in jj["names"]:
            i = jj["names"].index("b")
            row["b_joint"], row["b_joint_err"] = jj["value"][i], jj["sigma"][i]
        for kind, run in (("sep", main_runs[t]), ("joint", j)):
            bc = system(run, "bc")
            m = bc["states"][0]["model_GeV"] if bc else None
            row[f"bc_{kind}"] = (1e3 * (m - bc["states"][0]["pdg_GeV"])
                                 if m is not None else math.nan)
        out.append(row)
    return out


# ------------------------------------------------ annihilation widths --
def width_rows(main_runs, tags):
    """Every width with a PDG value: LO and NLO (with its scale band) over
    the PDG, per model."""
    out = []
    for t in tags:
        for sname in ("ccbar", "bbar"):
            for w in system(main_runs[t], sname).get("widths", []):
                if not w.get("pdg_keV"):
                    continue
                out.append({"tag": t, "label": plain(t), "system": sname,
                            "state": w["label"], "process": w["process"],
                            "lo": w["lo_keV"], "nlo": w["nlo_keV"],
                            "nlo_lo": w["nlo_band_lo_keV"], "nlo_hi": w["nlo_band_hi_keV"],
                            "pdg": w["pdg_keV"], "pdg_err": w["pdg_err_keV"],
                            "derived": w.get("pdg_derived", 0),
                            "nlo_err_par": w.get("nlo_err_par_keV", math.nan),
                            "nlo_err_th": w.get("nlo_err_th_keV", math.nan),
                            "pull": w.get("pull", math.nan),
                            "lo_ratio": w["lo_keV"] / w["pdg_keV"],
                            "nlo_ratio": w["nlo_keV"] / w["pdg_keV"]})
    return out


WIDTH_RATIOS = [  # (numerator, denominator, process): normalisation and the
    # QCD factor (nearly) cancel, so these test the wave functions
    ("psi(2S)", "J/psi(1S)", "e+ e-"), ("Upsilon(2S)", "Upsilon(1S)", "e+ e-"),
    ("Upsilon(3S)", "Upsilon(1S)", "e+ e-"), ("Upsilon(4S)", "Upsilon(1S)", "e+ e-"),
    ("chi_c2(1P)", "chi_c0(1P)", "gamma gamma")]


def width_ratio_rows(main_runs, tags):
    """Width ratios at LO (model) against the PDG, per model."""
    out = []
    for t in tags:
        ws = {}
        for sname in ("ccbar", "bbar"):
            for w in system(main_runs[t], sname).get("widths", []):
                ws[(w["label"], w["process"])] = w
        for num, den, proc in WIDTH_RATIOS:
            a, b = ws.get((num, proc)), ws.get((den, proc))
            if not a or not b or not a.get("pdg_keV") or not b.get("pdg_keV"):
                continue
            pr = a["pdg_keV"] / b["pdg_keV"]
            out.append({"tag": t, "label": plain(t), "ratio": f"{num} / {den}",
                        "model_lo": a["lo_keV"] / b["lo_keV"],
                        "model_nlo": a["nlo_keV"] / b["nlo_keV"], "pdg": pr,
                        "pdg_err": pr * math.hypot(a["pdg_err_keV"] / a["pdg_keV"],
                                                   b["pdg_err_keV"] / b["pdg_keV"])})
    return out


def common_floor_rows(floor2_runs, tags):
    """chi2/nu and its p-value per system and summed, every model at the same
    a-priori floor (the floor2 sweep): the conventional goodness-of-fit test
    with one error budget, fixed before the comparison."""
    out = []
    for t in tags:
        r = floor2_runs.get(t)
        if not r:
            continue
        row = {"tag": t, "label": plain(t)}
        tot, nut = 0.0, 0
        for sname in ("ccbar", "bbar"):
            s = system(r, sname)
            if s is None or "chi2" not in s:
                break
            row[f"{sname}_chi2nu"] = s["chi2"] / s["nu"]
            row[f"{sname}_p"] = float(stats.chi2.sf(s["chi2"], s["nu"]))
            tot += s["chi2"]
            nut += s["nu"]
        else:
            row["chi2nu"] = tot / nut
            row["p"] = float(stats.chi2.sf(tot, nut))
            out.append(row)
    out.sort(key=lambda r: -r["p"])
    return out


def decay_rows(main_runs, stats_runs, tags):
    """The decay benchmark per model and system: chi2/N, p, the Bayesian
    posterior predictive p (stats sweep) and the largest pull."""
    out = []
    for t in tags:
        for sname in ("ccbar", "bbar"):
            db = system(main_runs[t], sname)["analysis"].get("decay_benchmark")
            if not db:
                continue
            i = int(np.argmax(np.abs(db["pulls"])))
            st = stats_runs.get(t)
            sdb = system(st, sname)["analysis"].get("decay_benchmark", {}) if st else {}
            out.append({"tag": t, "label": plain(t), "system": sname, "n": db["n"],
                        "chi2n": db["chi2"] / db["n"], "p": db["p"],
                        "ppp": sdb.get("ppp", math.nan), "worst": db["labels"][i],
                        "worst_pull": db["pulls"][i]})
    return out


def bayes_rows(stats_runs, tags):
    """The Bayesian posterior per model and system: sigma_th (median, 68 %),
    posterior predictive p-values, convergence."""
    out = []
    for t in tags:
        r = stats_runs.get(t)
        if not r:
            continue
        for sname in ("ccbar", "bbar"):
            m = system(r, sname)["analysis"].get("mcmc")
            if not m or "sigma_th_median_GeV" not in m:
                continue
            out.append({"tag": t, "label": plain(t), "system": sname,
                        "sth_med": 1e3 * m["sigma_th_median_GeV"],
                        "sth_lo": 1e3 * m["sigma_th_q16_GeV"],
                        "sth_hi": 1e3 * m["sigma_th_q84_GeV"],
                        "ppp_chi2": m["ppp_chi2"], "ppp_max": m["ppp_max"],
                        "rhat_max": max(m["rhat"]), "neff_min": min(m["n_eff"])})
    return out


def mass_width_rows(main_runs, width_runs, tags):
    """Masses and widths fitted together: the mass floor, the width chi2, and
    the largest parameter shift from the mass-only fit (in its errors)."""
    out = []
    for t in tags:
        r = width_runs.get(t)
        if not r:
            continue
        for sname in ("ccbar", "bbar"):
            m = system(r, sname)["analysis"].get("mass_width_fit")
            if not m:
                continue
            i = int(np.argmax(np.abs(m["shift_sigma"])))
            out.append({"tag": t, "label": plain(t), "system": sname,
                        "floor_mass_only": 1e3 * system(main_runs[t], sname)["sigma_th_GeV"],
                        "floor_combined": 1e3 * m["sigma_th_GeV"],
                        "chi2w_n": m["chi2_widths"] / m["n_widths"], "p_w": m["p_widths"],
                        "shift_par": m["names"][i], "shift": m["shift_sigma"][i]})
    return out
