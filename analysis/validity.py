"""Physical and mathematical validity checks of the fits, from the sweep JSONs
(and the PDG file for the heavy-light mesons): spin structure against the PDG,
the quark-mass difference against heavy-quark symmetry, the string tension
against the lattice, the S-D mixing angle of psi(3770), and the wave
functions reconstructed from their logged parameters alone."""
import math
import sqlite3

import numpy as np
from scipy import optimize, special

from common import REPO, plain, system

PDG_FILE = REPO / "data" / "pdg-2026.0.sqlite"
# lattice string tension: sqrt(sigma) ~ 0.44 GeV (Bali, Phys. Rept. 343, 1 (2001))
SQRT_SIGMA_LATTICE = 0.44


# ------------------------------------------------ spin structure --
SPIN = {  # (system, label): what it tests
    "ccbar": {"hf1": ("J/psi(1S)", "eta_c(1S)"), "hf2": ("psi(2S)", "eta_c(2S)"),
              "P": ("chi_c0(1P)", "chi_c1(1P)", "chi_c2(1P)", "h_c(1P)")},
    "bbar": {"hf1": ("Upsilon(1S)", "eta_b(1S)"), "hf2": ("Upsilon(2S)", "eta_b(2S)"),
             "P": ("chi_b0(1P)", "chi_b1(1P)", "chi_b2(1P)", "h_b(1P)")},
}


def _spin_values(states, key, cov=None):
    """Splittings (MeV) and the fine-structure ratio; with the prediction
    covariance (GeV^2), their parametric errors (linear propagation)."""
    lab = [x["label"] for x in states]
    m = {x["label"]: 1e3 * x[key] for x in states if x[key] is not None}
    out = {}
    for sname, d in SPIN.items():
        if d["hf1"][0] not in m:
            continue
        lin = {"hf1": {d["hf1"][0]: 1, d["hf1"][1]: -1},
               "hf2": {d["hf2"][0]: 1, d["hf2"][1]: -1},
               "cog_h": {d["P"][0]: 1 / 9, d["P"][1]: 3 / 9, d["P"][2]: 5 / 9, d["P"][3]: -1}}
        for k, a in lin.items():
            out[k] = sum(c * m[n] for n, c in a.items())
        m0, m1, m2 = (m[k] for k in d["P"][:3])
        out["rho"] = (m2 - m1) / (m1 - m0)
        if cov is not None:
            C = 1e6 * np.asarray(cov)
            for k, a in lin.items():
                v = np.zeros(len(lab))
                for n, c in a.items():
                    v[lab.index(n)] = c
                out[f"{k}_err"] = float(np.sqrt(max(v @ C @ v, 0.0)))
            # rho = (m2 - m1)/(m1 - m0): gradient in (m0, m1, m2)
            g = np.zeros(len(lab))
            den = m1 - m0
            g[lab.index(d["P"][2])] = 1 / den
            g[lab.index(d["P"][1])] = -1 / den - (m2 - m1) / den ** 2
            g[lab.index(d["P"][0])] = (m2 - m1) / den ** 2
            out["rho_err"] = float(np.sqrt(max(g @ C @ g, 0.0)))
    return out


def spin_structure_rows(main_runs, tags):
    """1S and 2S hyperfine splittings, the 1P fine-structure ratio (M2 -
    M1)/(M1 - M0), and the 1P centre of gravity minus h (zero for a contact
    spin-spin force at first order), model and PDG, with the spin mode."""
    out = []
    for t in tags:
        r = main_runs[t]
        for sname in ("ccbar", "bbar"):
            s = system(r, sname)
            st = s["states"]
            mod = _spin_values(st, "model_GeV", s.get("pred_cov_GeV2"))
            pdg = _spin_values(st, "pdg_GeV")
            if not mod:
                continue
            out.append({"tag": t, "label": plain(t), "system": sname,
                        # C/D: the relativized spin terms are in H
                        "spin": ("GI, in H" if r["config"]["model"] in ("C", "D")
                                 else r["config"]["spin"]),
                        **{f"{k}_model" if not k.endswith("_err") else
                           f"{k[:-4]}_model_err": v for k, v in mod.items()},
                        **{f"{k}_pdg": v for k, v in pdg.items()}})
    return out


# ------------------------------------------------ heavy-quark masses --
def pdg_mass(mcid):
    """Summary-table mass (GeV) of a PDG Monte Carlo id."""
    with sqlite3.connect(PDG_FILE) as db:
        v, unit = db.execute(
            "SELECT d.value, d.unit_text FROM pdgparticle p JOIN pdgid i ON "
            "i.parent_pdgid = p.pdgid JOIN pdgdata d ON d.pdgid_id = i.id WHERE "
            "p.mcid = ? AND i.data_type = 'M' AND d.in_summary_table = 1",
            (mcid,)).fetchone()
    return v * {"MeV": 1e-3, "GeV": 1.0}[unit]


def quark_mass_window():
    """m_b - m_c from heavy-quark symmetry: the spin-averaged B and D masses,
    (M + 3 M*)/4, plus the kinetic correction -lambda1/2 (1/m_c - 1/m_b) for
    lambda1 in [-0.5, 0] GeV^2 (m_c = 1.5, m_b = 4.8 GeV): a plausibility
    window, not a measurement."""
    d = (pdg_mass(421) + 3 * pdg_mass(423)) / 4
    b = (pdg_mass(511) + 3 * pdg_mass(513)) / 4
    corr = 0.5 * 0.5 * (1 / 1.5 - 1 / 4.8)
    return b - d, b - d + corr


def quark_mass_rows(main_runs, joint_runs, tags):
    lo, hi = quark_mass_window()
    out = []
    for t in tags:
        sc, sb = system(main_runs[t], "ccbar"), system(main_runs[t], "bbar")
        mc, mb = sc["theta"][0], sb["theta"][0]
        row = {"tag": t, "label": plain(t), "dm_sep": mb - mc, "lo": lo, "hi": hi,
               "dm_sep_err": math.hypot(sc["analysis"]["sigma"][0], sb["analysis"]["sigma"][0])}
        j = joint_runs.get(t, {}).get("joint")
        if j:
            v = dict(zip(j["names"], j["value"]))
            e = dict(zip(j["names"], j["sigma"]))
            row["dm_joint"] = v.get("m1(b)", math.nan) - v.get("m1", math.nan)
            # the joint fit correlates the masses; the quadrature sum bounds it
            row["dm_joint_err"] = math.hypot(e.get("m1(b)", math.nan), e.get("m1", math.nan))
        out.append(row)
    return out


# ------------------------------------------------ wave functions from the JSON --
def basis_overlap(nu, norm, L):
    """S_ij of the normalised Gaussians N r^L exp(-nu r^2)."""
    s = nu[:, None] + nu[None, :]
    return (np.outer(norm, norm) * special.gamma(L + 1.5) / (2 * s ** (L + 1.5)))


def radial(nu, norm, c, L, r):
    """R_L(r) = sum_n c_n N_nL r^L exp(-nu_n r^2) on an array r (GeV^-1)."""
    return (r[:, None] ** L * np.exp(-np.outer(r * r, nu))) @ (np.asarray(c) * norm)


def wavefunction_check(run):
    """Reconstruct every logged wave function from the JSON alone (nu, N_nL,
    c): its norm, the radial value or derivative at the origin, and the
    radial nodes of the dominant component. Returns rows per state."""
    out = []
    for s in run["systems"]:
        if "basis" not in s:
            continue
        nu = np.array(s["basis"]["nu_GeV2"])
        norms = {int(k): np.array(v) for k, v in s["basis"]["norm_by_L"].items()}
        for x in s["states"]:
            w = x.get("wavefunction")
            if not w:
                continue
            tot, comps = 0.0, w["components"]
            for cpt in comps:
                c = np.array(cpt["c"])
                tot += c @ basis_overlap(nu, norms[cpt["L"]], cpt["L"]) @ c
            dom = max(comps, key=lambda q: q["weight"])
            # nodes within 4 r_rms where |R| > 1e-3 max|R|, as in qgem
            r = np.linspace(1e-4, 4 / math.sqrt(nu.min()), 20000)
            R = radial(nu, norms[dom["L"]], dom["c"], dom["L"], r)
            rr = math.sqrt(np.trapezoid(R ** 2 * r ** 4, r) / np.trapezoid(R ** 2 * r ** 2, r))
            r = np.linspace(rr * 4 / 4000, 4 * rr, 4000)
            R = radial(nu, norms[dom["L"]], dom["c"], dom["L"], r)
            keep = np.abs(R) > 1e-3 * np.abs(R).max()
            sg = np.sign(R[keep])
            inside = np.ones(int(keep.sum()), bool)
            nodes = int(np.sum(sg[inside][1:] != sg[inside][:-1]))
            out.append({"system": s["system"], "state": x["label"], "n": x["n"],
                        "norm_json": tot, "norm_logged": w["norm"],
                        "nodes_json": nodes, "nodes_logged": w["nodes"]})
    return out


def origin_values(run, sname):
    """{label: (R_S(0), R''_D(0), m_Q)} from the logged wave functions."""
    s = system(run, sname)
    nu = np.array(s["basis"]["nu_GeV2"])
    norms = {int(k): np.array(v) for k, v in s["basis"]["norm_by_L"].items()}
    out = {}
    for x in s["states"]:
        w = x.get("wavefunction")
        if not w:
            continue
        RS = sum(np.dot(q["c"], norms[0]) for q in w["components"] if q["L"] == 0)
        RD = sum(2 * np.dot(q["c"], norms[2]) for q in w["components"] if q["L"] == 2)
        out[x["label"]] = (RS, RD, s["theta"][0])
    return out


def mixing_angle(run):
    """The 2S-1D mixing angle that reproduces the PDG ratio Gamma_ee(psi(3770))
    / Gamma_ee(psi(2S)) with the model's own unmixed R_2S(0) and R''_1D(0):
    psi(2S) = cos t |2S> - sin t |1D>, psi(3770) = sin t |2S> + cos t |1D>;
    both share the QCD factor, so it cancels. Returns (theta_deg list, the
    model's own angle from its tensor mixing)."""
    ov = origin_values(run, "ccbar")
    s = system(run, "ccbar")
    w = {x["label"]: x for x in s["widths"]}
    st = {x["label"]: x for x in s["states"]}
    if "psi(2S)" not in ov or "psi(3770)" not in ov:
        return [], math.nan
    R2S, _, m = ov["psi(2S)"]
    _, R1D, _ = ov["psi(3770)"]
    D = 5 * R1D / (2 * math.sqrt(2) * m * m)
    M2, M3 = st["psi(2S)"]["pdg_GeV"], st["psi(3770)"]["pdg_GeV"]
    a, b = w["psi(3770)"], w["psi(2S)"]
    ratio = a["pdg_keV"] / b["pdg_keV"]
    ratio_err = ratio * math.hypot(a["pdg_err_keV"] / a["pdg_keV"],
                                   b["pdg_err_keV"] / b["pdg_keV"])

    def roots_for(target):
        def f(t):
            return (M2 / M3) ** 2 * (math.sin(t) * R2S + math.cos(t) * D) ** 2 \
                / (math.cos(t) * R2S - math.sin(t) * D) ** 2 - target
        grid = np.linspace(-math.pi / 3, math.pi / 3, 2001)
        vals = [f(t) for t in grid]
        return [optimize.brentq(f, x0, x1) for x0, x1, fa, fb in
                zip(grid, grid[1:], vals, vals[1:]) if fa * fb < 0]
    roots = roots_for(ratio)
    lo, hi = roots_for(ratio - ratio_err), roots_for(ratio + ratio_err)
    errs = []
    for k, t0 in enumerate(roots):  # +- from the PDG ratio error
        near = [min(x, key=lambda q: abs(q - t0)) for x in (lo, hi) if x]
        errs.append(max((abs(q - t0) for q in near), default=math.nan))
    own = 0.0
    for q in st["psi(3770)"]["wavefunction"]["components"]:
        if q["L"] == 0:
            own = math.degrees(math.asin(math.sqrt(q["weight"])))
    return [(math.degrees(t), math.degrees(e)) for t, e in zip(roots, errs)], own


def string_tension_rows(main_runs, joint_runs, tags):
    """b per system and joint against the lattice sqrt(sigma) (plausibility)."""
    out = []
    for t in tags:
        row = {"tag": t, "label": plain(t)}
        for sname in ("ccbar", "bbar"):
            a = system(main_runs[t], sname)["analysis"]
            if "b" in a["names"]:
                i = a["names"].index("b")
                row[f"sqrtb_{sname}"] = math.sqrt(a["value"][i])
                row[f"sqrtb_{sname}_err"] = a["sigma"][i] / (2 * math.sqrt(a["value"][i]))
        j = joint_runs.get(t, {}).get("joint")
        if j and "b" in j["names"]:
            i = j["names"].index("b")
            row["sqrtb_joint"] = math.sqrt(j["value"][i])
            row["sqrtb_joint_err"] = j["sigma"][i] / (2 * math.sqrt(j["value"][i]))
        out.append(row)
    return out


def summary(spin_rows, qrows, srows, wrows, mix, own):
    """Pass/fail lines of the audit (failures are kept)."""
    lines = []

    def med(rows, k):
        v = [r[k] for r in rows if k in r and math.isfinite(r[k])]
        return float(np.median(v)) if v else math.nan
    for sname in ("ccbar", "bbar"):
        rs = [r for r in spin_rows if r["system"] == sname]
        if not rs:
            continue
        for mode in sorted({r["spin"] for r in rs}):
            m = [r for r in rs if r["spin"] == mode]
            lines.append(
                f"{sname}, spin {mode} ({len(m)} models, medians ± median parametric error): "
                f"1S hyperfine {med(m, 'hf1_model'):.1f} ± {med(m, 'hf1_model_err'):.1f} MeV "
                f"(PDG {m[0]['hf1_pdg']:.1f}); fine-structure ratio {med(m, 'rho_model'):.2f} ± "
                f"{med(m, 'rho_model_err'):.2f} (PDG {m[0]['rho_pdg']:.2f}); cog − h "
                f"{med(m, 'cog_h_model'):.1f} ± {med(m, 'cog_h_model_err'):.1f} MeV "
                f"(PDG {m[0]['cog_h_pdg']:.2f})")
    if qrows:
        lo, hi = qrows[0]["lo"], qrows[0]["hi"]
        ins = sum(lo <= r["dm_sep"] <= hi for r in qrows)
        lines.append(f"m_b − m_c within the heavy-quark window [{lo:.3f}, {hi:.3f}] GeV: "
                     f"{ins} of {len(qrows)} separate fits (median {med(qrows, 'dm_sep'):.3f} ± "
                     f"{med(qrows, 'dm_sep_err'):.3f}); joint fits median "
                     f"{med(qrows, 'dm_joint'):.3f} ± {med(qrows, 'dm_joint_err'):.3f} GeV")
    if srows:
        lines.append(f"√b: median c c̄ {med(srows, 'sqrtb_ccbar'):.3f} ± "
                     f"{med(srows, 'sqrtb_ccbar_err'):.3f}, b b̄ {med(srows, 'sqrtb_bbar'):.3f} ± "
                     f"{med(srows, 'sqrtb_bbar_err'):.3f}, joint {med(srows, 'sqrtb_joint'):.3f} ± "
                     f"{med(srows, 'sqrtb_joint_err'):.3f} GeV (lattice √σ ≈ "
                     f"{SQRT_SIGMA_LATTICE} GeV)")
    if wrows:
        dn = max(abs(r["norm_json"] - 1) for r in wrows)
        bad = sum(r["nodes_json"] != r["n"] for r in wrows)
        lines.append(f"wave functions rebuilt from the JSON: max |norm − 1| = {dn:.1e}; "
                     f"radial nodes ≠ radial quantum number in {bad} of {len(wrows)} states")
    if mix is not None:
        lines.append(f"ψ(3770): mixing angle reproducing Γee(ψ(3770))/Γee(ψ(2S)) with the "
                     f"model's own R_2S(0), R″_1D(0): "
                     f"{', '.join(f'{a:.1f}° ± {e:.1f}°' for a, e in mix) or 'none'}; the tensor force "
                     f"alone gives {own:.1f}°")
    return lines
