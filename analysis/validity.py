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


def _spin_values(states, key):
    m = {x["label"]: 1e3 * x[key] for x in states if x[key] is not None}
    s = SPIN
    out = {}
    for sname, d in s.items():
        if d["hf1"][0] not in m:
            continue
        out["hf1"] = m[d["hf1"][0]] - m[d["hf1"][1]]
        out["hf2"] = m[d["hf2"][0]] - m[d["hf2"][1]]
        m0, m1, m2, h = (m[k] for k in d["P"])
        out["rho"] = (m2 - m1) / (m1 - m0)            # fine-structure ratio
        out["cog_h"] = (m0 + 3 * m1 + 5 * m2) / 9 - h  # P-wave hyperfine
    return out


def spin_structure_rows(main_runs, tags):
    """1S and 2S hyperfine splittings, the 1P fine-structure ratio (M2 -
    M1)/(M1 - M0), and the 1P centre of gravity minus h (zero for a contact
    spin-spin force at first order), model and PDG, with the spin mode."""
    out = []
    for t in tags:
        r = main_runs[t]
        for sname in ("ccbar", "bbar"):
            st = system(r, sname)["states"]
            mod, pdg = _spin_values(st, "model_GeV"), _spin_values(st, "pdg_GeV")
            if not mod:
                continue
            out.append({"tag": t, "label": plain(t), "system": sname,
                        "spin": r["config"]["spin"],
                        **{f"{k}_model": v for k, v in mod.items()},
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
        mc = system(main_runs[t], "ccbar")["theta"][0]
        mb = system(main_runs[t], "bbar")["theta"][0]
        row = {"tag": t, "label": plain(t), "dm_sep": mb - mc, "lo": lo, "hi": hi}
        j = joint_runs.get(t, {}).get("joint")
        if j:
            v = dict(zip(j["names"], j["value"]))
            row["dm_joint"] = v.get("m1(b)", math.nan) - v.get("m1", math.nan)
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
    target = w["psi(3770)"]["pdg_keV"] / w["psi(2S)"]["pdg_keV"]

    def f(t):
        return (M2 / M3) ** 2 * (math.sin(t) * R2S + math.cos(t) * D) ** 2 \
            / (math.cos(t) * R2S - math.sin(t) * D) ** 2 - target
    grid = np.linspace(-math.pi / 3, math.pi / 3, 2001)
    vals = [f(t) for t in grid]
    roots = [optimize.brentq(f, a, b) for a, b, fa, fb in
             zip(grid, grid[1:], vals, vals[1:]) if fa * fb < 0]
    own = 0.0
    for q in st["psi(3770)"]["wavefunction"]["components"]:
        if q["L"] == 0:
            own = math.degrees(math.asin(math.sqrt(q["weight"])))
    return [math.degrees(t) for t in roots], own


def string_tension_rows(main_runs, joint_runs, tags):
    """b per system and joint against the lattice sqrt(sigma) (plausibility)."""
    out = []
    for t in tags:
        row = {"tag": t, "label": plain(t)}
        for sname in ("ccbar", "bbar"):
            a = system(main_runs[t], sname)["analysis"]
            if "b" in a["names"]:
                row[f"sqrtb_{sname}"] = math.sqrt(a["value"][a["names"].index("b")])
        j = joint_runs.get(t, {}).get("joint")
        if j and "b" in j["names"]:
            row["sqrtb_joint"] = math.sqrt(j["value"][j["names"].index("b")])
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
                f"{sname}, spin {mode} ({len(m)} models): 1S hyperfine "
                f"{med(m, 'hf1_model'):.1f} MeV (PDG {m[0]['hf1_pdg']:.1f}); "
                f"fine-structure ratio {med(m, 'rho_model'):.2f} (PDG {m[0]['rho_pdg']:.2f}); "
                f"cog − h {med(m, 'cog_h_model'):.1f} MeV (PDG {m[0]['cog_h_pdg']:.2f})")
    if qrows:
        lo, hi = qrows[0]["lo"], qrows[0]["hi"]
        ins = sum(lo <= r["dm_sep"] <= hi for r in qrows)
        lines.append(f"m_b − m_c within the heavy-quark window [{lo:.3f}, {hi:.3f}] GeV: "
                     f"{ins} of {len(qrows)} separate fits (median {med(qrows, 'dm_sep'):.3f}); "
                     f"joint fits median {med(qrows, 'dm_joint'):.3f} GeV")
    if srows:
        lines.append(f"√b: median c c̄ {med(srows, 'sqrtb_ccbar'):.3f}, b b̄ "
                     f"{med(srows, 'sqrtb_bbar'):.3f}, joint {med(srows, 'sqrtb_joint'):.3f} GeV "
                     f"(lattice √σ ≈ {SQRT_SIGMA_LATTICE} GeV)")
    if wrows:
        dn = max(abs(r["norm_json"] - 1) for r in wrows)
        bad = sum(r["nodes_json"] != r["n"] for r in wrows)
        lines.append(f"wave functions rebuilt from the JSON: max |norm − 1| = {dn:.1e}; "
                     f"radial nodes ≠ radial quantum number in {bad} of {len(wrows)} states")
    if mix is not None:
        lines.append(f"ψ(3770): mixing angle reproducing Γee(ψ(3770))/Γee(ψ(2S)) with the "
                     f"model's own R_2S(0), R″_1D(0): "
                     f"{', '.join(f'{a:.1f}°' for a in mix) or 'none'}; the tensor force "
                     f"alone gives {own:.1f}°")
    return lines
