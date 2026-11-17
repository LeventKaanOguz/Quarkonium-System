#!/usr/bin/env python3
"""Tables of the v1.2 sweeps (rep, oneloop, cx, screen) from results/*/*.json.
Writes results/v1.2-summary.md and prints it. Usage: python3 analysis/summarize_v12.py"""
import json, math
from pathlib import Path
RES = Path(__file__).resolve().parent.parent / "results"
out = []
def P(s=""): out.append(s)
def load(sweep, tag): return json.load(open(RES / sweep / f"{tag}.json"))
def floors(r): return [1e3 * s["sigma_th_GeV"] for s in r["systems"][:2]]
def loo(s):
    m = s["analysis"].get("loo_exact_miss_GeV")
    return 1e3 * math.sqrt(sum(x * x for x in m) / len(m)) if m else float("nan")
def mu_of(s):
    a = s["analysis"]
    if "mu_s" in a["names"]:
        i = a["names"].index("mu_s"); return a["value"][i], a["sigma"][i]
    return s["theta"][s["theta_names"].index("mu_s")], None

P("# v1.2 results\n")
P("## rep: representatives against v1.1.0 (`main`)\n")
P("| model | max rel dtheta | max dM (MeV) | B_c miss | par | interp |\n|---|---|---|---|---|---|")
for t in ["A-const-coff-first", "B-apt-coff-first", "C-const-coff-s", "D-gi-coff-s"]:
    a, b = load("main", t)["systems"], load("rep", t)["systems"]
    dth = max(abs(x - y) / max(abs(x), 1e-12) for k in range(2) for x, y in zip(a[k]["theta"], b[k]["theta"]))
    dm = 1e3 * max(abs(x["model_GeV"] - y["model_GeV"]) for k in range(2) for x, y in zip(a[k]["states"], b[k]["states"]))
    bc = b[2]["states"][0]
    P(f"| {t} | {dth:.1e} | {dm:.1e} | {1e3*(bc['model_GeV']-bc['pdg_GeV']):+.1f} | {1e3*bc['sigma_par_GeV']:.1f} | {1e3*bc['sigma_interp_GeV']:.1f} |")

P("\n## oneloop\n")
P("| tag | floor cc | floor bb | LOO cc | LOO bb | chi_c0 miss | chi_c0 gg (keV) | Ups(1S) ee (keV) | B_c miss |\n|---|---|---|---|---|---|---|---|---|")
for t in ["A", "B"]:
    for suf in ["tree", "ol1", "ol2", "ol3"]:
        tag = f"{t}-const-coff-first-{suf}"; r = load("oneloop", tag); S = r["systems"]
        st = {x["label"]: x for x in S[0]["states"]}
        w = {x["label"] + "/" + x["process"]: x for x in S[0]["widths"]}
        wb = {x["label"] + "/" + x["process"]: x for x in S[1]["widths"]}
        c0 = st["chi_c0(1P)"]; bc = S[2]["states"][0]
        P(f"| {tag} | {floors(r)[0]:.2f} | {floors(r)[1]:.2f} | {loo(S[0]):.2f} | {loo(S[1]):.2f} | {1e3*(c0['model_GeV']-c0['pdg_GeV']):+.1f} | "
          f"{w['chi_c0(1P)/gamma gamma']['nlo_keV']:.2f} | {wb['Upsilon(1S)/e+ e-']['nlo_keV']:.2f} | {1e3*(bc['model_GeV']-bc['pdg_GeV']):+.0f} |")

P("\n## cx: mixed real + complex-range basis (25 + 6) against the real-range fit (37)\n")
P("| model | system | floor real | floor cx | max dM over all PDG states (MeV) |\n|---|---|---|---|---|")
for t in ["A-const-coff-first", "B-apt-coff-first", "A-apt-coff-smeared", "B-const-coff-smeared"]:
    a, b = load("main", t)["systems"], load("cx", t + "-cx")["systems"]
    for k, nm in enumerate(["ccbar", "bbar"]):
        dm = 1e3 * max(abs(x["model_GeV"] - y["model_GeV"]) for x, y in zip(a[k]["states"], b[k]["states"]))
        P(f"| {t} | {nm} | {1e3*a[k]['sigma_th_GeV']:.3f} | {1e3*b[k]['sigma_th_GeV']:.3f} | {dm:.4f} |")

P("\n## screen: the screening mass as an optimization parameter\n")
P("| model | run | mu_s cc (GeV) | mu_s bb (GeV) | floor cc | floor bb |\n|---|---|---|---|---|---|")
for t in ["B-apt-coff-first", "D-gi-coff-s"]:
    for suf in ["fit", "shared"] + [f"mu{m:g}" for m in (0.02, 0.05, 0.1, 0.2, 0.3, 0.5)]:
        r = load("screen", f"{t}-{suf}"); S = r["systems"]
        m = [mu_of(s) for s in S[:2]]
        f = lambda x: f"{x[0]:.4f}" + (f" ± {x[1]:.3f}" if x[1] else "")
        extra = ""
        if suf == "shared" and "joint" in r:
            jn = r["joint"]
            k = jn["names"].index("mu_s") if "mu_s" in jn["names"] else None
            val = f"{jn['value'][k]:.4f} ± {jn['sigma'][k]:.3f}" if k is not None else "n/a"
            extra = f" joint mu_s = {val}; delta chi2 = {jn['dchi2']:.2f} (P = {jn['p']:.2f})"
        P(f"| {t} | {suf} | {f(m[0])} | {f(m[1])} | {floors(r)[0]:.2f} | {floors(r)[1]:.2f} |{extra}")
text = "\n".join(out) + "\n"
(RES / "v1.2-summary.md").write_text(text)
print(text)
