#!/usr/bin/env python3
"""Rerun the Markov chain of a `stats` run after a change of the width code, and
put the new widths and decay benchmark into results/stats/<tag>.json.

The `stats` sweep runs a bootstrap, leave-one-out refits and a Markov-chain Monte
Carlo. Only the decay benchmark depends on the width code: its posterior
predictive p-value is computed from the widths at every draw. So only the chain is
rerun (same configuration, seed and number of steps; no bootstrap, no LOO) into
results/stats_mcmc/<tag>.json. The chain does not share random numbers with the
bootstrap, so the rerun gives back the stored chain; this is checked to 1e-6
(chi^2, parameters, R-hat, n_eff, ppp of chi^2 and the posterior summary)
before anything is written. Then the `widths` and `decay_benchmark` blocks of every
system are replaced, and the bootstrap, jackknife and LOO results are kept.

    python3 analysis/rerun_mcmc.py A-const-coff-smeared C-const-coff-s ...
    python3 analysis/rerun_mcmc.py --replace-chain D-free-coff

A chain that does not converge (R-hat far above 1) is chaotic: rounding at the
1e-16 level changes which steps are accepted, and the rerun is a different chain.
--replace-chain then takes the whole `mcmc` block from the rerun as well, and
marks it with "rerun": true. Fit, bootstrap and LOO are still checked to be the
stored ones (chi^2 and parameters).

Used in v1.1.0 for the 14 `stats` runs with S-D mixing in H (S-D sign fix): the
12 of the physical set were reproduced; the contrast models D-gi-coff and
D-free-coff, whose chains do not converge (R-hat up to 2.5, and ~1e12 in c c-bar
for D-free-coff), needed --replace-chain.
"""
import json
import subprocess
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import sweep as S  # noqa: E402


def mcmc_args(tag, threads):
    """The stats-sweep command line of `tag` without --bootstrap and --loo."""
    args = dict(S.configurations())[tag] + S.stats_args(tag) + S.seeds(tag, ["main"])
    out, i = [], 0
    while i < len(args):
        if args[i] == "--bootstrap":
            i += 2
        elif args[i] == "--loo":
            i += 1
        elif args[i] == "--threads":
            out += ["--threads", str(threads)]
            i += 2
        else:
            out.append(args[i])
            i += 1
    return out


def close(a, b, tol):
    if isinstance(a, dict):
        return a.keys() == b.keys() and all(close(a[k], b[k], tol) for k in a)
    if isinstance(a, list):
        return len(a) == len(b) and all(close(x, y, tol) for x, y in zip(a, b))
    if isinstance(a, float) or isinstance(b, float):
        return abs(a - b) <= tol * max(1.0, abs(a))
    return a == b


def same_chain(so, sn, tol=1e-6):
    """The same fit and chain to 1e-6: the refit may stop a few 1e-7 away from the
    stored minimum (the fit tolerance), which moves the chain at that level only; a
    chain with other random numbers differs in R-hat and n_eff at the 1e-2 level."""
    mo, mn = so["analysis"]["mcmc"], sn["analysis"]["mcmc"]
    return (abs(so["chi2"] - sn["chi2"]) < 1e-6 and close(so["theta"], sn["theta"], tol)
            and close(mo["rhat"], mn["rhat"], tol) and close(mo["n_eff"], mn["n_eff"], tol)
            and close(mo["ppp_chi2"], mn["ppp_chi2"], tol)
            and close(mo["summary"], mn["summary"], tol))


def rerun(tag, threads=2, replace_chain=False):
    d = S.RES / "stats_mcmc"
    d.mkdir(exist_ok=True)
    js = d / f"{tag}.json"
    if not js.exists():
        t0 = time.time()
        p = subprocess.run([str(S.BUILT), *mcmc_args(tag, threads), "--out", str(js)],
                           capture_output=True, text=True)
        if p.returncode:
            return f"FAILED ({p.returncode}) after {time.time() - t0:.0f}s"
    new = json.loads(js.read_text())
    path = S.RES / "stats" / f"{tag}.json"
    old = json.loads(path.read_text())
    pairs = [(so, sn) for so, sn in zip(old["systems"], new["systems"]) if "analysis" in so]
    same_fit = all(so["system"] == sn["system"] and abs(so["chi2"] - sn["chi2"]) < 1e-6
                   and close(so["theta"], sn["theta"], 1e-6) for so, sn in pairs)
    if not same_fit or not (replace_chain or all(same_chain(so, sn) for so, sn in pairs)):
        return "MISMATCH (not written)"
    msg = []
    for so, sn in pairs:
        a, b = so["analysis"].get("decay_benchmark", {}), sn["analysis"].get("decay_benchmark", {})
        msg.append(f"{so['system']}: ppp {a.get('ppp', float('nan')):.3f} -> "
                   f"{b.get('ppp', float('nan')):.3f}")
        so["widths"] = sn["widths"]
        so["analysis"]["decay_benchmark"] = b
        if replace_chain:
            so["analysis"]["mcmc"] = dict(sn["analysis"]["mcmc"], rerun=True)
    path.write_text(json.dumps(old, indent=2))
    return "ok  " + "; ".join(msg)


if __name__ == "__main__":
    replace = "--replace-chain" in sys.argv
    for t in [a for a in sys.argv[1:] if not a.startswith("--")]:
        print(f"{t:26s} {rerun(t, replace_chain=replace)}", flush=True)
