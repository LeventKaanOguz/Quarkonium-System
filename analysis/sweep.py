#!/usr/bin/env python3
"""Run qgem over every model configuration, in parallel.

Each run writes results/<sweep>/<tag>.json (machine-readable) and <tag>.txt
(the printed table). A run whose JSON exists is skipped, so an interrupted
sweep resumes where it stopped.

    python3 analysis/sweep.py all       # every sweep below, in order
    python3 analysis/sweep.py stop      # stop a running sweep (and its fits)

    lm      Levenberg-Marquardt, 9 starts
    minuit  Minuit2 Migrad from the default start, Hesse and Minos
    main    the final fits: seeded with the minima of lm and minuit for the
            configuration and the models nested in it; exact leave-one-out
            with the floor re-solved in every refit
    floor0  the conventional fit, PDG errors only (sigma_th = 0), physical set
    floor   common fixed floor (10 MeV): sensitivity of the same
    floor2  per-system a-priori floors (the reference model's truncation
            estimate): chi2/nu as a test, AIC, AICc, BIC, evidence
    joint   c c-bar + b b-bar with one string tension b: universality, B_c
    joint2  the same with the flavour-independent coupling or smearing shared
    cv      leave-multiplet-out, repeated 4-fold, learning curves
    stats   bootstrap, MCMC, jackknife (the physical set and two contrasts)
    nested  basis convergence of every state on nested sizes 10 ... 145
    bopt    each model's own variational (r1, rmax)
    bgrid   the (r1, rmax) plane: level sum and fit
    shape   study: state-dependent floors (shapes from the reference model's
            <v^2>, <p^4>/8m^3, r_rms) on the physical set
    rel     study: semi-relativistic kinetic energy for A/B (running couplings)

lm and minuit come first, main reads them, and the rest read main. The qgem
binary is copied once per invocation (or taken from $QGEM_BIN), so rebuilding
during a sweep cannot affect it.
"""
import argparse
import itertools
import json
import math
import os
import re
import shutil
import signal
import subprocess
import sys
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
RES = ROOT / "results"
BUILT = ROOT / "build" / "release" / "qgem"
PIDDIR = RES

SYSTEMS = ("ccbar", "bbar")
# representative models: the best of each family and the notes' worked fit
REPRESENTATIVE = ["A-apt-coff-smeared", "A-frozen-coff-smeared",
                  "A-gi-coff-smeared", "A-const-coff-first", "B-gi-coff-smeared",
                  "C-gi-coff-s", "D-const-coff-s", "D-gi-coff-s"]
# degenerate contrasts for the error methods (the sigma0 parametrisation of GI)
CONTRASTS = ["D-gi-coff", "D-free-coff"]
BGRID_SET = ["A-apt-coff-smeared", "D-const-coff-s"]
# the best uniform-floor model: its per-state diagnostics define the floor
# shapes, frozen and shared by every model (one error model for all)
REFERENCE = "A-apt-coff-smeared"
SHAPES = ("v2", "p4", "r")


def configurations():
    """Every model the CLI allows: (tag, args)."""
    out = []
    for model, running, c, spin in itertools.product(
            "AB", ["const", "gi", "frozen", "apt"], ["off", "on"],
            ["none", "first", "smeared"]):
        out.append((f"{model}-{running}-c{c}-{spin}",
                    ["--model", model, "--running", running, "--c", c,
                     "--spin", spin]))
    for model, running, c in itertools.product("CD", ["const", "gi", "free"],
                                                ["off", "on"]):
        out.append((f"{model}-{running}-c{c}",
                    ["--model", model, "--running", running, "--c", c]))
    # screened models with the screening mass fitted
    for running, spin in itertools.product(["const", "gi", "frozen", "apt"],
                                           ["first", "smeared"]):
        out.append((f"B-{running}-coff-{spin}-fitmu",
                    ["--model", "B", "--running", running, "--spin", spin,
                     "--fit-mu-s"]))
    for running in ["const", "gi", "free"]:
        out.append((f"D-{running}-coff-fitmu",
                    ["--model", "D", "--running", running, "--fit-mu-s"]))
    # Modified GI: the smearing width enters only through
    # sigma12^2 = sigma0^2 (...) + s^2 (2 m1 m2 / (m1 + m2))^2, so within one
    # system sigma0 and s are degenerate. Fix sigma0 at its lower bound and fit
    # the slope s instead: identifiable parameters, and sigma12 may go below
    # 1.55 m, which the original parametrisation forbids.
    for tag, args in list(out):
        if tag[0] in "CD" and "-coff" in tag:
            out.append((tag + "-s", args + ["--fix", "sigma0=0.1", "--free", "s"]))
    return out


def physical_set():
    """c = 0, spin terms on, Fisher matrix not degenerate in either system
    (read from the main sweep)."""
    out = []
    for t, _ in configurations():
        js = RES / "main" / f"{t}.json"
        if not js.exists():
            continue
        r = json.loads(js.read_text())
        cfg = r["config"]
        conds = [(s.get("analysis", {}).get("condition") or math.inf)
                 for s in r["systems"] if s["system"] in ("ccbar", "bbar")]
        if (not cfg["c_free"] and cfg.get("spin") != "none" and conds
                and max(conds) < 1e10):
            out.append(t)
    return [t for t in out if not duplicate_fit(t, out)]


def duplicate_fit(tag, tags):
    """True if tag + '-s' is in `tags` and gives the same c c-bar and b b-bar
    fit. For equal masses sigma12^2 = sigma0^2 + s^2 m^2, so (sigma0 free,
    s = 1.55) and (sigma0 = 0.1, s free) coincide whenever the optimum has
    sigma12 > 1.55 m; the s form, which spans every sigma12 > 0.1 GeV, is kept.
    The B_c (unequal masses) still differs; it is not fitted."""
    if tag + "-s" not in tags:
        return False
    floors = []
    for t in (tag, tag + "-s"):
        r = json.loads((RES / "main" / f"{t}.json").read_text())
        floors.append([s.get("sigma_th_GeV", math.nan) for s in r["systems"]
                       if s["system"] in ("ccbar", "bbar")])
    return all(abs(a - b) <= 1e-5 * max(abs(a), abs(b)) for a, b in zip(*floors))


def nulls(tag):
    """The configurations nested in `tag` (one parameter or group fewer)."""
    if tag.endswith("-s"):                  # modified GI: same nesting
        return [n + "-s" for n in nulls(tag[:-2]) if "-coff" in n]
    out = []
    parts = tag.split("-")
    if "con" in parts:                      # additive constant c -> 0
        out.append("-".join("coff" if p == "con" else p for p in parts))
    if tag.endswith("-fitmu"):              # mu_s fixed, or mu_s -> 0 (A, C)
        base = tag[: -len("-fitmu")]
        out += [base, {"B": "A", "D": "C"}[base[0]] + base[1:]]
    if tag[0] in "CD" and "-free-" in tag:  # free GI amplitudes -> GI values
        out.append(tag.replace("-free-", "-gi-"))
    return out


def depth(tag):
    """0 for a configuration with no nested model, else 1 + the deepest null."""
    return 1 + max(map(depth, nulls(tag))) if nulls(tag) else 0


def seeds(tag, sources):
    """--start arguments: the minima stored in results/<source>/ for `tag` and
    for the models nested in it. A nested model's minimum is a point of the
    extended model, so starting there guarantees chi2(extended) <= chi2(null)."""
    args, seen = [], []
    family = [tag] + nulls(tag)
    if tag.endswith("-s"):                  # also the sigma0 parametrisation
        family += [tag[:-2]] + nulls(tag[:-2])
    for t in family:
        for source in sources:
            src, _, suffix = source.partition(":")   # "floor2:-ccbar"
            js = RES / src / f"{t}{suffix}.json"
            if not js.exists():
                continue
            for s in json.loads(js.read_text())["systems"]:
                an = s.get("analysis", {})
                if s.get("system") not in ("ccbar", "bbar") or "value" not in an:
                    continue
                kv = dict(zip(an["names"], an["value"]))
                if tag.endswith("-s") and "sigma0" in kv:
                    # same sigma12 with sigma0 -> 0.1: s' = sqrt(sigma12^2 - 0.01)/m
                    m = kv["m1"]
                    s12 = kv.pop("sigma0") ** 2 + (1.55 * m) ** 2
                    kv["s"] = math.sqrt(max(s12 - 0.01, 1e-6)) / m
                if "fitmu" in tag.split("-") and "mu_s" not in kv:
                    kv["mu_s"] = 1.2e-4     # B -> A, D -> C at the boundary
                if any(sy == s["system"] and kv.keys() == old.keys() and
                       all(abs(kv[k] - old[k]) <= 1e-3 * abs(old[k]) + 1e-9
                           for k in kv) for sy, old in seen):
                    continue                # the same minimum twice
                seen.append((s["system"], kv))
                args += ["--start", s["system"] + ":" +
                         ",".join(f"{k}={v:.10g}" for k, v in kv.items())]
    return args


def floor_shapes():
    """Relative floors w_i from the reference model, per system: proportional
    to <v^2>_i ("v2"), to the neglected kinetic term <p^4>_i/8m^3 ("p4") or to
    r_rms,i ("r"); normalised to mean w^2 = 1 over the fitted states, so that
    sigma_th stays the rms floor. Every state (fitted or not) gets a weight."""
    js = RES / "main" / f"{REFERENCE}.json"
    if not js.exists():
        return {}
    key = {"v2": "v2", "p4": "p4_term_GeV", "r": "r_rms_fm"}
    out = {sh: {} for sh in SHAPES}
    for s in json.loads(js.read_text())["systems"]:
        if s["system"] not in ("ccbar", "bbar"):
            continue
        for sh in SHAPES:
            raw = {x["label"]: x[key[sh]] for x in s["states"]}
            fit = [x[key[sh]] for x in s["states"] if x["in_fit"]]
            scale = math.sqrt(sum(v * v for v in fit) / len(fit))
            out[sh][s["system"]] = {k: v / scale for k, v in raw.items()}
    return out


def shape_args(shape):
    sh = floor_shapes().get(shape, {})
    return [a for sy, w in sh.items() for a in
            ("--floor-weights", sy + ":" + ",".join(f"{k}={v:.6g}"
                                                    for k, v in w.items()))]


def rel_configurations():
    """A/B with a running (smeared) coupling and spin terms, c = 0."""
    return [(t, a) for t, a in configurations()
            if t[0] in "AB" and t.split("-")[1] in ("gi", "frozen", "apt")
            and "-coff-" in t and t.split("-")[3] in ("first", "smeared")
            and not t.endswith("fitmu")]


def apriori_floors():
    """The a-priori floors (MeV) per system: the reference model's truncation
    estimate, i.e. the part of the omitted -<p^4>/8m^3 term its parameters
    cannot absorb (qgem's truncation test). A common floor fixed before the
    model comparison, not tuned to it."""
    js = RES / "main" / f"{REFERENCE}.json"
    if not js.exists():
        raise SystemExit(f"a-priori floors: {js} missing (run the main sweep first)")
    r = json.loads(js.read_text())
    out = {s["system"]: 1e3 * s["analysis"]["truncation"]["sigma_a_GeV"]
           for s in r["systems"] if s["system"] in SYSTEMS
           and "truncation" in s.get("analysis", {})}
    if set(out) != set(SYSTEMS):
        raise SystemExit(f"a-priori floors: no truncation test in {js} for "
                         f"{sorted(set(SYSTEMS) - set(out))}")
    return out


def joint_shared(tag):
    """Flavour-independent parameters beyond b, where the model has one: the
    normalisation of a running NR coupling (alpha_s(r) is one function of r)
    and the GI smearing slope s. None for constant couplings."""
    fam, run = tag[0], tag.split("-")[1]
    if fam in "AB" and run in ("apt", "frozen", "gi"):
        return "b,alpha_s"
    if fam in "CD" and tag.endswith("-s"):
        return "b,s"
    return None


def restarts(tag, nr, gi):
    return ["--multistart", str(gi if tag[0] in "CD" else nr)]


def stats_args(tag):
    """Bootstrap/MCMC sizes: a relativized fit costs ~100x more."""
    if tag[0] in "CD":
        return ["--bootstrap", "100", "--mcmc", "1500", "--threads", "8", "--loo"]
    return ["--bootstrap", "400", "--mcmc", "6000", "--threads", "4", "--loo"]


def _cfg():
    return configurations()


SWEEPS = {
    "lm": lambda: [(t, a + ["--multistart", "8"]) for t, a in _cfg()],
    "minuit": lambda: [(t, a + ["--fitter", "minuit2", "--minuit-check"])
                       for t, a in _cfg()],
    "main": lambda: [(t, a + restarts(t, 8, 3) + ["--loo"] +
                      seeds(t, ["lm", "minuit"])) for t, a in _cfg()],
    # the conventional fit: PDG errors only (sigma_th = 0), physical set
    "floor0": lambda: [(t, a + ["--floor", "0"] + restarts(t, 4, 0) +
                        seeds(t, ["main"])) for t, a in _cfg() if t in physical_set()],
    "floor": lambda: [(t, a + ["--floor", "10"] + restarts(t, 4, 0) +
                       seeds(t, ["main"])) for t, a in _cfg()],
    "floor2": lambda: (lambda fl: [
        (f"{t}-{sy}", a + ["--systems", sy, "--floor", f"{fl[sy]:.6g}"] +
         restarts(t, 4, 0) + seeds(t, ["main", f"floor2:-{sy}"]))
        for t, a in _cfg() for sy in SYSTEMS])(apriori_floors()),
    "joint": lambda: [(t, a + ["--shared", "b"] + restarts(t, 0, 0) +
                       seeds(t, ["main"])) for t, a in _cfg() if t in physical_set()],
    "joint2": lambda: [(t, a + ["--shared", joint_shared(t)] + restarts(t, 0, 0) +
                        seeds(t, ["main"])) for t, a in _cfg()
                       if t in physical_set() and joint_shared(t)],
    "cv": lambda: [(t, a + ["--lmo", "--kfold", "4", "--learning-curve"] +
                    restarts(t, 0, 0) + seeds(t, ["main"])) for t, a in _cfg()],
    "stats": lambda: (lambda keep: [(t, a + stats_args(t) + seeds(t, ["main"]))
                                    for t, a in _cfg() if t in keep])(
        set(physical_set() + CONTRASTS)),
    "nested": lambda: [(t, a + ["--nested-check", "0.1"] + restarts(t, 0, 0) +
                        seeds(t, ["main"])) for t, a in _cfg()],
    "bopt": lambda: [(t, a + ["--optimize-basis"] + restarts(t, 0, 0) +
                      seeds(t, ["main"])) for t, a in _cfg() if t in REPRESENTATIVE],
    "bgrid": lambda: [(f"{t}-r{r1}-R{rm}", a + ["--basis", f"37,{r1},{rm}"] +
                       restarts(t, 0, 0) + seeds(t, ["main"]))
                      for t, a in _cfg() if t in BGRID_SET
                      for r1 in (0.01, 0.02, 0.05, 0.1, 0.2)
                      for rm in (8, 15, 30, 60, 120)],
    "shape": lambda: [(f"{t}-{sh}", a + shape_args(sh) + restarts(t, 4, 0) +
                       ["--loo", "--lmo"] + seeds(t, ["main"]))
                      for t, a in _cfg() if t in physical_set()
                      for sh in SHAPES],
    "rel": lambda: [(f"{t}-rel", a + ["--kinetic", "rel"] + restarts(t, 8, 0) +
                     ["--loo", "--lmo", "--nested-check", "0.1"] +
                     seeds(t, ["main"])) for t, a in rel_configurations()],
}
ORDER = ["lm", "minuit", "main", "floor0", "floor", "floor2", "joint", "joint2",
         "cv", "stats", "nested", "bopt", "bgrid", "shape", "rel"]


def run_one(qgem, sweep, tag, args):
    d = RES / sweep
    js, txt = d / f"{tag}.json", d / f"{tag}.txt"
    if js.exists():
        return tag, "cached", 0.0
    t0 = time.time()
    p = subprocess.run([str(qgem), *args, "--out", str(js)],
                       capture_output=True, text=True)
    txt.write_text(p.stdout + ("\n# stderr:\n" + p.stderr if p.stderr else ""))
    return tag, "ok" if p.returncode == 0 else f"FAILED ({p.returncode})", \
        time.time() - t0


def jobs_of(sweep, only):
    jobs = [j for j in SWEEPS[sweep]() if re.search(only, j[0])]
    # qgem defaults to Minuit2 in a ROOT build; every sweep except `minuit`
    # uses seeded Levenberg-Marquardt explicitly
    if sweep != "minuit":
        jobs = [(t, ["--fitter", "lm"] + args) for t, args in jobs]
    return sorted(jobs, key=lambda j: j[0][0] not in "CD")  # slow ones first


def run_sweep(qgem, sweep, n_jobs, only=""):
    (RES / sweep).mkdir(parents=True, exist_ok=True)
    jobs = jobs_of(sweep, only)
    if sweep == "floor2":
        # seeds include floor2 results of the nested models: run by depth, and
        # rebuild the job list (it reads those seeds) for every wave
        base = {t for t, _ in configurations()}
        level = lambda t: depth(next(b for b in base if t.startswith(b + "-")
                                     and t[len(b) + 1:] in SYSTEMS))
        waves = [[t for t, _ in jobs if level(t) == d] for d in range(4)]
    else:
        waves = [[t for t, _ in jobs]]
    done, total = 0, len(jobs)
    for wave in waves:
        if not wave:
            continue
        current = dict(jobs_of(sweep, only))
        with ThreadPoolExecutor(n_jobs) as ex:
            futs = [ex.submit(run_one, qgem, sweep, t, current[t]) for t in wave]
            for f in as_completed(futs):
                done += 1
                tag, status, dt = f.result()
                print(f"{sweep:7s} [{done:3d}/{total}] {tag:34s} {status:10s} "
                      f"{dt:7.1f}s", flush=True)


def stop():
    """Stop every running sweep invocation (and its fits)."""
    files = sorted(PIDDIR.glob("sweep*.pid"))
    if not files:
        sys.exit("no sweep is running")
    for f in files:
        pgid = int(f.read_text())
        try:
            os.killpg(pgid, signal.SIGTERM)
            print(f"stopped {f.name} (process group {pgid})")
        except ProcessLookupError:
            print(f"{f.name}: had already finished")
        f.unlink(missing_ok=True)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawTextHelpFormatter)
    ap.add_argument("sweep", choices=ORDER + ["all", "stop"])
    ap.add_argument("-j", "--jobs", type=int, default=max(1, os.cpu_count() - 2))
    ap.add_argument("--only", default="",
                    help="regular expression: run only matching tags")
    a = ap.parse_args()
    if a.sweep == "stop":
        return stop()
    pidfile = PIDDIR / f"sweep-{a.sweep}.pid"
    if pidfile.exists():
        sys.exit(f"`sweep.py {a.sweep}` is already running "
                 f"(pid {pidfile.read_text().strip()}); `sweep.py stop` first")
    src = Path(os.environ.get("QGEM_BIN", BUILT))
    if not src.exists():
        sys.exit(f"build first: {src} not found")
    RES.mkdir(exist_ok=True)
    if os.getpgrp() != os.getpid():    # own process group: `stop` ends it
        try:                           # (already a leader under setsid)
            os.setpgrp()
        except PermissionError:
            pass
    pidfile.write_text(str(os.getpgrp()))
    qgem = RES / f".qgem-{os.getpid()}"
    shutil.copy2(src, qgem)            # frozen for the whole invocation
    try:
        for sweep in (ORDER if a.sweep == "all" else [a.sweep]):
            run_sweep(qgem, sweep, a.jobs, a.only)
    finally:
        qgem.unlink(missing_ok=True)
        pidfile.unlink(missing_ok=True)


if __name__ == "__main__":
    main()
