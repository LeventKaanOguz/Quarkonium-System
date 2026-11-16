#!/usr/bin/env python3
"""Re-evaluate the annihilation widths of stored fits after a change in the width code.

The widths are predictions: no fit depends on them (the selection uses the mass-only fits), so a
change in src/obs/observables.cpp or src/fit/widths.cpp leaves every fitted parameter, covariance
and mass untouched. This script re-runs qgem for each stored fit with the floor held at the stored
value and the minimum as its only start, checks that chi^2 and the parameters come back, and
replaces only the "widths" block and the "decay_benchmark" block (a function of the widths) of the stored JSON. Nothing else is modified.

    python3 analysis/refresh_widths.py main,floor2 'A-.*-smeared|B-.*-smeared|[CD]-.*-s$'

The first argument lists the sweeps, the second is a regular expression on the tag.
Used for the sign of the D-wave term of width_ee_SD (v1.1.0)."""
import json
import math
import re
import subprocess
import sys
import tempfile
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import sweep  # noqa: E402

ROOT = sweep.ROOT
RES = sweep.RES
BIN = sweep.BUILT


def base_args(tag):
    return dict(sweep.configurations())[tag]


def start_arg(system, an):
    kv = dict(zip(an["names"], an["value"]))
    return ["--start", system + ":" + ",".join(f"{k}={v:.12g}" for k, v in kv.items())]


def refresh(path, tag, system, floor_GeV):
    stored = json.loads(path.read_text())
    s_old = next(s for s in stored["systems"] if s["system"] == system)
    with tempfile.TemporaryDirectory() as d:
        out = Path(d) / "x.json"
        args = [str(BIN), "--fitter", "lm", *base_args(tag), "--systems", system,
                "--floor", f"{floor_GeV * 1e3:.12g}", *start_arg(system, s_old["analysis"]),
                "--out", str(out)]
        p = subprocess.run(args, capture_output=True, text=True)
        if p.returncode:
            raise SystemExit(f"{tag} {system}: qgem failed\n{p.stderr[-500:]}")
        new = json.loads(out.read_text())
    s_new = next(s for s in new["systems"] if s["system"] == system)
    old_widths = json.loads(json.dumps(s_old["widths"]))
    dchi = abs(s_new["chi2"] - s_old["chi2"])
    dth = max(abs(a - b) / max(abs(b), 1e-12) for a, b in zip(s_new["theta"], s_old["theta"]))
    ok = dchi <= 1e-4 * max(1.0, abs(s_old["chi2"])) and dth <= 1e-4
    if ok:
        s_old["widths"] = s_new["widths"]
        if "decay_benchmark" in s_new["analysis"]:
            s_old["analysis"]["decay_benchmark"] = s_new["analysis"]["decay_benchmark"]
        path.write_text(json.dumps(stored))
    return ok, dchi, dth, old_widths, s_new


def main():
    sweeps = sys.argv[1].split(",")
    pat = re.compile(sys.argv[2])
    for sw in sweeps:
        for js in sorted((RES / sw).glob("*.json")):
            tag = js.stem
            system = None
            if sw == "floor2":
                tag, _, system = tag.rpartition("-")
            if not pat.search(tag):
                continue
            stored = json.loads(js.read_text())
            for s in stored["systems"]:
                if s["system"] not in ("ccbar", "bbar") or (system and s["system"] != system):
                    continue
                if not s["widths"]:
                    continue
                ok, dchi, dth, so, sn = refresh(js, tag, s["system"], s["sigma_th_GeV"])
                moved = [(a["label"], a["process"], a["nlo_keV"], b["nlo_keV"])
                         for a, b in zip(so, sn["widths"])
                         if abs(a["nlo_keV"] - b["nlo_keV"]) > 1e-9 * max(1, abs(a["nlo_keV"]))]
                print(f"{sw:7s} {tag:24s} {s['system']:5s} chi2 {dchi:.1e} theta {dth:.1e} "
                      f"{'ok' if ok else 'MISMATCH (not written)'}  widths changed: {len(moved)}", flush=True)


if __name__ == "__main__":
    main()
