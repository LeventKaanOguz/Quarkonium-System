#!/usr/bin/env python3
"""Compare two results trees (for example an archived release and a fresh rerun of `sweep.py all`).

    python3 analysis/compare_runs.py results-archive-extracted/results results [--sweeps main,joint,...]

For every JSON file present in both trees (results/<sweep>/<tag>.json) every numeric leaf is compared. Reported per
sweep: the number of files, the files missing on either side, and the largest relative difference of the fitted
parameters (`theta`), of the masses (`states[].model_GeV`, in MeV), of the theory floors (`sigma_th_GeV`) and
of every other numeric field. Differences below --tol (relative, default 1e-6) count as identical.
"""
import argparse, json, math
from pathlib import Path


def leaves(x, path=""):
    if isinstance(x, dict):
        for k, v in x.items():
            yield from leaves(v, f"{path}/{k}")
    elif isinstance(x, list):
        for i, v in enumerate(x):
            yield from leaves(v, f"{path}[{i}]")
    elif isinstance(x, (int, float)) and not isinstance(x, bool):
        yield path, float(x)


def rel(a, b):
    if a == b:
        return 0.0
    if math.isnan(a) and math.isnan(b):
        return 0.0
    return abs(a - b) / max(abs(a), abs(b), 1e-12)


def classify(path):
    if "/theta[" in path or path.endswith("/theta"):
        return "parameters"
    if "model_GeV" in path:
        return "masses"
    if "sigma_th_GeV" in path:
        return "floors"
    return "other"


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("old"); ap.add_argument("new"); ap.add_argument("--sweeps", default="")
    ap.add_argument("--tol", type=float, default=1e-6)
    a = ap.parse_args()
    old, new = Path(a.old), Path(a.new)
    want = set(a.sweeps.split(",")) if a.sweeps else None
    sweeps = sorted(p.name for p in new.iterdir() if p.is_dir() and (want is None or p.name in want))
    print(f"{'sweep':<10}{'files':>6}{'missing':>8}{'differ':>7}   max rel. difference: parameters / masses / floors / other  (worst file)")
    tot_files = tot_diff = 0
    for sw in sweeps:
        nf = nd = 0
        miss = []
        worst = {"parameters": 0.0, "masses": 0.0, "floors": 0.0, "other": 0.0}
        wf = ("", 0.0)
        for f in sorted((new / sw).glob("*.json")):
            g = old / sw / f.name
            if not g.exists():
                miss.append(f.name); continue
            nf += 1
            try:
                A = dict(leaves(json.load(open(g)))); B = dict(leaves(json.load(open(f))))
            except Exception as e:
                miss.append(f"{f.name} ({e})"); continue
            fm = 0.0; differ = False
            for k in set(A) & set(B):
                r = rel(A[k], B[k])
                c = classify(k)
                worst[c] = max(worst[c], r)
                fm = max(fm, r)
                differ = differ or r > a.tol
            if set(A) != set(B):
                differ = True
            if differ:
                nd += 1
            if fm > wf[1]:
                wf = (f.name, fm)
        tot_files += nf; tot_diff += nd
        print(f"{sw:<10}{nf:>6}{len(miss):>8}{nd:>7}   " + " / ".join(f"{worst[c]:.1e}" for c in ("parameters", "masses", "floors", "other")) + f"   ({wf[0]})")
        for m in miss[:3]:
            print("            not in the old tree:", m)
    print(f"total {tot_files} files compared, {tot_diff} differ by more than {a.tol:g} (relative)")


if __name__ == "__main__":
    main()
