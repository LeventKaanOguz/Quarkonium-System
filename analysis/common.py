"""Shared helpers of the analysis scripts: paths, loading, labels, output.

$QGEM_RESULTS and $QGEM_REPORT override the input and output directories
(results/ and docs/report/), e.g. to analyse an archived sweep.
"""
import csv
import json
import math
import os
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
RES = Path(os.environ.get("QGEM_RESULTS", REPO / "results"))
REPORT = Path(os.environ.get("QGEM_REPORT", REPO / "docs" / "report"))
FIG = REPORT / "figures"
DATA = REPORT / "data"


def load(sweep):
    out = {}
    for f in sorted((RES / sweep).glob("*.json")):
        try:
            out[f.stem] = json.loads(f.read_text())
        except json.JSONDecodeError:
            print("unreadable:", f)
    return out


def system(run, name):
    for s in run["systems"]:
        if s["system"] == name:
            return s
    return None


def describe(tag):
    """Human label of a configuration tag, e.g. A-gi-coff-first."""
    p = tag.split("-")
    model, running = p[0], p[1]
    rn = {"const": "const", "gi": "GI run", "frozen": "frozen",
          "apt": "APT", "free": "free run"}[running]
    s = f"{model} {rn}"
    if "con" in p:
        s += " +c"
    for sp in ("none", "smeared"):
        if sp in p:
            s += f" {sp}"
    if "fitmu" in p:
        s += " fit#mu_{s}"
    if p[-1] == "s" and model in "CD":
        s += " (s free)"
    return s


def plain(tag):
    return describe(tag).replace("#mu_{s}", "mu_s")


def tlabel(name):
    """PDG name -> ROOT TLatex: eta_c(1S) -> #eta_{c}(1S)."""
    import re
    m = re.match(r"^([A-Za-z]+)(?:_([A-Za-z0-9]+))?(.*)$", name)
    if not m:
        return name
    base, sub, rest = m.groups()
    greek = {"eta": "#eta", "chi": "#chi", "psi": "#psi", "Upsilon": "#Upsilon"}
    if base == "J/psi" or name.startswith("J/psi"):
        return name.replace("psi", "#psi")
    b = greek.get(base, base)
    return b + (f"_{{{sub}}}" if sub else "") + rest


def rms(xs):
    xs = [x for x in xs if x is not None and math.isfinite(x)]
    return math.sqrt(sum(x * x for x in xs) / len(xs)) if xs else float("nan")


def write_csv(name, rows, keys=None):
    if not rows:
        return
    keys = keys or list(dict.fromkeys(k for r in rows for k in r))  # union, in order
    DATA.mkdir(parents=True, exist_ok=True)
    with open(DATA / name, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=keys, extrasaction="ignore")
        w.writeheader()
        for r in rows:
            w.writerow({k: (f"{v:.6g}" if isinstance(v, float) else v)
                        for k, v in r.items() if k in keys})


def md_table(rows, cols, heads, fmt):
    s = "| " + " | ".join(heads) + " |\n|" + "|".join("---" for _ in heads) + "|\n"
    for r in rows:
        s += "| " + " | ".join(f(r.get(c)) if callable(f) else
                               fmt_val(r.get(c, float("nan")), f)
                               for c, f in zip(cols, fmt)) + " |\n"
    return s


def fmt_val(v, f):
    if isinstance(v, float) and not math.isfinite(v):
        return "–"
    return format(v, f) if f else str(v)
