#!/usr/bin/env python3
"""Checks of examples/medium_dissociation.py against known results (run: python3 examples/test_medium_dissociation.py).

1. Coulomb limit (b -> 0, mu = 0): E_n = -m_r (C_F alpha)^2 / (2 n^2) for the nS levels (n = 1, 2) and the 2P level.
2. Critical screening of the 1S state of the Yukawa potential: the bound state disappears at mu_c = 1.1906 m_r C_F alpha
   (the standard value, mu_c a_B = 1.19).
3. Variational bound: the first-order energy in the vacuum Cornell wave function of a real fit is not below the exact
   screened level, for every mu (needs results/main/A-const-coff-first.json, skipped if absent).
"""
import json, math, os, sys
sys.path.insert(0, os.path.dirname(__file__))
import numpy as np
from medium_dissociation import ExactSolver, grid, radial_from_json, expectation, potential, CF

mr, al, b = 0.7, 0.3, 1e-12
s0 = ExactSolver(0, al, b, mr, nb=120, r1=0.02, rmax=800.0, npts=3500)
s1 = ExactSolver(1, al, b, mr, nb=80, r1=0.02, rmax=120.0)
ry = mr * (CF * al) ** 2 / 2
e = s0.levels(0.0)
assert abs(e[0] / (-ry) - 1) < 1e-4, e[0]
assert abs(e[1] / (-ry / 4) - 1) < 2e-3, e[1]
assert abs(s1.levels(0.0)[0] / (-ry / 4) - 1) < 2e-3
mu_c = None
lo, hi = 0.0, 2.0
for _ in range(40):                      # bisection on the sign of the lowest level of the Yukawa potential
    mid = 0.5 * (lo + hi)
    if s0.levels(mid * mr * CF * al)[0] < 0: lo = mid
    else: hi = mid
assert abs(0.5 * (lo + hi) - 1.1906) < 0.004   # the finite box makes it a lower bound: 1.1778 (rmax 120), 1.1888 (800), 0.5 * (lo + hi)
print("Coulomb levels and the Yukawa critical screening (%.4f) agree" % (0.5 * (lo + hi)))

fit = "results/main/A-const-coff-first.json"
if os.path.exists(fit):
    sysd = json.load(open(fit))["systems"][0]
    th = dict(zip(sysd["theta_names"], sysd["theta"])); mr = th["m1"] / 2
    st = next(x for x in sysd["states"] if x["label"] == "J/psi(1S)")
    r, w = grid(); R, dR = radial_from_json(sysd["basis"], st["wavefunction"]["components"][0], r)
    ex = ExactSolver(0, th["alpha_s"], th["b"], mr)
    for mu in (0.1, 0.3, 0.6, 1.0):
        e1 = (lambda t: t[0] + t[1])(expectation(R, dR, 0, r, w, mr, potential(r, mu, th["alpha_s"], th["b"])))
        assert e1 >= ex.levels(mu)[0] - 1e-6, (mu, e1, ex.levels(mu)[0])
    print("first-order energy is an upper bound of the exact level")
print("all checks passed")
