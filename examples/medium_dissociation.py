#!/usr/bin/env python3
"""Quarkonium in a hot medium from the GEM wave functions of a qgem fit.

Application of the exported wave functions (qgem --wavefunctions, or any JSON of results/): for every S and P state
of a fit, the radial function R(r) = sum_n c_n N_nL r^L exp(-nu_n r^2) is rebuilt from the coefficients in the JSON, and
the vacuum state is put into the Debye-screened static potential of the Karsch-Mehr-Satz form,

    V_T(r; mu) = -C_F alpha_s exp(-mu r)/r + (b/mu) (1 - exp(-mu r)),        mu = Debye screening mass,

with the alpha_s and b of the fit (so that V_T -> -C_F alpha_s/r + b r for mu -> 0, the vacuum Cornell potential).
The continuum starts at 2 m_Q + V_T(infinity) = 2 m_Q + b/mu. The binding energy of a state is

    B(mu) = b/mu - E(mu),      E = (kinetic + potential) energy above 2 m_Q,

and the state is dissolved where B(mu*) = 0. Two evaluations of E(mu) are made:

  first order  E1(mu) = <T> + <V_T(mu)> in the vacuum wave function of the JSON (no new diagonalization). By the
               variational principle this is an upper bound of the exact lowest level of its channel, so B1 <= B
               and mu*_1 <= mu*_exact for the lowest state of every L: the estimate is a lower bound of mu*;
  exact        the same Hamiltonian diagonalized in a new, larger Gaussian basis (the numerical check).

The spin-dependent terms of the vacuum fit are not carried into the medium: the Hamiltonian is T + V_T only, and the
states are labelled by (L, number of nodes). The kinetic term is p^2/(2 m_r), as in Model A. Only a Model A fit with
constant coupling is accepted (alpha_s and b enter as numbers). Temperatures are indicative: at leading order
mu = sqrt(N_c/3 + n_f/6) g T, and --g and --nf set g and n_f.

Usage:  python3 examples/medium_dissociation.py results/main/A-const-coff-first.json [--out DIR] [--g 2.0] [--nf 3]
Writes DIR/medium-<system>.json, DIR/medium-<system>.txt (also printed) and DIR/medium-<system>.pdf.
Needs numpy, scipy and matplotlib.
"""
import argparse, json, math, sys
from pathlib import Path
import numpy as np
from scipy.linalg import eigh
from scipy.optimize import brentq

CF = 4.0 / 3.0


def grid(rmin=1e-3, rmax=80.0, n=4000):
    """Gauss-Legendre nodes in ln r; returns r and the weights of  int f(r) r^2 dr."""
    x, w = np.polynomial.legendre.leggauss(n)
    a, b = math.log(rmin), math.log(rmax)
    lr = 0.5 * (b - a) * x + 0.5 * (b + a)
    r = np.exp(lr)
    return r, 0.5 * (b - a) * w * r ** 3     # dr = r dx, times r^2


def radial_from_json(basis, comp, r):
    """R(r) and R'(r) of one component of a JSON wave function."""
    L = comp["L"]
    nu = np.array(basis["nu_GeV2"]); N = np.array(basis["norm_by_L"][str(L)]); c = np.array(comp["c"])
    g = np.exp(-nu[None, :] * r[:, None] ** 2)
    phi = g * (N * c)[None, :] * (r[:, None] ** L)
    R = phi.sum(1)
    dR = (phi * (L / r[:, None] - 2.0 * nu[None, :] * r[:, None])).sum(1)
    return R, dR


def potential(r, mu, alpha, b):
    if mu <= 0:
        return -CF * alpha / r + b * r
    return -CF * alpha * np.exp(-mu * r) / r + (b / mu) * (-np.expm1(-mu * r))


def expectation(R, dR, L, r, w, mr, V):
    """(<T>, <V>) of a radial function, normalized on the grid."""
    n2 = np.sum(w * R * R)
    T = np.sum(w * (dR * dR + L * (L + 1) * R * R / r ** 2)) / (2.0 * mr) / n2
    return T, np.sum(w * R * R * V) / n2


class ExactSolver:
    """Lowest levels of T + V_T in a geometric Gaussian basis for one L (generalized eigenproblem, canonical
    orthogonalization). The overlap and kinetic matrices are built once; V_T is integrated for every mu on the grid."""

    def __init__(self, L, alpha, b, mr, nb=70, r1=0.01, rmax=60.0, nlev=4, npts=1500):
        self.L, self.alpha, self.b, self.nlev = L, alpha, b, nlev
        self.r, self.w = grid(1e-3, 3.0 * rmax, npts)
        r, w = self.r, self.w
        nu = 1.0 / (r1 * (rmax / r1) ** (np.arange(nb) / (nb - 1))) ** 2
        phi = np.exp(-nu[None, :] * r[:, None] ** 2) * r[:, None] ** L
        dphi = phi * (L / r[:, None] - 2.0 * nu[None, :] * r[:, None])
        S = phi.T @ (w[:, None] * phi)
        d = 1.0 / np.sqrt(np.diag(S)); S = S * d[:, None] * d[None, :]
        self.phi = phi * d[None, :]; dphi = dphi * d[None, :]
        T = (dphi.T @ (w[:, None] * dphi) + L * (L + 1) * self.phi.T @ ((w / r ** 2)[:, None] * self.phi)) / (2.0 * mr)
        lam, U = eigh(S)
        keep = lam > 1e-12
        self.X = U[:, keep] / np.sqrt(lam[keep])
        self.Tp = self.X.T @ T @ self.X
        self.cache = {}

    def levels(self, mu):
        key = round(mu, 9)
        if key not in self.cache:
            Vm = self.phi.T @ ((self.w * potential(self.r, mu, self.alpha, self.b))[:, None] * self.phi)
            self.cache[key] = np.linalg.eigvalsh(self.Tp + self.X.T @ Vm @ self.X)[: self.nlev]
        return self.cache[key]


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("fit"); ap.add_argument("--system", default="ccbar,bbar")
    ap.add_argument("--out", default="results/medium"); ap.add_argument("--g", type=float, default=2.0)
    ap.add_argument("--nf", type=int, default=3)
    ap.add_argument("--mumax", type=float, default=2.5)
    a = ap.parse_args()
    data = json.load(open(a.fit)); cfg = data["config"]
    if cfg["model"] != "A" or cfg["running"] != "const":
        sys.exit("medium_dissociation.py needs a Model A fit with constant coupling (alpha_s, b are used as numbers)")
    out = Path(a.out); out.mkdir(parents=True, exist_ok=True)
    mgrid = np.linspace(0.05, a.mumax, 100)
    fac = math.sqrt(1.0 + a.nf / 6.0) * a.g
    for sysd in data["systems"]:
        name = sysd["system"]
        if name not in a.system.split(","):
            continue
        th = dict(zip(sysd["theta_names"], sysd["theta"]))
        m1, m2, al, b = th["m1"], th["m2"], th["alpha_s"], th["b"]
        mr = m1 * m2 / (m1 + m2)
        lines = [f"# medium_dissociation {a.fit} system={name}: m={m1:.4f} GeV, alpha_s={al:.4f}, b={b:.4f} GeV^2, "
                 f"g={a.g}, n_f={a.nf}",
                 f"{'state':<14}{'L':>2}{'nodes':>6}{'r_rms fm':>9}{'E vac':>8}{'E exact':>8}{'mu*_1':>8}{'mu*':>8}"
                 f"{'1/(r_rms mu*)':>14}{'T* (GeV)':>9}"]
        res = []
        solvers = {L_: ExactSolver(L_, al, b, mr) for L_ in (0, 1)}
        for st in sysd["states"]:
            comps = st["wavefunction"]["components"]
            if len(comps) != 1 or comps[0]["L"] > 1 or not st["in_fit"]:
                continue       # S and P states of the fit set (single-component, first-order spin treatment)
            comp = comps[0]; L = comp["L"]; nodes = st["wavefunction"]["nodes"]
            r, w = grid()
            R, dR = radial_from_json(sysd["basis"], comp, r)
            n2 = np.sum(w * R * R)
            r2 = np.sum(w * R * R * r * r) / n2
            Tm, V0 = expectation(R, dR, L, r, w, mr, potential(r, 0.0, al, b))

            def e1(mu):
                return Tm + expectation(R, dR, L, r, w, mr, potential(r, mu, al, b))[1]

            def ex(mu):
                return solvers[L].levels(mu)[nodes]

            def bind(f):
                return lambda mu: b / mu - f(mu)

            # first-order and exact melting masses: first sign change of B(mu) on the scan, refined by Brent
            def root(f):
                vals = [bind(f)(m) for m in mgrid]
                for i in range(len(mgrid) - 1):
                    if vals[i] > 0 >= vals[i + 1]:
                        return brentq(bind(f), mgrid[i], mgrid[i + 1], xtol=1e-5)
                return float("nan")

            m1st = root(e1)
            mex = root(ex)
            Evac = Tm + V0
            Eex0 = solvers[L].levels(0.0)[nodes]
            rrms = math.sqrt(r2) * 0.1973269804
            T_star = mex / fac if mex == mex else float("nan")
            size = 1.0 / (math.sqrt(r2) * mex) if mex == mex else float("nan")
            lines.append(f"{st['label']:<14}{L:>2}{nodes:>6}{rrms:>9.3f}{Evac:>8.3f}{Eex0:>8.3f}{m1st:>8.3f}{mex:>8.3f}"
                         f"{size:>14.2f}{T_star:>9.3f}")
            res.append(dict(label=st["label"], L=L, nodes=nodes, r_rms_fm=rrms, E_vac_GeV=Evac, E_vac_exact_GeV=Eex0,
                            mu_star_first_order_GeV=m1st, mu_star_exact_GeV=mex, inv_r_rms_mu_star=size,
                            T_star_GeV=T_star,
                            B_first_order=[b / m - e1(m) for m in mgrid[::5]],
                            B_exact=[b / m - ex(m) for m in mgrid[::5]]))
        text = "\n".join(lines)
        print(text)
        (out / f"medium-{name}.txt").write_text(text + "\n")
        json.dump(dict(fit=a.fit, system=name, mu_GeV=list(mgrid[::5]), states=res, g=a.g, nf=a.nf,
                       params=dict(m=m1, alpha_s=al, b=b)),
                  open(out / f"medium-{name}.json", "w"), indent=1)
        try:
            import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
            fig, ax = plt.subplots(figsize=(4.2, 3.2))
            for k, s_ in enumerate(res):
                col = f"C{k % 10}"
                ax.plot(mgrid[::5], np.array(s_["B_exact"]) * 1e3, color=col, label=s_["label"])
                ax.plot(mgrid[::5], np.array(s_["B_first_order"]) * 1e3, color=col, ls="--")
            ax.axhline(0, color="k", lw=0.6); ax.set_xlabel(r"$\mu$ (GeV)"); ax.set_ylabel("binding energy (MeV)")
            ax.set_ylim(-300, 1500); ax.legend(fontsize=6, ncol=2); fig.tight_layout()
            fig.savefig(out / f"medium-{name}.pdf"); plt.close(fig)
        except Exception as e:    # the numbers do not need the figure
            print("figure skipped:", e)


if __name__ == "__main__":
    main()
