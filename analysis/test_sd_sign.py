"""Independent check of the sign of the D-wave term in the S-D mixed leptonic width.

width_ee_SD (src/obs/observables.cpp) uses A = R_S(0) - 5 R''_D(0) / (2 sqrt2 m^2) for the two
components of one eigenvector in the Condon-Shortley coupling, where <3D1|S12|3S1> = +sqrt8.
This script rebuilds both numbers from scratch, with no reference to the C++ code:

  1. the CS spin-angle function of L = 2, S = 1, J = 1 (Clebsch-Gordan coefficients, spherical
     harmonics, spin-1 vectors e_s) against the Cartesian structure [r_i r_j - delta_ij/3] e_j;
  2. <3D1|S12|3S1> by quadrature on the sphere, with S12 = 2[3 (S.r)^2 - S^2];
  3. the Dirac current u-bar gamma^i v(-p) = sigma^i (1 + p^2/4m^2) - p^i sigma.p / 2m^2 (checked
     numerically against explicit four-spinors) acting on r^2 [r_i r_j - delta_ij r^2/3] e_j,
     whose p-space form is -d_i d_j at the origin.

Run: python3 analysis/test_sd_sign.py (also registered with CTest)."""
import math

import numpy as np

def sph_harm_y(l, m, th, ph):
    """Y_lm for l = 0, 2 (Condon-Shortley phase), written out so that no scipy version matters."""
    if l == 0:
        return 1 / math.sqrt(4 * math.pi) + 0j
    c, s_ = math.cos(th), math.sin(th)
    return {0: 0.25 * math.sqrt(5 / math.pi) * (3 * c * c - 1) + 0j,
            1: -0.5 * math.sqrt(15 / (2 * math.pi)) * s_ * c * np.exp(1j * ph),
            -1: 0.5 * math.sqrt(15 / (2 * math.pi)) * s_ * c * np.exp(-1j * ph),
            2: 0.25 * math.sqrt(15 / (2 * math.pi)) * s_ * s_ * np.exp(2j * ph),
            -2: 0.25 * math.sqrt(15 / (2 * math.pi)) * s_ * s_ * np.exp(-2j * ph)}[m]


CG = {(2, -1, 1): math.sqrt(3 / 5), (1, 0, 1): -math.sqrt(3 / 10), (0, 1, 1): math.sqrt(1 / 10)}
E = {1: -np.array([1, 1j, 0]) / math.sqrt(2), 0: np.array([0, 0, 1.0 + 0j]),
     -1: np.array([1, -1j, 0]) / math.sqrt(2)}


def spin1_matrices():
    eps = np.zeros((3, 3, 3))
    for a, b, c in [(0, 1, 2), (1, 2, 0), (2, 0, 1)]:
        eps[a, b, c], eps[a, c, b] = 1, -1
    return [-1j * eps[a] for a in range(3)]         # (S_a)_bc = -i eps_abc; S_z e_{+1} = +e_{+1}


def d_cs(th, ph):
    """CS spin-angle function (L=2, S=1 -> J=1, M=1) as a Cartesian vector."""
    return sum(CG[(m, s, 1)] * sph_harm_y(2, m, th, ph) * E[s]
               for (m, s, M) in CG if M == 1)


def test_cs_function_is_minus_quadrupole_structure():
    th, ph = 0.7, 1.1
    r = np.array([math.sin(th) * math.cos(ph), math.sin(th) * math.sin(ph), math.cos(th)])
    q_e = r * np.dot(r, E[1]) - E[1] / 3
    ratio = d_cs(th, ph) / q_e
    assert np.allclose(ratio, -math.sqrt(9 / 2) / math.sqrt(4 * math.pi), atol=1e-12)


def test_tensor_matrix_element_is_plus_sqrt8():
    sx = spin1_matrices()
    xs, ws = np.polynomial.legendre.leggauss(40)
    nph, tot = 80, 0
    for x, w in zip(xs, ws):
        th = math.acos(x)
        for k in range(nph):
            ph = 2 * math.pi * (k + 0.5) / nph
            r = np.array([math.sin(th) * math.cos(ph), math.sin(th) * math.sin(ph), x])
            sr = sum(r[a] * sx[a] for a in range(3))
            s12 = 2 * (3 * sr @ sr - 2 * np.eye(3))
            tot += w * 2 * math.pi / nph * np.vdot(d_cs(th, ph), s12 @ (sph_harm_y(0, 0, th, ph) * E[1]))
    assert abs(tot - 2 * math.sqrt(2)) < 1e-12


def test_dirac_current_expansion():
    m, p = 1.0, np.array([0.13, -0.07, 0.11])
    en = math.sqrt(m * m + p @ p)
    sig = [np.array([[0, 1], [1, 0]]), np.array([[0, -1j], [1j, 0]]), np.diag([1.0, -1.0])]
    sp = sum(p[a] * sig[a] for a in range(3))
    # u(p) = (chi, sigma.p/(E+m) chi), v(-p) = (-sigma.p/(E+m) eta, eta); u-bar gamma^i v = u^dag alpha^i v
    for i in range(3):
        alpha = np.block([[np.zeros((2, 2)), sig[i]], [sig[i], np.zeros((2, 2))]])
        for chi in np.eye(2):
            for eta in np.eye(2):
                u = np.concatenate([chi, sp @ chi / (en + m)])
                v = np.concatenate([-sp @ eta / (en + m), eta])
                exact = np.vdot(u, alpha @ v)
                approx = np.vdot(chi, (sig[i] * (1 + p @ p / (4 * m * m)) - p[i] * sp / (2 * m * m)) @ eta)
                assert abs(exact - approx) < 3e-4, (exact, approx)    # O(p^4/m^4)


def test_d_wave_amplitude_sign():
    """A_i = (1/2m^2) d_i d_b Psi_b(0) for Psi_b = (R''/2) N_D (r_b r_c - delta_bc r^2/3) e_c; N_D < 0 in CS."""
    n_d = -math.sqrt(9 / 2) / math.sqrt(4 * math.pi)
    e = E[1]

    def psi(rv):                                                   # R''/2 = 1
        return n_d * (rv * np.dot(rv, e) - e * np.dot(rv, rv) / 3)

    h = 1e-3
    amp = np.zeros(3, dtype=complex)
    for i in range(3):
        for b in range(3):
            ei, eb = np.eye(3)[i] * h, np.eye(3)[b] * h
            dd = (psi(ei + eb)[b] - psi(ei - eb)[b] - psi(-ei + eb)[b] + psi(-ei - eb)[b]) / (4 * h * h)
            amp[i] += dd
    amp *= 1 / (2.0)                                               # 1/(2 m^2) with m = 1; Tr[sigma sigma] = 2 below
    amp *= 2
    s_amp = 2 / math.sqrt(4 * math.pi) * e                         # 2 R_S(0)/sqrt(4 pi) e_i with R_S = 1
    ratio = np.vdot(e, amp)[()] / np.vdot(e, s_amp)[()]            # = A_D / A_S per (R''/2 = 1, R_S = 1)
    # in terms of R'' = 2: A_D/A_S = -5 R'' / (2 sqrt2 m^2 R_S) = -5 * 2 / (2 sqrt2)
    assert abs(ratio.real - (-5 * 2 / (2 * math.sqrt(2)))) < 1e-5


if __name__ == "__main__":
    for name, fn in sorted(globals().items()):
        if name.startswith("test_"):
            fn()
            print("ok", name)
