// M3: observables against the notes' Test 5 (12.2) and exact identities.
#include <gtest/gtest.h>

#include <cmath>
#include <numbers>

#include "models/alpha_s.hpp"
#include "models/nr_model.hpp"
#include "obs/observables.hpp"

using namespace qgem;

namespace {
const BasisSpec kSpec{.n = 40, .r1 = 0.02, .rmax = 50.0, .eta = 1e-13};
Vec theta_of(double m1, double m2, double als, double b, double c, double sig,
             double mu = 0.1) {
  Vec t(kNRParams);
  t << m1, m2, als, b, c, sig, mu;
  return t;
}
Channel single(int L, int S, int J) {
  Channel c;
  c.J = J;
  c.comps = {{L, S}};
  return c;
}
}  // namespace

TEST(SixJ, KnownValuesAndOrthogonality) {
  EXPECT_NEAR(sixj(1, 1, 1, 1, 1, 1), 1.0 / 6.0, 1e-14);
  EXPECT_NEAR(std::fabs(sixj(0, 1, 1, 1, 1, 1)), 1.0 / 3.0, 1e-14);
  EXPECT_NEAR(sixj(2, 2, 2, 2, 2, 2), -3.0 / 70.0, 1e-14);
  // orthogonality: sum_x (2x+1)(2j+1) {a b x; c d j}{a b x; c d j'} = delta_jj'
  for (int j = 0; j <= 3; ++j) {
    double s = 0;
    for (int x = 0; x <= 4; ++x)
      s += (2 * x + 1) * (2 * j + 1) * std::pow(sixj(1, 2, x, 2, 1, j), 2);
    EXPECT_NEAR(s, j <= 2 ? 1.0 : 0.0, 1e-12)
        << j;  // triads (1,1,j), (2,2,j): j = 0..2
  }
}

TEST(Observables, E1AngularFactorsOfTheNotes) {
  // n3S1 -> n'3PJ: C = (2J+1)/9 ; n3PJ -> n'3S1: C = 1/3 ; n3D1 -> n'3P0,1,2:
  // 2/9, 1/6, 1/90
  for (int Jf : {0, 1, 2})
    EXPECT_NEAR(e1_angular_factor(0, 1, 1, Jf, 1), (2 * Jf + 1) / 9.0, 1e-14)
        << Jf;
  for (int Ji : {0, 1, 2})
    EXPECT_NEAR(e1_angular_factor(1, Ji, 0, 1, 1), 1.0 / 3.0, 1e-14) << Ji;
  EXPECT_NEAR(e1_angular_factor(2, 1, 1, 0, 1), 2.0 / 9.0, 1e-14);
  EXPECT_NEAR(e1_angular_factor(2, 1, 1, 1, 1), 1.0 / 6.0, 1e-14);
  EXPECT_NEAR(e1_angular_factor(2, 1, 1, 2, 1), 1.0 / 90.0, 1e-14);
}

TEST(Test5, ObservablesOfTheSpinIndependentCornellModel) {
  NRModel model({.screened = false, .spin = SpinMode::None});
  MatrixCache cache(kSpec);
  const Vec th = theta_of(1.5, 1.5, 0.5, 0.18, 0.0, 1.1);
  const SolvedChannel s0 = solve_channel(model, single(0, 1, 1), th, cache);
  const SolvedChannel p0 = solve_channel(model, single(1, 1, 1), th, cache);
  const SolvedChannel s0s = solve_channel(model, single(0, 0, 0), th, cache);
  const double eQ = 2.0 / 3.0, mQ = 1.5;
  const double M1S = s0.level(0).M, M2S = s0.level(1).M, M1P = p0.level(0).M;
  EXPECT_NEAR(M1S, 3.2164462, 1e-6);
  EXPECT_NEAR(M2S, 3.8672945, 1e-6);
  EXPECT_NEAR(M1P, 3.6993328, 1e-6);

  EXPECT_NEAR(R0_squared(s0, 0, cache), 1.3796907, 1e-6);
  EXPECT_NEAR(R0_squared(s0, 1, cache), 0.8293325, 1e-6);
  EXPECT_NEAR(1e6 * width_ee(R0_squared(s0, 0, cache), M1S, eQ), 12.6252, 1e-3);
  EXPECT_NEAR(1e6 * width_ee(R0_squared(s0, 1, cache), M2S, eQ), 5.2496, 1e-3);
  EXPECT_NEAR(1e3 * decay_constant(R0_squared(s0, 0, cache), M1S), 640.01,
              5e-3);
  EXPECT_NEAR(1e3 * decay_constant(R0_squared(s0, 1, cache), M2S), 452.53,
              5e-3);
  EXPECT_NEAR(
      1e6 * width_gammagamma_S(R0_squared(s0s, 0, cache), s0s.level(0).M, eQ),
      16.8336, 2e-3);
  EXPECT_NEAR(
      1e6 * width_gammagamma_S(R0_squared(s0s, 1, cache), s0s.level(1).M, eQ),
      6.9994, 2e-3);

  // E1: 2 3S1 -> 1 3PJ gamma  (ratios 1 : 3 : 5), 1 3PJ -> 1 3S1 gamma
  const double r21 = std::fabs(radial_overlap(p0, 0, s0, 1, cache, 1));
  const double k21 = photon_energy(M2S, M1P);
  EXPECT_NEAR(r21, 2.614908, 2e-6);
  EXPECT_NEAR(1e3 * k21, 164.314, 5e-3);
  const double expect21[] = {14.5753, 43.7260, 72.8767};
  for (int J = 0; J <= 2; ++J)
    EXPECT_NEAR(1e6 * width_E1(e1_angular_factor(0, 1, 1, J, 1), eQ, k21, r21),
                expect21[J], 2e-3)
        << J;
  const double r10 = radial_overlap(s0, 0, p0, 0, cache, 1);
  const double k10 = photon_energy(M1P, M1S);
  EXPECT_NEAR(std::fabs(r10), 1.940281,
              2e-6);  // overall eigenvector phase is arbitrary
  EXPECT_NEAR(1e3 * k10, 451.370, 5e-3);
  for (int J = 0; J <= 2; ++J)
    EXPECT_NEAR(1e6 * width_E1(e1_angular_factor(1, J, 0, 1, 1), eQ, k10,
                               std::fabs(r10)),
                499.032, 5e-2)
        << J;

  // M1: 1 3S1 -> 1 1S0 gamma with the model's own masses of the
  // spin-independent calculation: the notes use k = 90.9611 MeV (from Delta_hf
  // = 92.2753 MeV of Test 2) and <j0(kr/2)> = 0.99864790
  const double k = 0.0909611;
  const double j0 = radial_overlap_j0(s0, 0, s0, 0, cache, k);
  EXPECT_NEAR(j0, 0.99864790, 2e-7);
  EXPECT_NEAR(1e6 * width_M1(eQ, mQ, k, j0), 1.4426, 1e-3);

  // <v^2> per quark, 1S
  EXPECT_NEAR(p2_expectation(s0, 0, cache) / (mQ * mQ), 0.2806, 5e-5);

  // <p^4> from <(E-V)^2>: 1S, 2S
  const Basis& B = cache.basis(0);
  const double kappa = 2.0 / 3.0, b = 0.18, mr = 0.75;
  const Mat V = -kappa * r_power(B, B, -1) + b * r_power(B, B, 1);
  const Mat V2 = kappa * kappa * r_power(B, B, -2) -
                 2 * kappa * b * overlap(B) + b * b * r_power(B, B, 2);
  EXPECT_NEAR(p4_from_potential(s0.level(0).c, s0.level(0).E, V, V2, mr),
              1.16717, 2e-5);
  EXPECT_NEAR(p4_from_potential(s0.level(1).c, s0.level(1).E, V, V2, mr),
              1.43262, 2e-5);
}

TEST(Observables, HydrogenicWavefunctionAtOriginAndKatoCusp) {
  // Pure Coulomb: |R(0)|^2 = 4 (mr kappa)^3 (1S);  P wave R'(0)^2 for 2P: (mr
  // kappa)^5 / 24 ... check 1S only with a basis that resolves the cusp.
  const double mr = 0.75, kappa = 2.0 / 3.0;
  const Basis B = make_basis(0, {.n = 60, .r1 = 0.003, .rmax = 80.0});
  const Vec c = solve_gep(kinetic_nr(B, mr) - kappa * r_power(B, B, -1),
                          overlap(B), 1e-13)
                    .C.col(0);
  const double R0 = c.dot(origin_coefficients(B));
  EXPECT_NEAR(R0 * R0 / (4 * std::pow(mr * kappa, 3)), 1.0, 2e-3);
  // <r> = 3/(2 mr kappa) for 1S hydrogen-like; <1/r> = mr kappa (exact via
  // virial/Feynman-Hellmann)
  EXPECT_NEAR(c.dot(r_power(B, B, -1) * c) / (mr * kappa), 1.0, 1e-6);
  EXPECT_NEAR(c.dot(r_power(B, B, 1) * c) * 2 * mr * kappa / 3, 1.0, 1e-5);
}

TEST(Observables, OriginDerivativesOfPAndDStates) {
  // Linear + Coulomb: compare R'(0) of 1P and R''(0) of 1D with the
  // finite-difference of the radial function.
  NRModel model({.screened = false, .spin = SpinMode::None});
  MatrixCache cache(kSpec);
  const Vec th = theta_of(1.5, 1.5, 0.5, 0.18, 0.0, 1.1);
  const SolvedChannel p = solve_channel(model, single(1, 0, 1), th, cache);
  const SolvedChannel d = solve_channel(model, single(2, 0, 2), th, cache);
  const Basis &B1 = cache.basis(1), &B2 = cache.basis(2);
  const Vec c1 = p.level(0).c, c2 = d.level(0).c;
  const double h =
      1e-5;  // R = A r^L (1 + O(r^2)): the finite-difference error is O(h^2)
  EXPECT_NEAR(origin_derivative(p, 0, cache, 1), radial_value(B1, c1, h) / h,
              1e-8 * std::fabs(origin_derivative(p, 0, cache, 1)));
  EXPECT_NEAR(origin_derivative(d, 0, cache, 2),
              2 * radial_value(B2, c2, h) / (h * h),
              1e-7 * std::fabs(origin_derivative(d, 0, cache, 2)));
  // normalization: Int r^2 R^2 = 1 = <r^0>
  EXPECT_NEAR(moment(p, 0, cache, 0), 1.0, 1e-12);
}

TEST(Widths, QcdCouplingFactorsAndMixedLeptonicWidth) {
  // MS-bar coupling: the input at M_Z, asymptotic freedom, two-loop values
  EXPECT_NEAR(alpha_s_msbar(91.1876), 0.1180, 1e-10);
  const double ab = alpha_s_msbar(4.18), ac = alpha_s_msbar(1.5);
  EXPECT_GT(ab, 0.215);
  EXPECT_LT(ab, 0.235);
  EXPECT_GT(ac, 0.32);
  EXPECT_LT(ac, 0.38);
  EXPECT_GT(alpha_s_msbar(1.0), ac);
  // one-loop factors
  const double a = 0.3, pi = std::numbers::pi;
  EXPECT_DOUBLE_EQ(qcd_factor("ee", a), 1 - 16 * a / (3 * pi));
  EXPECT_DOUBLE_EQ(qcd_factor("gg_P2", a), qcd_factor("ee", a));
  EXPECT_NEAR(qcd_factor("gg_S", a), 1 - 3.3767 * a / pi, 1e-4);
  EXPECT_NEAR(qcd_factor("gg_P0", a), 1 + 0.1788 * a / pi, 1e-4);
  EXPECT_THROW(qcd_factor("xx", a), std::invalid_argument);
  // the mixed S-D leptonic width reduces to the pure ones
  const double M = 3.77, m = 1.5, e = 2.0 / 3, R0 = 0.8, Rpp = 0.3;
  EXPECT_NEAR(width_ee_SD(R0, 0, M, m, e), width_ee(R0 * R0, M, e), 1e-18);
  EXPECT_NEAR(width_ee_SD(0, Rpp, M, m, e), width_ee_D(Rpp * Rpp, M, m, e),
              1e-18);
}

TEST(Observables, WaveFunctionNormNodesAndOriginValue) {
  // Cornell S and P waves: unit norm, n radial nodes for the n-th level, and
  // the coefficients reproduce R(0) through the basis
  NRModel model({.screened = false, .spin = SpinMode::None});
  MatrixCache cache(kSpec);
  const Vec th = theta_of(1.5, 1.5, 0.5, 0.18, 0.0, 1.1);
  const SolvedChannel s0 = solve_channel(model, single(0, 1, 1), th, cache);
  const SolvedChannel p0 = solve_channel(model, single(1, 1, 1), th, cache);
  for (int n = 0; n < 4; ++n) {
    const WaveFunction w = wave_function(s0, n, cache);
    ASSERT_EQ(w.c.size(), 1u);
    EXPECT_EQ(w.L[0], 0);
    EXPECT_NEAR(w.norm, 1.0, 1e-10) << n;
    EXPECT_EQ(w.nodes, n) << n;
    const double R0 = w.c[0].dot(origin_coefficients(cache.basis(0)));
    EXPECT_NEAR(R0 * R0, R0_squared(s0, n, cache), 1e-12) << n;
  }
  for (int n = 0; n < 3; ++n)
    EXPECT_EQ(wave_function(p0, n, cache).nodes, n) << n;
}

TEST(Observables, VirialTheoremOfTheCornellEigenstates) {
  // exact eigenstates of p^2/2m_r - 4 a/3r + b r obey 2<T> = <r V'> =
  // (4a/3)<1/r> + b<r>; the GEM states do to the basis' completeness
  // (observed: <= 3e-6 for n <= 1, 7e-5 for the third D level)
  NRModel model({.screened = false, .spin = SpinMode::None});
  MatrixCache cache(kSpec);
  const double m = 1.5, a = 0.5, b = 0.18, mr = m / 2;
  const Vec th = theta_of(m, m, a, b, 0.0, 1.1);
  for (int L : {0, 1, 2}) {
    const SolvedChannel s =
        solve_channel(model, single(L, 1, L + 1), th, cache);
    for (int n = 0; n < 3; ++n) {
      const double T = p2_expectation(s, n, cache) / (2 * mr);
      const double rV =
          4 * a / 3 * moment(s, n, cache, -1) + b * moment(s, n, cache, 1);
      // first order in the wave-function error (the energy is second
      // order): 1e-6 for ground states, 1e-4 up to the third D level
      EXPECT_NEAR(2 * T / rV, 1.0, n == 0 ? 1e-6 : 1e-4)
          << "L = " << L << ", n = " << n;
    }
  }
}
