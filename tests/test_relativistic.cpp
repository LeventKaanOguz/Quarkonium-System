// M8/M9: relativistic kinetic energy (notes Test 3) and non-local sandwiches
// f(p) W(r) f(p) (Test 4).
#include <gtest/gtest.h>

#include <cmath>
#include <numbers>

#include "linalg/gep.hpp"
#include "me/momentum.hpp"
#include "me/radial.hpp"
#include "models/nr_model.hpp"
#include "obs/observables.hpp"
#include "solve/spectrum.hpp"

using namespace qgem;

namespace {
const double m = 1.5, kappa = 2.0 / 3.0, tau = 2.0, b = 0.18;
double E(double p) { return std::sqrt(p * p + m * m); }
}  // namespace

TEST(Momentum, DiagonalMatrixElementsOfPowersAreExact) {
  for (int l : {0, 1, 2}) {
    const Basis B = make_basis(l, {.n = 8, .r1 = 0.1, .rmax = 8.0});
    const Mat p2m = momentum_diagonal(B, [](double p) { return p * p; });
    EXPECT_LT((p2m - p2(B)).cwiseAbs().maxCoeff(),
              1e-11 * p2(B).cwiseAbs().maxCoeff())
        << l;  // Lemma 7.6 via p-space
    const Mat one = momentum_diagonal(B, [](double) { return 1.0; });
    EXPECT_LT((one - overlap(B)).cwiseAbs().maxCoeff(), 1e-14);
    // <p^4> = <p^2 R | p^2 R> = A_i A_j M0 - 4 (A_i nu_j^2 + A_j nu_i^2) M2 +
    // 16 nu_i^2 nu_j^2 M4   (eq. 10.21)
    const Mat p4m =
        momentum_diagonal(B, [](double p) { return p * p * p * p; });
    for (int i = 0; i < 8; ++i)
      for (int j = 0; j < 8; ++j) {
        const double Ai = 2 * B.nu(i) * (2 * l + 3),
                     Aj = 2 * B.nu(j) * (2 * l + 3), a = B.nu(i) + B.nu(j);
        auto Mk = [&](int k) {
          return std::tgamma(l + 1.5 + 0.5 * k) /
                 (2 * std::pow(a, l + 1.5 + 0.5 * k));
        };
        const double ref =
            std::exp(0.5 * (B.lognorm2(i) + B.lognorm2(j))) *
            (Ai * Aj * Mk(0) -
             4 * (Ai * B.nu(j) * B.nu(j) + Aj * B.nu(i) * B.nu(i)) * Mk(2) +
             16 * B.nu(i) * B.nu(i) * B.nu(j) * B.nu(j) * Mk(4));
        EXPECT_NEAR(p4m(i, j) / ref, 1.0, 1e-10) << l << " " << i << "," << j;
      }
  }
}

TEST(Test3, RelativisticKineticEnergyInAnErfCoulombPlusLinearPotential) {
  const Basis B =
      make_basis(0, {.n = 40, .r1 = 0.02, .rmax = 50.0, .eta = 1e-13});
  const Mat S = overlap(B);
  const Mat K = momentum_diagonal(B, [](double p) { return 2 * E(p); }, 80);
  const Mat V = -kappa * erf_over_r(B, B, tau) + b * r_power(B, B, 1);
  const Vec M = solve_gep(K + V, S, 1e-13).E;
  EXPECT_NEAR(M(0), 3.2060876, 3e-7);
  EXPECT_NEAR(M(1), 3.8167930, 3e-7);
  EXPECT_NEAR(M(2), 4.2405517, 3e-7);
  // non-relativistic comparison (same V): E = M - 2m with p^2 / (2 m_r), m_r =
  // m/2
  const Vec Enr = solve_gep(kinetic_nr(B, m / 2) + V, S, 1e-13).E;
  EXPECT_NEAR(Enr(0), 0.2498198, 3e-7);
  EXPECT_NEAR(Enr(1), 0.8870021, 3e-7);
  EXPECT_NEAR(Enr(2), 1.3401780, 3e-7);
  // the relativistic levels lie below 2m + E_NR: 2E(p) <= 2m + p^2/m
  for (int k = 0; k < 3; ++k) EXPECT_LT(M(k), 2 * m + Enr(k));
}

TEST(Test4, SandwichedCoulombByClosureAndByTheTwoStepMethod) {
  // H = 2E(p) + f(p) [-kappa erf(tau r)/r] f(p) + b r,  f = (1 +
  // p^2/E^2)^{1/2};  S waves
  auto f = [](double p) { return std::sqrt(1 + p * p / (E(p) * E(p))); };
  const double dvr[] = {3.0335728, 3.6863728, 4.1247389};
  // closure, n = 30
  {
    const Basis B =
        make_basis(0, {.n = 30, .r1 = 0.02, .rmax = 40.0, .eta = 1e-12});
    const Mat S = overlap(B);
    const Mat X = orthonormalizer(S, 1e-12);
    const Mat K = momentum_diagonal(
                  B, [](double p) { return 2 * E(p); }, 80),
              F = momentum_diagonal(B, f, 80);
    const Mat W = -kappa * erf_over_r(B, B, tau);
    Mat Hp = X.transpose() * K * X + sandwich_closure(F, W, F, X, X) +
             b * X.transpose() * r_power(B, B, 1) * X;
    Hp = 0.5 * (Hp + Hp.transpose());
    Eigen::SelfAdjointEigenSolver<Mat> es(Hp);
    // notes: GEM n = 30 closure 3.0335728, 3.6863728, 4.1247409 (their basis
    // range is not stated; the 3S level moves by ~1e-6 with rmax, see the scan
    // in docs/VALIDATION.md); this basis gives 4.1247402
    const double expect_closure[] = {3.0335728, 3.6863728, 4.1247409};
    for (int k = 0; k < 3; ++k) {
      EXPECT_NEAR(es.eigenvalues()(k), expect_closure[k], 1.5e-6)
          << "closure " << k;
      EXPECT_NEAR(es.eigenvalues()(k), dvr[k], 5e-6)
          << "closure vs DVR " << k;  // documented ~2e-6 closure error
    }
  }
  // a well-placed basis: the closure approaches the DVR value
  {
    const Basis B =
        make_basis(0, {.n = 30, .r1 = 0.1, .rmax = 30.0, .eta = 1e-12});
    const Mat S = overlap(B), X = orthonormalizer(S, 1e-12);
    const Mat K = momentum_diagonal(
                  B, [](double p) { return 2 * E(p); }, 80),
              F = momentum_diagonal(B, f, 80);
    Mat Hp = X.transpose() * K * X +
             sandwich_closure(F, -kappa * erf_over_r(B, B, tau), F, X, X) +
             b * X.transpose() * r_power(B, B, 1) * X;
    Hp = 0.5 * (Hp + Hp.transpose());
    Eigen::SelfAdjointEigenSolver<Mat> es(Hp);
    for (int k = 0; k < 3; ++k)
      EXPECT_NEAR(es.eigenvalues()(k), dvr[k], 3e-7)
          << "closure, converged basis " << k;
  }
  // two-step, n = 30: exact matrix elements of f W f (no closure approximation)
  {
    const Basis B =
        make_basis(0, {.n = 30, .r1 = 0.02, .rmax = 40.0, .eta = 1e-12});
    const Mat S = overlap(B);
    const Mat K = momentum_diagonal(B, [](double p) { return 2 * E(p); }, 80);
    const Mat W2 = sandwich_two_step(
        B, B, f, [](double r) { return -kappa * std::erf(tau * r) / r; }, 45.0);
    EXPECT_LT((W2 - W2.transpose()).cwiseAbs().maxCoeff(),
              1e-8 * W2.cwiseAbs().maxCoeff());
    const Mat H = K + 0.5 * (W2 + W2.transpose()) + b * r_power(B, B, 1);
    const Vec Ev = solve_gep(H, S, 1e-12).E;
    const double expect_two_step[] = {3.0335727, 3.6863728,
                                      4.1247408};  // notes, two-step n = 30
    for (int k = 0; k < 3; ++k) {
      EXPECT_NEAR(Ev(k), expect_two_step[k], 1.5e-6) << "two-step " << k;
      EXPECT_NEAR(Ev(k), dvr[k], 5e-6) << "two-step vs DVR " << k;
    }
  }
}

TEST(SemiRelativistic, MassDerivativeFirstOrderLimitAndHerbstBound) {
  // sum_i (sqrt(p^2 + m_i^2) - m_i) in Model A with the Godfrey-Isgur
  // (smeared) running coupling
  NRModel rel({.running = Running::Gaussian,
               .spin = SpinMode::None,
               .relativistic_kinetic = true});
  NRModel nr({.running = Running::Gaussian, .spin = SpinMode::None});
  const Channel ch{.J = 0, .comps = {{0, 0}}};
  MatrixCache cache({.n = 20, .r1 = 0.05, .rmax = 15.0, .eta = 1e-12});
  Vec th(kNRParams);
  th << 1.5, 1.5, 0.6, 0.18, 0.0, 1.2, 0.1;
  // 1. dH/dm1 against a central finite difference
  const ChannelMatrices m = rel.build(ch, th, cache, true);
  const double h = 1e-5;
  Vec tp = th, tm = th;
  tp(kM1) += h;
  tm(kM1) -= h;
  const Mat fd =
      (rel.build(ch, tp, cache, false).H - rel.build(ch, tm, cache, false).H) /
      (2 * h);
  EXPECT_LT((m.dH[kM1] - fd).cwiseAbs().maxCoeff(), 1e-6 * fd.norm());
  // 2. small velocities: the shift from the relativistic kinetic energy
  // equals the first-order mass-velocity term <p^4>/4m^3 (weak coupling,
  // <v^2> ~ 0.005; the residue is O(v^2) relative)
  Vec weak = th;
  weak(kM1) = weak(kM2) = 10.0;
  weak(kAlphaS) = 0.1;
  const SolvedChannel s_nr = solve_channel(nr, ch, weak, cache);
  const double shift =
      s_nr.level(0).M - solve_channel(rel, ch, weak, cache).level(0).M;
  const double first_order =
      p4_expectation(s_nr, 0, cache) / (4 * std::pow(10.0, 3));
  EXPECT_GT(shift, 0.0);  // the relativistic kinetic energy is softer
  EXPECT_NEAR(shift / first_order, 1.0, 0.05);
  // 3. a point Coulomb potential is refused (unbounded for kappa > 2/pi)
  NRModel point({.running = Running::Constant,
                 .spin = SpinMode::None,
                 .relativistic_kinetic = true});
  EXPECT_THROW(point.build(ch, th, cache, false), std::invalid_argument);
}
