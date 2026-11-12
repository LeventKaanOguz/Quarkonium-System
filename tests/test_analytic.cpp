// Exact analytic matrix elements against the independent numerical route
// (Gauss-Legendre of the pointwise functions, `local`) and against closed-form
// limits.
#include <gtest/gtest.h>

#include <cmath>
#include <numbers>

#include "me/analytic.hpp"
#include "me/coulomb.hpp"
#include "me/radial.hpp"
#include "me/smeared.hpp"
#include "reference_quadrature.hpp"

using namespace qgem;

namespace {
double maxrel(const Mat& a, const Mat& b) {
  return ((a - b).cwiseAbs().array() /
          (b.cwiseAbs().array() + 1e-12 * b.cwiseAbs().maxCoeff()))
      .maxCoeff();
}
const BasisSpec kS{.n = 8, .r1 = 0.04, .rmax = 25.0};
}  // namespace

TEST(Analytic, InverseQuadraticIntegralMatchesQuadrature) {
  for (double a : {1e-3, 0.4, 7.0, 3e3})
    for (double tau : {0.3, 2.0, 20.0})
      for (int two_s : {1, 2, 3, 4, 5, 6, 9}) {
        const double ref = test::integrate(
            [&](double t) { return std::pow(a + t * t, -0.5 * two_s); }, 0, tau,
            800, 16);
        EXPECT_NEAR(inverse_quadratic_integral(two_s, a, tau) / ref, 1.0, 1e-12)
            << a << " " << tau << " " << two_s;
      }
}

TEST(Analytic, ErfMatrixIsTheClosedFormAndMatchesQuadrature) {
  for (int l : {0, 1, 2}) {
    const Basis B = make_basis(l, kS);
    EXPECT_LT(
        maxrel(erf_over_r(B, B, 1.7),
               local(B, B, [](double r) { return std::erf(1.7 * r) / r; })),
        1e-11)
        << l;
  }
  const Basis B0 = make_basis(0, kS), B2 = make_basis(2, kS),
              B1 = make_basis(1, kS);
  EXPECT_LT(
      maxrel(erf_over_r(B0, B1, 1.7),
             local(B0, B1, [](double r) { return std::erf(1.7 * r) / r; })),
      1e-11);  // l_bar = 1/2
  EXPECT_LT(
      maxrel(erf_over_r(B0, B2, 2.3),
             local(B0, B2, [](double r) { return std::erf(2.3 * r) / r; })),
      1e-11);
}

TEST(Analytic, SmearedCoulombSpinFunctions) {
  for (double tau : {0.9, 2.4})
    for (auto [lb, lk] :
         std::vector<std::pair<int, int>>{{1, 1}, {2, 2}, {0, 2}, {1, 3}}) {
      const Basis A = make_basis(lb, kS), K = make_basis(lk, kS);
      EXPECT_LT(
          maxrel(
              gprime_over_r_matrix(A, K, tau),
              local(A, K, [&](double r) { return coulomb_gp_over_r(tau, r); })),
          1e-10)
          << tau << " " << lb << lk;
      EXPECT_LT(
          maxrel(tensor_matrix(A, K, tau),
                 local(A, K, [&](double r) { return coulomb_tensor(tau, r); })),
          1e-9)
          << tau << " " << lb << lk;
    }
}

TEST(Analytic, SmearedLinearPotential) {
  for (double sig : {0.8, 3.1})
    for (int l : {0, 1, 2}) {
      const Basis B = make_basis(l, kS);
      EXPECT_LT(
          maxrel(smeared_r_matrix(B, B, sig),
                 local(B, B, [&](double r) { return smeared_r(sig, r); })),
          1e-10)
          << sig << l;
      if (l >= 1) {
        EXPECT_LT(maxrel(smeared_r_prime_over_r_matrix(B, B, sig),
                         local(B, B,
                               [&](double r) {
                                 return smeared_r_prime_over_r(sig, r);
                               })),
                  1e-8)
            << sig << l;
      }
    }
}

TEST(Analytic, SmearedExponentialByTheTransferLemmas) {
  for (double sig : {1.2, 3.0})
    for (double mu : {0.1, 0.8})
      for (int l : {0, 1, 2, 3}) {
        const Basis B = make_basis(l, kS);
        const Mat ref =
            local(B, B, [&](double r) { return smeared_exp(sig, mu, r); });
        EXPECT_LT(maxrel(smeared_exp_matrix(B, sig, mu), ref), 1e-9)
            << sig << " " << mu << " " << l;
        if (l >= 1) {
          const Mat refp = local(B, B, [&](double r) {
            return smeared_exp_prime(sig, mu, r) / r;
          });
          EXPECT_LT(maxrel(smeared_exp_prime_over_r_matrix(B, sig, mu), refp),
                    1e-8)
              << sig << " " << mu << " " << l;
        }
      }
}
