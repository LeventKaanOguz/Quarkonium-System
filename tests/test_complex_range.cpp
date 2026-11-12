// M11: complex-range Gaussians (notes 7.10): matrix elements by analytic
// continuation against quadrature of the real basis functions, and the table
// of 7.10.5 (a mixed real + complex-range basis resolves highly excited S
// states).
#include <gtest/gtest.h>

#include <cmath>
#include <numbers>

#include "basis/basis.hpp"
#include "linalg/gep.hpp"
#include "me/radial.hpp"
#include "reference_fd.hpp"
#include "reference_quadrature.hpp"

using namespace qgem;

namespace {
// reference: N r^l e^{-eta r^2} {cos, sin}(omega eta r^2), N from quadrature of
// r^2 f^2
double fval(const Basis& B, int i, double r) {
  const double g = std::exp(-B.nu(i) * r * r),
               ph = B.omega(i) * B.nu(i) * r * r;
  const double f =
      B.kind[i] == Fn::Real
          ? g
          : (B.kind[i] == Fn::Cos ? g * std::cos(ph) : g * std::sin(ph));
  return std::pow(r, B.L) * f;
}
double ref_integral(const Basis& B, int i, int j,
                    const std::function<double(double)>& w) {
  const double nu = std::min(B.nu(i), B.nu(j)), rc = std::sqrt(45.0 / (2 * nu));
  return test::integrate(
      [&](double r) { return r * r * fval(B, i, r) * w(r) * fval(B, j, r); }, 0,
      rc, 600, 24);
}
}  // namespace

TEST(ComplexRange, NormalizationOverlapAndMomentsMatchQuadrature) {
  for (int l : {0, 1, 2}) {
    const Basis B = make_basis(l, {.n = 3,
                                   .r1 = 0.3,
                                   .rmax = 4.0,
                                   .nc = 3,
                                   .rho1 = 0.8,
                                   .rhomax = 10.0,
                                   .omega = 0.6});
    ASSERT_EQ(B.size(), 3 + 6);
    const Mat S = overlap(B), R1 = r_power(B, B, 1), Rm1 = r_power(B, B, -1),
              G = gaussian(B, B, 0.7);
    for (int i = 0; i < B.size(); ++i)
      for (int j = 0; j < B.size(); ++j) {
        const double ni = B.norm(i) * B.norm(j);
        EXPECT_NEAR(S(i, j),
                    ni * ref_integral(B, i, j, [](double) { return 1.0; }),
                    1e-10)
            << l << " " << i << "," << j;
        EXPECT_NEAR(R1(i, j),
                    ni * ref_integral(B, i, j, [](double r) { return r; }),
                    1e-9)
            << l << " " << i << "," << j;
        EXPECT_NEAR(Rm1(i, j),
                    ni * ref_integral(B, i, j, [](double r) { return 1 / r; }),
                    1e-9)
            << l << " " << i << "," << j;
        EXPECT_NEAR(
            G(i, j),
            ni * ref_integral(B, i, j,
                              [](double r) { return std::exp(-0.7 * r * r); }),
            1e-10)
            << l << " " << i << "," << j;
      }
    for (int i = 0; i < B.size(); ++i) EXPECT_NEAR(S(i, i), 1.0, 1e-12);
    EXPECT_LT((S - S.transpose()).cwiseAbs().maxCoeff(), 1e-13);
  }
}

TEST(ComplexRange, KineticMatrixIsSymmetricAndMatchesTheGradientForm) {
  const int l = 1;
  const Basis B = make_basis(l, {.n = 2,
                                 .r1 = 0.4,
                                 .rmax = 3.0,
                                 .nc = 2,
                                 .rho1 = 1.0,
                                 .rhomax = 6.0,
                                 .omega = 0.5});
  const Mat P = p2(B);
  EXPECT_LT((P - P.transpose()).cwiseAbs().maxCoeff(),
            1e-12 * P.cwiseAbs().maxCoeff());
  // <p^2> = Int r^2 [f_i' f_j' + l(l+1)/r^2 f_i f_j] dr, derivatives by central
  // differences of the basis functions
  const double h = 1e-5;
  auto fn = [&](int i, double r) { return B.norm(i) * fval(B, i, r); };
  auto dfn = [&](int i, double r) {
    return (fn(i, r + h) - fn(i, r - h)) / (2 * h);
  };
  for (int i = 0; i < B.size(); ++i)
    for (int j = 0; j <= i; ++j) {
      const double rc = std::sqrt(45.0 / (2 * std::min(B.nu(i), B.nu(j))));
      const double ref = test::integrate(
          [&](double r) {
            return r * r * dfn(i, r) * dfn(j, r) +
                   l * (l + 1) * fn(i, r) * fn(j, r);
          },
          1e-4, rc, 800, 24);
      EXPECT_NEAR(P(i, j), ref, 2e-6 * std::max(1.0, std::fabs(ref)))
          << i << "," << j;
    }
}

TEST(ComplexRange, RealSetAloneIsUnchangedAndEnlargingNeverRaisesAnyLevel) {
  const BasisSpec real{.n = 12, .r1 = 0.05, .rmax = 20.0, .eta = 1e-13};
  BasisSpec mixed = real;
  mixed.nc = 4;
  mixed.rho1 = 1.0;
  mixed.rhomax = 12.0;
  mixed.omega = 0.5;
  const double mr = 0.75, kappa = 2.0 / 3.0, b = 0.18;
  auto levels = [&](const BasisSpec& sp) {
    const Basis B = make_basis(0, sp);
    return solve_gep(kinetic_nr(B, mr) - kappa * r_power(B, B, -1) +
                         b * r_power(B, B, 1),
                     overlap(B), sp.eta)
        .E;
  };
  const Vec a = levels(real), m = levels(mixed);
  // adding functions is a nested enlargement of the variational space: Cauchy
  // interlacing, every level can only drop
  for (int k = 0; k < 6; ++k) EXPECT_LE(m(k), a(k) + 1e-10) << k;
  // and the mixed basis is no worse than the real one anywhere the real one is
  // poor
  EXPECT_LT(m(5), a(5) - 1e-6);
}

TEST(ComplexRange, MixedBasisResolvesHighlyExcitedSWavesAsInThe710Table) {
  const double mr = 0.75, kappa = 2.0 / 3.0, b = 0.18;
  auto vA = [&](double r) { return -kappa / r + b * r; };
  Vec ref(8);
  for (int k = 0; k < 8; ++k) ref(k) = test::fd_level(vA, 0, mr, 45.0, 9000, k);
  auto levels = [&](const BasisSpec& sp) {
    const Basis B = make_basis(0, sp);
    return solve_gep(kinetic_nr(B, mr) - kappa * r_power(B, B, -1) +
                         b * r_power(B, B, 1),
                     overlap(B), sp.eta)
        .E;
  };
  // notes: real (0.05, 40), n = 40: errors up to 1e-3 .. 2e-3 for 7S, 8S;
  //        real 20 (0.02, 30) + complex 10 x 2 (1, 20), omega = 0.5: errors <=
  //        7e-7
  const Vec real = levels({.n = 40, .r1 = 0.05, .rmax = 40.0, .eta = 1e-13});
  const Vec mixed = levels({.n = 20,
                            .r1 = 0.02,
                            .rmax = 30.0,
                            .eta = 1e-13,
                            .nc = 10,
                            .rho1 = 1.0,
                            .rhomax = 20.0,
                            .omega = 0.5});
  for (int k = 0; k < 8; ++k) {
    EXPECT_GE(mixed(k), ref(k) - 2e-6)
        << "HUM, level "
        << k;  // upper bounds (up to the FD reference accuracy)
    EXPECT_LT(mixed(k) - ref(k), 5e-6) << "mixed basis, level " << k;
  }
  EXPECT_GT(real(7) - ref(7), 1e-4);  // the plain real basis fails for 8S
  EXPECT_GT(real(7) - ref(7), 100 * (mixed(7) - ref(7)));
  EXPECT_GT(real(6) - ref(6), 20 * std::max(mixed(6) - ref(6), 1e-7));
}
