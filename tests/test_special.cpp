// M1a: special functions against independent quadrature / known values.
#include <gtest/gtest.h>

#include <cmath>
#include <numbers>
#include <vector>

#include "math/quadrature.hpp"
#include "math/special.hpp"
#include "reference_quadrature.hpp"

using namespace qgem;

TEST(Special, ErfcxMatchesDefinitionAndAsymptotics) {
  EXPECT_DOUBLE_EQ(erfcx(0.0), 1.0);
  EXPECT_NEAR(erfcx(1.0), std::exp(1.0) * std::erfc(1.0), 1e-15);
  // both sides of the switch to the asymptotic series, against long-double
  // evaluation
  for (double x : {25.9, 26.0, 26.1, 40.0, 100.0}) {
    const long double ref = std::exp(static_cast<long double>(x) * x) *
                            std::erfc(static_cast<long double>(x));
    EXPECT_NEAR(erfcx(x) / static_cast<double>(ref), 1.0, 1e-13) << x;
  }
  const double x =
      30.0;  // asymptotic branch, reference by series to many terms
  EXPECT_NEAR(erfcx(x),
              1.0 / (x * std::sqrt(std::numbers::pi)) *
                  (1 - 1 / (2 * x * x) + 3 / (4 * x * x * x * x)),
              1e-9);
  EXPECT_NEAR(erfcx(-1.0), 2 * std::exp(1.0) - erfcx(1.0), 1e-14);
}

TEST(Special, MasterIntegralClosedForm) {
  // Int r^2 e^{-a r^2} dr = sqrt(pi)/(4 a^{3/2});   Int r^0 e^{-a r^2} dr =
  // sqrt(pi)/(2 sqrt a)
  EXPECT_NEAR(master_integral(2, 3.0) /
                  (std::sqrt(std::numbers::pi) / (4 * std::pow(3.0, 1.5))),
              1.0, 1e-14);
  EXPECT_NEAR(master_integral(0, 0.7) /
                  (std::sqrt(std::numbers::pi) / (2 * std::sqrt(0.7))),
              1.0, 1e-14);
  EXPECT_NEAR(
      master_integral(-0.5, 2.0),
      qgem::test::integrate(
          [](double r) { return std::pow(r, -0.5) * std::exp(-2 * r * r); }, 0,
          6, 4000, 24),
      2e-3);
  EXPECT_THROW(master_integral(-1.0, 1.0), std::invalid_argument);
}

TEST(Special, ScreenedMomentsAgainstQuadrature) {
  // x = mu / (2 sqrt a) from 0.01 to 1000: covers both the upward and the
  // downward branch.
  for (double a : {1e-3, 0.3, 4.0, 2e4}) {
    for (double x : {0.01, 0.3, 0.49, 0.51, 1.0, 3.0, 8.0, 50.0, 1000.0}) {
      const double mu = 2 * x * std::sqrt(a);
      std::vector<double> Q(12);
      screened_moments(11, a, mu, Q);
      const double rc = (-mu + std::sqrt(mu * mu + 4 * a * 60.0)) /
                        (2 * a);  // exponent -a r^2 - mu r = -60
      for (int k : {0, 1, 2, 5, 11}) {
        const double ref = qgem::test::integrate(
            [&](double r) {
              return std::pow(r, k) * std::exp(-a * r * r - mu * r);
            },
            0, rc * 1.6, 600, 24);
        EXPECT_NEAR(Q[k] / ref, 1.0, 2e-12)
            << "a=" << a << " x=" << x << " k=" << k;
      }
    }
  }
}

TEST(Quadrature, GaussLegendreIntegratesPolynomialsExactly) {
  const auto& q = gauss_legendre01(8);
  for (int k = 0; k <= 15; ++k) {
    double s = 0;
    for (int i = 0; i < 8; ++i) s += q.w(i) * std::pow(q.x(i), k);
    EXPECT_NEAR(s, 1.0 / (k + 1), 1e-15) << k;
  }
}

TEST(Quadrature, GaussLaguerreMomentsAreGammaFunctions) {
  // sum w x^k = Gamma(alpha + k + 1) for k <= 2n-1 (n = 40, alpha = 1.5, k <=
  // 10: ~1e-14).
  const auto& q = gauss_laguerre(40, 1.5);
  for (int k = 0; k <= 10; ++k) {
    double s = 0;
    for (int i = 0; i < 40; ++i) s += q.w(i) * std::pow(q.x(i), k);
    EXPECT_NEAR(s / std::tgamma(2.5 + k), 1.0, 1e-13) << k;
  }
}

TEST(Special, SphericalBesselFunctions) {
  // j_l against the explicit closed forms for l <= 2 over several decades,
  // including where std::sph_bessel gives up
  for (double x : {1e-6, 1e-3, 0.1, 0.9, 3.9, 4.1, 7.0, 50.0, 1e3, 3e4}) {
    const long double xl = x, s = std::sin(xl), c = std::cos(xl);
    const long double j0 = s / xl, j1 = s / (xl * xl) - c / xl,
                      j2 = (3 / (xl * xl) - 1) * s / xl - 3 * c / (xl * xl);
    EXPECT_NEAR(sph_bessel_j(0, x), static_cast<double>(j0),
                1e-13 * std::max(1.0, 1 / x))
        << x;
    if (x > 1e-2) {
      EXPECT_NEAR(sph_bessel_j(1, x), static_cast<double>(j1), 1e-12) << x;
      EXPECT_NEAR(sph_bessel_j(2, x), static_cast<double>(j2), 1e-12) << x;
    }
  }
  EXPECT_NEAR(sph_bessel_j(1, 1e-4), 1e-4 / 3, 1e-12);
  EXPECT_NEAR(sph_bessel_j(2, 1e-3), 1e-6 / 15 * (1 - 1e-6 / 14),
              1e-20);  // x^2/15 (1 - x^2/14 + ...)
  EXPECT_DOUBLE_EQ(sph_bessel_j(3, 0.0), 0.0);
  // agrees with the standard library where that works, for l = 3, 4
  for (double x : {0.3, 2.0, 9.0, 20.0})
    for (int l : {3, 4})
      EXPECT_NEAR(sph_bessel_j(l, x), std::sph_bessel(l, x), 1e-13)
          << l << " " << x;
}
