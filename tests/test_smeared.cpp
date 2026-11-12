// Smeared confining functions (notes 6.5.4-6.5.5) against the radial
// convolution formula (6.9) by quadrature:
//   f~(r) = (sigma / (sqrt(pi) r)) Int_{-inf}^{inf} du  u f(|u|) exp(-sigma^2
//   (u - r)^2).
#include <gtest/gtest.h>

#include <cmath>
#include <functional>
#include <numbers>

#include "me/smeared.hpp"
#include "reference_quadrature.hpp"

using namespace qgem;

namespace {
double smear_by_quadrature(const std::function<double(double)>& f, double sigma,
                           double r) {
  const double w = 12.0 / sigma;
  const auto g = [&](double u) {
    return u * f(std::fabs(u)) * std::exp(-sigma * sigma * (u - r) * (u - r));
  };
  // split at u = 0 where |u| has a kink
  return sigma / (std::sqrt(std::numbers::pi) * r) *
         (test::integrate(g, r - w, 0.0, 200, 16) +
          test::integrate(g, 0.0, r + w, 200, 16));
}
}  // namespace

TEST(Smeared, LinearPotentialAndItsDerivative) {
  for (double sigma : {0.7, 1.8, 4.0})
    for (double r : {0.05, 0.4, 1.0, 3.0, 9.0}) {
      const double ref =
          smear_by_quadrature([](double x) { return x; }, sigma, r);
      EXPECT_NEAR(smeared_r(sigma, r), ref, 1e-10 * (1 + ref))
          << sigma << " " << r;
      const double h = 1e-5;
      const double d =
          (smeared_r(sigma, r + h) - smeared_r(sigma, r - h)) / (2 * h);
      EXPECT_NEAR(smeared_r_prime_over_r(sigma, r) * r, d, 1e-7)
          << sigma << " " << r;
    }
  // limits (6.23): r~(0) = 2/(sqrt(pi) sigma), r~'/r -> 4 sigma/(3 sqrt(pi)); r
  // -> inf: r + 1/(2 sigma^2 r)
  const double s = 1.3, sp = std::sqrt(std::numbers::pi);
  EXPECT_NEAR(smeared_r(s, 1e-9), 2 / (sp * s), 1e-8);
  EXPECT_NEAR(smeared_r_prime_over_r(s, 1e-7), 4 * s / (3 * sp), 1e-10);
  EXPECT_NEAR(smeared_r(s, 30.0), 30.0 + 1 / (2 * s * s * 30.0), 1e-12);
  // series/direct switch at s r = 1 is continuous
  EXPECT_NEAR(smeared_r_prime_over_r(s, (1 - 1e-9) / s),
              smeared_r_prime_over_r(s, (1 + 1e-9) / s), 1e-8);
}

TEST(Smeared, ExponentialAndItsDerivativeInTheStableForm) {
  for (double sigma : {0.8, 1.8, 3.5})
    for (double mu : {0.05, 0.4, 1.5})
      for (double r : {0.1, 0.5, 2.0, 6.0, 15.0}) {
        const double ref = smear_by_quadrature(
            [&](double x) { return std::exp(-mu * x); }, sigma, r);
        EXPECT_NEAR(smeared_exp(sigma, mu, r), ref, 1e-10)
            << sigma << " " << mu << " " << r;
        const double h = 1e-5;
        const double d =
            (smeared_exp(sigma, mu, r + h) - smeared_exp(sigma, mu, r - h)) /
            (2 * h);
        EXPECT_NEAR(smeared_exp_prime(sigma, mu, r), d, 2e-8)
            << sigma << " " << mu << " " << r;
      }
  // checks of notes 6.5.5: sigma -> infinity recovers e^{-mu r}; mu -> 0
  // recovers 1; large r: e^{mu^2/4 s^2} e^{-mu r}
  EXPECT_NEAR(smeared_exp(600.0, 0.4, 1.0), std::exp(-0.4),
              1e-6);  // corrections ~ mu^2/(4 sigma^2)
  EXPECT_NEAR(smeared_exp(1.5, 1e-9, 2.0), 1.0, 1e-8);
  {
    const double sg = 2.0, m = 0.6, r = 20.0,
                 a = m / (2 * sg * sg);  // tail: e^{beta^2} (1 - a/r) e^{-mu r}
    EXPECT_NEAR(smeared_exp(sg, m, r) / std::exp(-m * r),
                std::exp(m * m / (4 * sg * sg)) * (1 - a / r), 1e-9);
  }
  // no overflow for large sigma r and large mu/sigma (the reason for the erfcx
  // form)
  EXPECT_TRUE(std::isfinite(smeared_exp(5.0, 3.0, 40.0)));
  EXPECT_TRUE(std::isfinite(smeared_exp_prime(5.0, 3.0, 40.0)));
}
