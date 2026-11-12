#include "me/smeared.hpp"

#include <algorithm>
#include <cmath>
#include <numbers>

#include "math/special.hpp"

namespace qgem {
namespace {
constexpr double kSqrtPi = 1.7724538509055160273;
constexpr double kRFloor =
    1e-6;  // smeared functions are regular at 0: clamp to avoid 0/0
}  // namespace

double smeared_r(double s, double r) {
  r = std::max(r, kRFloor);
  const double z = s * r;
  return std::exp(-z * z) / (kSqrtPi * s) +
         (r + 1.0 / (2.0 * s * s * r)) * std::erf(z);
}

double smeared_r_prime_over_r(double s, double r) {
  r = std::max(r, kRFloor);
  const double z = s * r;
  if (z < 1.0) {
    // r~' = (4/sqrt(pi)) sum_m (-1)^m z^{2m+1} / (m! (2m+1)(2m+3))   =>  r~'/r
    // = s (4/sqrt(pi)) sum_m (-1)^m z^{2m}/(...)
    double term = 1.0, sum = 0.0;
    for (int m = 0; m < 40; ++m) {
      sum += term / ((2.0 * m + 1.0) * (2.0 * m + 3.0));
      term *= -z * z / (m + 1);
    }
    return s * 4.0 / kSqrtPi * sum;
  }
  const double prime = (1.0 - 1.0 / (2.0 * z * z)) * std::erf(z) +
                       std::exp(-z * z) / (kSqrtPi * z);
  return prime / r;
}

double smeared_exp(double s, double mu, double r) {
  r = std::max(r, 1e-5);
  const double a = mu / (2 * s * s), beta = mu / (2 * s), z = s * r;
  const double first =
      (r - a) * std::exp(beta * beta - mu * r) * std::erfc(beta - z);
  const double second = (r + a) * std::exp(-z * z) * erfcx(beta + z);
  return (first + second) / (2.0 * r);
}

double smeared_exp_prime(double s, double mu, double r) {
  r = std::max(r, 1e-5);
  const double a = mu / (2 * s * s), beta = mu / (2 * s), z = s * r;
  const double e1 = std::exp(beta * beta - mu * r) * std::erfc(beta - z);
  const double e2 = std::exp(-z * z) * erfcx(beta + z);
  const double F = ((r - a) * e1 + (r + a) * e2) / (2.0 * r);
  return -F / r - mu / (kSqrtPi * s * r) * std::exp(-z * z) +
         ((1.0 - mu * (r - a)) * e1 + (1.0 + mu * (r + a)) * e2) / (2.0 * r);
}

}  // namespace qgem
