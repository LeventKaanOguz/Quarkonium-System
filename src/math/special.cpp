#include "math/special.hpp"

#include <cmath>
#include <numbers>
#include <stdexcept>
#include <vector>

namespace qgem {

double erfcx(double x) {
  if (x < 0.0)
    return 2.0 * std::exp(x * x) - erfcx(-x);  // overflows only for x < -26.6
  if (x < 26.0) return std::exp(x * x) * std::erfc(x);
  // Asymptotic series e^{x^2} erfc(x) ~ 1/(x sqrt(pi)) * sum (-1)^n (2n-1)!! /
  // (2x^2)^n.
  const double y = 1.0 / (2.0 * x * x);
  const double s =
      1.0 -
      y * (1.0 -
           3.0 * y *
               (1.0 -
                5.0 * y *
                    (1.0 - 7.0 * y * (1.0 - 9.0 * y * (1.0 - 11.0 * y)))));
  return s / (x * std::sqrt(std::numbers::pi));
}

double master_integral(double p, double a) {
  if (!(a > 0.0))
    throw std::invalid_argument("master_integral: a must be positive");
  if (!(p > -1.0)) throw std::invalid_argument("master_integral: need p > -1");
  const double s = 0.5 * (p + 1.0);
  return std::exp(std::lgamma(s) - std::log(2.0) - s * std::log(a));
}

void screened_moments(int kmax, double a, double mu, std::span<double> Q) {
  if (kmax < 0 || Q.size() < static_cast<std::size_t>(kmax) + 1)
    throw std::invalid_argument("screened_moments: output too small");
  if (!(a > 0.0) || mu < 0.0)
    throw std::invalid_argument("screened_moments: need a > 0, mu >= 0");

  const double sqa = std::sqrt(a);
  const double x = mu / (2.0 * sqa);
  // q_k = a^{(k+1)/2} Q_k depends on x only.
  std::vector<double> q(kmax + 2);
  q[0] = 0.5 * std::sqrt(std::numbers::pi) * erfcx(x);
  q[1] = 0.5 * (1.0 - 2.0 * x * q[0]);
  if (x < 0.5) {
    for (int k = 0; k + 2 <= kmax + 1; ++k)
      q[k + 2] = 0.5 * (k + 1) * q[k] - x * q[k + 1];
  } else {
    // q_k is the minimal solution of the three-term recursion: recur rho_k =
    // q_{k+1}/q_k downward, rho_k = (k+1) / (2x + 2 rho_{k+1}), which only adds
    // positive numbers.
    const double t = 20.0 / x + 2.0 * kmax;
    const int K = static_cast<int>(0.5 * t * t) + 10;
    std::vector<double> rho(K + 1, 0.0);
    for (int k = K - 1; k >= 0; --k)
      rho[k] = (k + 1.0) / (2.0 * x + 2.0 * rho[k + 1]);
    for (int k = 1; k <= kmax; ++k) q[k] = rho[k - 1] * q[k - 1];
  }
  double ak = sqa;  // a^{(k+1)/2}
  for (int k = 0; k <= kmax; ++k) {
    Q[k] = q[k] / ak;
    ak *= sqa;
  }
}

double sph_bessel_j(int l, double x) {
  if (l < 0) throw std::invalid_argument("sph_bessel_j: l >= 0");
  if (x == 0.0) return l == 0 ? 1.0 : 0.0;
  if (x < static_cast<double>(l) + 1.5 && x < 4.0) {
    // j_l(x) = x^l / (2l+1)!! * sum_k (-x^2/2)^k / (k!
    // (2l+3)(2l+5)...(2l+2k+1))
    double dfact = 1.0;
    for (int i = 1; i <= l; ++i) dfact *= (2 * i + 1);
    double term = 1.0, sum = 1.0;
    for (int k = 1; k < 60; ++k) {
      term *= -0.5 * x * x / (k * (2.0 * l + 2 * k + 1));
      sum += term;
      if (std::fabs(term) < 1e-17 * std::fabs(sum)) break;
    }
    return std::pow(x, l) / dfact * sum;
  }
  const double s = std::sin(x), c = std::cos(x);
  double j0 = s / x, j1 = s / (x * x) - c / x;
  if (l == 0) return j0;
  for (int n = 1; n < l; ++n) {
    const double j2 = (2.0 * n + 1.0) / x * j1 - j0;
    j0 = j1;
    j1 = j2;
  }
  return j1;
}

}  // namespace qgem
