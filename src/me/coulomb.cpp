#include "me/coulomb.hpp"

#include <cmath>
#include <numbers>

namespace qgem {
namespace {
constexpr double kSqrtPi = 1.7724538509055160273;

// sum_k (-1)^k z^{2k} / (k! (2k+3)) and sum_k (-1)^k z^{2k} 2k/(k! (2k+3)), z
// < 1.
void series(double z, double& h_sum, double& t_sum) {
  const double w = -z * z;
  double term = 1.0, hs = 0.0, ts = 0.0;  // term = (-z^2)^k / k!
  for (int k = 0; k < 40; ++k) {
    hs += term / (2.0 * k + 3.0);
    ts += term * 2.0 * k / (2.0 * k + 3.0);
    term *= w / (k + 1);
  }
  h_sum = hs;
  t_sum = ts;
}
}  // namespace

double coulomb_gp_over_r(double tau, double r) {
  const double z = tau * r;
  const double t3 = tau * tau * tau;
  if (z < 1.0) {
    double hs, ts;
    series(z, hs, ts);
    return -t3 * 4.0 / kSqrtPi * hs;
  }
  const double h =
      2.0 / kSqrtPi * std::exp(-z * z) / (z * z) - std::erf(z) / (z * z * z);
  return t3 * h;
}

double coulomb_tensor(double tau, double r) {
  const double z = tau * r;
  const double t3 = tau * tau * tau;
  if (z < 1.0) {
    double hs, ts;
    series(z, hs, ts);
    return t3 * 4.0 / kSqrtPi * ts;
  }
  const double h =
      2.0 / kSqrtPi * std::exp(-z * z) / (z * z) - std::erf(z) / (z * z * z);
  return t3 * (3.0 * h + 4.0 / kSqrtPi * std::exp(-z * z));
}

}  // namespace qgem
