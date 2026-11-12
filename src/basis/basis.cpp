#include "basis/basis.hpp"

#include <cmath>
#include <complex>
#include <stdexcept>

namespace qgem {
namespace {
double log_norm2(int l, double nu) {
  return std::log(2.0) + (l + 1.5) * std::log(2.0 * nu) - std::lgamma(l + 1.5);
}
}  // namespace

double Basis::norm(int i) const { return std::exp(0.5 * lognorm2(i)); }

Basis make_basis(int L, const BasisSpec& spec) {
  if (L < 0) throw std::invalid_argument("make_basis: L >= 0");
  if (spec.n < 1 || !(spec.r1 > 0) || !(spec.rmax >= spec.r1))
    throw std::invalid_argument("make_basis: need n >= 1, 0 < r1 <= rmax");
  if (spec.nc < 0 ||
      (spec.nc > 0 && !(spec.rho1 > 0 && spec.rhomax >= spec.rho1)))
    throw std::invalid_argument(
        "make_basis: need 0 < rho1 <= rhomax for the complex-range set");
  Basis b;
  b.L = L;
  b.spec = spec;
  const int N = spec.n + 2 * spec.nc;
  b.nu.resize(N);
  b.omega = Vec::Zero(N);
  b.kind.assign(N, Fn::Real);
  b.lognorm2.resize(N);
  const double q =
      spec.n > 1 ? std::pow(spec.rmax / spec.r1, 1.0 / (spec.n - 1)) : 1.0;
  for (int k = 0; k < spec.n; ++k) {
    const double rn = spec.r1 * std::pow(q, k);
    b.nu(k) = 1.0 / (rn * rn);
    b.lognorm2(k) = log_norm2(L, b.nu(k));
  }
  if (spec.nc > 0) {
    const double qc =
        spec.nc > 1 ? std::pow(spec.rhomax / spec.rho1, 1.0 / (spec.nc - 1))
                    : 1.0;
    for (int k = 0; k < spec.nc; ++k) {
      const double rho = spec.rho1 * std::pow(qc, k), eta = 1.0 / (rho * rho);
      for (int part = 0; part < 2; ++part) {  // cos then sin
        const int i = spec.n + 2 * k + part;
        b.nu(i) = eta;
        b.omega(i) = spec.omega;
        b.kind[i] = part == 0 ? Fn::Cos : Fn::Sin;
        // N^-2 = <f|f> from the bilinear rules (7.33) with the closed form
        // B0(a) = Gamma(l+3/2) / (2 a^{l+3/2}) at complex exponent sums
        const std::complex<double> z(eta, eta * spec.omega);
        const double pre = std::exp(std::lgamma(L + 1.5) - std::log(2.0));
        const auto B0 = [&](std::complex<double> a) {
          return pre * std::pow(a, -(L + 1.5));
        };
        const double bb = B0(z + z).real(), bc = B0(z + std::conj(z)).real();
        const double S = 0.5 * (part == 0 ? (bb + bc) : (bc - bb));
        b.lognorm2(i) = -std::log(S);
      }
    }
  }
  return b;
}

double radial_value(const Basis& b, const Vec& c, double r) {
  double sum = 0;
  for (int i = 0; i < b.size(); ++i) {
    const double g = std::exp(-b.nu(i) * r * r),
                 ph = b.omega(i) * b.nu(i) * r * r;
    const double f =
        b.kind[i] == Fn::Real
            ? g
            : (b.kind[i] == Fn::Cos ? g * std::cos(ph) : g * std::sin(ph));
    sum += c(i) * b.norm(i) * f;
  }
  return sum * std::pow(r, b.L);
}

Vec origin_coefficients(const Basis& b) {
  // coefficient of r^L at the origin: real and cos functions give N, sin
  // functions vanish faster
  Vec v(b.size());
  for (int i = 0; i < b.size(); ++i)
    v(i) = (b.kind[i] == Fn::Sin) ? 0.0 : b.norm(i);
  return v;
}

}  // namespace qgem
