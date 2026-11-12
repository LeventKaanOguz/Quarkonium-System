#include "me/analytic.hpp"

#include <array>
#include <cmath>
#include <numbers>
#include <stdexcept>
#include <vector>

#include "math/special.hpp"

namespace qgem {
namespace {
constexpr double kSqrtPi = 1.7724538509055160273;

double lbar(const Basis& a, const Basis& b) { return 0.5 * (a.L + b.L); }
double nn(const Basis& bra, const Basis& ket, int i, int j) {
  return std::exp(0.5 * (bra.lognorm2(i) + ket.lognorm2(j)));
}

// Int_0^inf r^{2m-1} e^{-a r^2} erf(tau r) dr = Gamma(m+1/2)/sqrt(pi)
// J_{m+1/2}(a, tau), m >= 0 (integer or half-integer)
double erf_moment(double m, double a, double tau) {
  return std::exp(std::lgamma(m + 0.5)) / kSqrtPi *
         inverse_quadratic_integral(static_cast<int>(std::lround(2 * m + 1)), a,
                                    tau);
}
// Int_0^inf r^{2m} e^{-c r^2} dr
double gauss_moment(double m, double c) {
  return std::exp(std::lgamma(m + 0.5) - std::log(2.0)) *
         std::pow(c, -(m + 0.5));
}

template <class F>
Mat fill(const Basis& bra, const Basis& ket, F&& f) {
  if (!bra.is_real() || !ket.is_real())
    throw std::logic_error("analytic matrix elements need real-range bases");
  Mat M(bra.size(), ket.size());
  for (int i = 0; i < bra.size(); ++i)
    for (int j = 0; j < ket.size(); ++j)
      M(i, j) = nn(bra, ket, i, j) * f(bra.nu(i) + ket.nu(j));
  return M;
}

// ---- truncated Taylor series in eps (a = a0 + eps)
// -------------------------------------------------------------------
struct Jet {
  std::vector<double> c;
  explicit Jet(int n) : c(n + 1, 0.0) {}
};
Jet operator*(const Jet& x, const Jet& y) {
  Jet r(static_cast<int>(x.c.size()) - 1);
  for (std::size_t i = 0; i < x.c.size(); ++i)
    for (std::size_t j = 0; i + j < x.c.size(); ++j)
      r.c[i + j] += x.c[i] * y.c[j];
  return r;
}
// (x0 + eps)^p
Jet jet_pow(double x0, double p, int n) {
  Jet r(n);
  double coef = 1.0;
  for (int m = 0; m <= n; ++m) {
    r.c[m] = coef * std::pow(x0, p - m);
    coef *= (p - m) / (m + 1);
  }
  return r;
}
// l-th derivative with respect to a of sigma^{2q} (a+s2)^{-q'}... generic:
// d^l/da^l of  H(a) * Phi(a'(a)), H(a) = (s2/(a+s2))^{p}, a' = a s2/(a+s2),
// Phi^{(m)}(a') = coefs[m] (m-th derivative at a'0).
double transfer_derivative(int l, double a, double s2, double p,
                           const std::vector<double>& phi_derivs) {
  const Jet hh = [&] {
    Jet j = jet_pow(a + s2, -p, l);
    const double f = std::pow(s2, p);
    for (double& x : j.c) x *= f;
    return j;
  }();
  // a'(eps) = s2 (1 - s2/(a+s2+eps)) = s2 - s2^2 (a+s2+eps)^{-1}
  Jet ap = jet_pow(a + s2, -1.0, l);
  for (double& x : ap.c) x *= -s2 * s2;
  ap.c[0] += s2;
  const double ap0 = ap.c[0];
  Jet delta = ap;
  delta.c[0] = 0.0;
  // Phi(a'(eps)) = sum_m Phi^{(m)}(a'0)/m! delta^m
  Jet phi(l), pw(l);
  pw.c[0] = 1.0;
  double fact = 1.0;
  for (int m = 0; m <= l; ++m) {
    if (m > 0) fact *= m;
    for (int k = 0; k <= l; ++k) phi.c[k] += phi_derivs[m] / fact * pw.c[k];
    pw = pw * delta;
  }
  (void)ap0;
  const Jet prod = hh * phi;
  double lfact = 1.0;
  for (int m = 2; m <= l; ++m) lfact *= m;
  return prod.c[l] * lfact * ((l % 2) ? -1.0 : 1.0);  // (-d/da)^l
}
}  // namespace

double inverse_quadratic_integral(int two_s, double a, double tau) {
  if (two_s < 1)
    throw std::invalid_argument("inverse_quadratic_integral: s >= 1/2");
  if (!(a > 0))
    throw std::invalid_argument("inverse_quadratic_integral: a > 0");
  const double sa = std::sqrt(a);
  // J_{s+1} = [tau (a+tau^2)^{-s} + (2s-1) J_s] / (2 s a)
  int k = (two_s % 2 == 1) ? 1 : 2;  // base: 2s = 1 (asinh) or 2s = 2 (atan)
  double J = (k == 1) ? std::asinh(tau / sa) : std::atan(tau / sa) / sa;
  while (k < two_s) {
    const double s = 0.5 * k;
    J = (tau * std::pow(a + tau * tau, -s) + (2 * s - 1) * J) / (2 * s * a);
    k += 2;
  }
  return J;
}

Mat gprime_over_r_matrix(const Basis& bra, const Basis& ket, double tau) {
  const double lb = lbar(bra, ket);
  // N'N [ (tau/sqrt(pi)) Gamma(l+1/2)/(a+tau^2)^{l+1/2} -
  // (Gamma(l+1/2)/sqrt(pi)) J_{l+1/2} ]
  return fill(bra, ket, [&](double a) {
    const double g = std::exp(std::lgamma(lb + 0.5)) / kSqrtPi;
    return g * (tau * std::pow(a + tau * tau, -(lb + 0.5)) -
                inverse_quadratic_integral(
                    static_cast<int>(std::lround(2 * lb + 1)), a, tau));
  });
}

Mat tensor_matrix(const Basis& bra, const Basis& ket, double tau) {
  const double lb = lbar(bra, ket);
  const Mat g = gprime_over_r_matrix(bra, ket, tau);
  Mat G = Mat(bra.size(), ket.size());
  for (int i = 0; i < bra.size(); ++i)
    for (int j = 0; j < ket.size(); ++j) {
      const double a = bra.nu(i) + ket.nu(j);
      G(i, j) = nn(bra, ket, i, j) * gauss_moment(lb + 1.0, a + tau * tau);
    }
  return 3.0 * g + 4.0 * std::pow(tau, 3) / kSqrtPi * G;
}

Mat smeared_r_matrix(const Basis& bra, const Basis& ket, double s) {
  const double lb = lbar(bra, ket);
  return fill(bra, ket, [&](double a) {
    return gauss_moment(lb + 1.0, a + s * s) / (kSqrtPi * s) +
           erf_moment(lb + 2.0, a, s) +
           erf_moment(lb + 1.0, a, s) / (2 * s * s);
  });
}

Mat smeared_r_prime_over_r_matrix(const Basis& bra, const Basis& ket,
                                  double s) {
  const double lb = lbar(bra, ket);
  return fill(bra, ket, [&](double a) {
    return erf_moment(lb + 1.0, a, s) - erf_moment(lb, a, s) / (2 * s * s) +
           gauss_moment(lb, a + s * s) / (kSqrtPi * s);
  });
}

Mat smeared_exp_matrix(const Basis& b, double s, double mu) {
  const double s2 = s * s;
  const int l = b.L;
  Mat M(b.size(), b.size());
  std::vector<double> Q(2 + 2 * l + 1), dq(l + 1);
  for (int i = 0; i < b.size(); ++i)
    for (int j = 0; j < b.size(); ++j) {
      const double a = b.nu(i) + b.nu(j), ap = a * s2 / (a + s2);
      screened_moments(2 + 2 * l, ap, mu, Q);
      for (int m = 0; m <= l; ++m)
        dq[m] = ((m % 2) ? -1.0 : 1.0) *
                Q[2 + 2 * m];  // d^m/da'^m Q_2 = (-1)^m Q_{2+2m}
      M(i, j) = nn(b, b, i, j) * transfer_derivative(l, a, s2, 1.5, dq);
    }
  return M;
}

Mat smeared_exp_prime_over_r_matrix(const Basis& b, double s, double mu) {
  const double s2 = s * s;
  const int l = b.L;
  if (l < 1)
    return Mat::Zero(
        b.size(),
        b.size());  // the spin-orbit operator vanishes for l = 0 (not needed)
  Mat M(b.size(), b.size());
  std::vector<double> Q(3 + 2 * l + 1), dq(l);
  for (int i = 0; i < b.size(); ++i)
    for (int j = 0; j < b.size(); ++j) {
      const double a = b.nu(i) + b.nu(j), ap = a * s2 / (a + s2);
      screened_moments(3 + 2 * l, ap, mu, Q);
      // J_f(a') = Int r^3 f'(r) e^{-a' r^2} dr = -mu Q_3(a', mu);  d^m/da'^m J
      // = -mu (-1)^m Q_{3+2m}
      for (int m = 0; m < l; ++m)
        dq[m] = -mu * ((m % 2) ? -1.0 : 1.0) * Q[3 + 2 * m];
      M(i, j) = nn(b, b, i, j) * transfer_derivative(l - 1, a, s2, 2.5, dq);
    }
  return M;
}

}  // namespace qgem
