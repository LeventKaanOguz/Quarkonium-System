#include "me/radial.hpp"

#include <algorithm>
#include <cmath>
#include <complex>
#include <numbers>
#include <stdexcept>
#include <vector>

#include "math/quadrature.hpp"
#include "math/special.hpp"
#include "me/analytic.hpp"

namespace qgem {
namespace {

using Cx = std::complex<double>;

template <class Elem>
Mat fill(const Basis& bra, const Basis& ket, Elem&& f) {
  Mat M(bra.size(), ket.size());
  for (int i = 0; i < bra.size(); ++i)
    for (int j = 0; j < ket.size(); ++j) M(i, j) = f(i, j);
  return M;
}

double lbar(const Basis& bra, const Basis& ket) {
  return 0.5 * (bra.L + ket.L);
}

Cx exponent(const Basis& b, int i) {
  return Cx(b.nu(i), b.kind[i] == Fn::Real ? 0.0 : b.omega(i) * b.nu(i));
}

// Matrix of a real operator in a basis that may contain complex-range
// functions, from the unnormalized bilinear form B(z_i, z_j) = Int r^2 [r^{l_i}
// e^{-z_i r^2}] O [r^{l_j} e^{-z_j r^2}] dr (no complex conjugation), by the
// rules of notes (7.33):
//   <c|O|c> = 1/2 Re[B(zi,zj) + B(zi,zj*)],   <s|O|s> = 1/2 Re[B(zi,zj*) -
//   B(zi,zj)], <c|O|s> = -1/2 Im[B(zi,zj) - B(zi,zj*)];
// real Gaussians are the omega = 0 case of "c".
template <class B>
Mat bilinear_matrix(const Basis& bra, const Basis& ket, B&& bf) {
  return fill(bra, ket, [&](int i, int j) {
    const Cx zi = exponent(bra, i), zj = exponent(ket, j);
    const bool si = bra.kind[i] == Fn::Sin, sj = ket.kind[j] == Fn::Sin;
    double v;
    if (!si && !sj) {
      v = 0.5 * (bf(zi, zj) + bf(zi, std::conj(zj))).real();
    } else if (si && sj) {
      v = 0.5 * (bf(zi, std::conj(zj)) - bf(zi, zj)).real();
    } else if (!si && sj) {
      v = -0.5 * (bf(zi, zj) - bf(zi, std::conj(zj))).imag();
    } else {  // <s|O|c>: use the symmetry of the bilinear form
      v = -0.5 * (bf(zj, zi) - bf(zj, std::conj(zi))).imag();
    }
    return std::exp(0.5 * (bra.lognorm2(i) + ket.lognorm2(j))) * v;
  });
}

Cx master(double s, Cx a) {
  return std::exp(std::lgamma(s) - std::log(2.0)) * std::pow(a, -s);
}

// Composite Gauss-Legendre of N' f_i(r) V(r) f_j(r) r^2 evaluated with the
// actual real/cos/sin basis functions.
Mat local_by_values(const Basis& bra, const Basis& ket,
                    const std::function<double(double)>& V) {
  const double nu_min = std::min(bra.nu.minCoeff(), ket.nu.minCoeff());
  const double rc = std::sqrt(40.0 / (2 * nu_min));
  const int panels = 48;
  const QuadRule& q = gauss_legendre01(24);
  std::vector<double> r, w;
  for (int p = 0; p < panels; ++p) {
    // quadratic spacing: dense near the origin where the narrow functions live
    const double a = rc * std::pow(static_cast<double>(p) / panels, 2),
                 b = rc * std::pow(static_cast<double>(p + 1) / panels, 2);
    for (int k = 0; k < q.x.size(); ++k) {
      r.push_back(a + (b - a) * q.x(k));
      w.push_back(q.w(k) * (b - a));
    }
  }
  const int n = static_cast<int>(r.size());
  auto tab = [&](const Basis& B) {
    Mat F(B.size(), n);
    for (int i = 0; i < B.size(); ++i)
      for (int m = 0; m < n; ++m) {
        const double g = std::exp(-B.nu(i) * r[m] * r[m]),
                     ph = B.omega(i) * B.nu(i) * r[m] * r[m];
        const double f =
            B.kind[i] == Fn::Real
                ? g
                : (B.kind[i] == Fn::Cos ? g * std::cos(ph) : g * std::sin(ph));
        F(i, m) = B.norm(i) * std::pow(r[m], B.L) * f;
      }
    return F;
  };
  const Mat Fb = tab(bra), Fk = tab(ket);
  Vec wv(n);
  for (int m = 0; m < n; ++m) wv(m) = w[m] * r[m] * r[m] * V(r[m]);
  return Fb * wv.asDiagonal() * Fk.transpose();
}

}  // namespace

Mat overlap(const Basis& b) {
  if (!b.is_real())
    return bilinear_matrix(
        b, b, [&](Cx zi, Cx zj) { return master(b.L + 1.5, zi + zj); });
  const double e = b.L + 1.5;
  return fill(b, b, [&](int i, int j) {
    if (i == j) return 1.0;
    return std::pow(2.0 * std::sqrt(b.nu(i) * b.nu(j)) / (b.nu(i) + b.nu(j)),
                    e);
  });
}

Mat p2(const Basis& b) {
  if (!b.is_real()) {
    // p^2 u_j = [2 z_j (2l+3) - 4 z_j^2 r^2] u_j  (Lemma 7.6 holds for complex
    // z_j)
    return bilinear_matrix(b, b, [&](Cx zi, Cx zj) {
      const Cx a = zi + zj;
      return 2.0 * zj * (2.0 * b.L + 3.0) * master(b.L + 1.5, a) -
             4.0 * zj * zj * master(b.L + 2.5, a);
    });
  }
  const Mat S = overlap(b);
  return fill(b, b, [&](int i, int j) {
    return 2.0 * (2 * b.L + 3) * b.nu(i) * b.nu(j) / (b.nu(i) + b.nu(j)) *
           S(i, j);
  });
}

Mat kinetic_nr(const Basis& b, double mr) { return p2(b) / (2.0 * mr); }

Mat r_power(const Basis& bra, const Basis& ket, double k) {
  if (!(bra.L + ket.L + 3 + k > 0.0))
    throw std::invalid_argument("r_power: divergent moment");
  const double lb = lbar(bra, ket), s = lb + 0.5 * (k + 3.0);
  if (!bra.is_real() || !ket.is_real())
    return bilinear_matrix(bra, ket,
                           [&](Cx zi, Cx zj) { return master(s, zi + zj); });
  const double g = std::lgamma(s) - std::log(2.0);
  return fill(bra, ket, [&](int i, int j) {
    const double a = bra.nu(i) + ket.nu(j);
    return std::exp(0.5 * (bra.lognorm2(i) + ket.lognorm2(j)) + g -
                    s * std::log(a));
  });
}

Mat gaussian(const Basis& bra, const Basis& ket, double beta, double k) {
  if (beta < 0) throw std::invalid_argument("gaussian: beta >= 0");
  if (!(bra.L + ket.L + 3 + k > 0.0))
    throw std::invalid_argument("gaussian: divergent moment");
  const double lb = lbar(bra, ket), s = lb + 0.5 * (k + 3.0);
  if (!bra.is_real() || !ket.is_real())
    return bilinear_matrix(
        bra, ket, [&](Cx zi, Cx zj) { return master(s, zi + zj + beta); });
  const double g = std::lgamma(s) - std::log(2.0);
  return fill(bra, ket, [&](int i, int j) {
    const double a = bra.nu(i) + ket.nu(j) + beta;
    return std::exp(0.5 * (bra.lognorm2(i) + ket.lognorm2(j)) + g -
                    s * std::log(a));
  });
}

Mat erf_over_r(const Basis& bra, const Basis& ket, double tau) {
  if (!(tau > 0)) throw std::invalid_argument("erf_over_r: tau > 0");
  if (!bra.is_real() || !ket.is_real())
    return local_by_values(bra, ket,
                           [&](double r) { return std::erf(tau * r) / r; });
  // <erf(tau r)/r> = N'N Gamma(l_bar+3/2)/sqrt(pi) J_{l_bar+3/2}(a, tau):
  // exact, from the elementary integral J_s = Int_0^tau (a+t^2)^{-s} dt (upward
  // recursion from asinh / atan).
  const double lb = lbar(bra, ket),
               lg = std::lgamma(lb + 1.5) - 0.5 * std::log(std::numbers::pi);
  const int two_s = static_cast<int>(std::lround(2 * lb + 3));
  return fill(bra, ket, [&](int i, int j) {
    const double a = bra.nu(i) + ket.nu(j);
    return std::exp(0.5 * (bra.lognorm2(i) + ket.lognorm2(j)) + lg) *
           inverse_quadratic_integral(two_s, a, tau);
  });
}

Mat exp_power(const Basis& bra, const Basis& ket, double mu, int k) {
  const int idx = bra.L + ket.L + 2 + k;
  if (idx < 0) throw std::invalid_argument("exp_power: divergent moment");
  if (!bra.is_real() || !ket.is_real())
    return local_by_values(
        bra, ket, [&](double r) { return std::pow(r, k) * std::exp(-mu * r); });
  std::vector<double> Q(idx + 1);
  return fill(bra, ket, [&](int i, int j) {
    screened_moments(idx, bra.nu(i) + ket.nu(j), mu, Q);
    return std::exp(0.5 * (bra.lognorm2(i) + ket.lognorm2(j))) * Q[idx];
  });
}

Mat local(const Basis& bra, const Basis& ket,
          const std::function<double(double)>& V, int nq) {
  if (!bra.is_real() || !ket.is_real()) return local_by_values(bra, ket, V);
  const QuadRule& q = gauss_legendre01(nq);
  const int p = bra.L + ket.L + 2;
  const bool same = bra.L == ket.L && bra.size() == ket.size() &&
                    (bra.nu - ket.nu).norm() == 0.0;
  Mat M(bra.size(), ket.size());
  // the nodes depend on the pair only through r_c = sqrt(40/a); symmetric bases
  // fill one triangle
  for (int i = 0; i < bra.size(); ++i)
    for (int j = same ? i : 0; j < ket.size(); ++j) {
      const double a = bra.nu(i) + ket.nu(j), rc = std::sqrt(40.0 / a);
      double sum = 0;
      for (int m = 0; m < nq; ++m) {
        const double r = rc * q.x(m);
        sum += q.w(m) * std::pow(r, p) * std::exp(-a * r * r) * V(r);
      }
      M(i, j) = std::exp(0.5 * (bra.lognorm2(i) + ket.lognorm2(j))) * rc * sum;
      if (same) M(j, i) = M(i, j);
    }
  return M;
}

}  // namespace qgem
