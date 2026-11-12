#include "me/momentum.hpp"

#include <cmath>
#include <numbers>
#include <stdexcept>
#include <vector>

#include "math/quadrature.hpp"
#include "math/special.hpp"
#include "me/radial.hpp"

namespace qgem {

Mat momentum_diagonal(const Basis& b, const std::function<double(double)>& f,
                      int nq) {
  if (!b.is_real())
    throw std::logic_error(
        "momentum_diagonal: complex-range bases are not supported for "
        "relativistic kinetic terms");
  const QuadRule& q = gauss_laguerre(nq, b.L + 0.5);
  const Mat S = overlap(b);
  const double g = std::tgamma(b.L + 1.5);
  Mat M(b.size(), b.size());
  for (int i = 0; i < b.size(); ++i)
    for (int j = i; j < b.size(); ++j) {  // symmetric in i <-> j
      const double beta = (b.nu(i) + b.nu(j)) / (4.0 * b.nu(i) * b.nu(j));
      double sum = 0;
      for (int k = 0; k < nq; ++k) sum += q.w(k) * f(std::sqrt(q.x(k) / beta));
      M(i, j) = M(j, i) = S(i, j) / g * sum;
    }
  return M;
}

Mat sandwich_closure(const Mat& F, const Mat& W, const Mat& G, const Mat& Xb,
                     const Mat& Xk) {
  return (Xb.transpose() * F * Xb) * (Xb.transpose() * W * Xk) *
         (Xk.transpose() * G * Xk);
}

namespace {

// chi_n(r_q) for all n on the nodes: sqrt(2/pi) Int_0^inf p^2 j_l(p r) f(p)
// R~_n(p) dp R~_nl(p) = N_nl (2 nu)^{-l-3/2} p^l exp(-p^2/(4 nu)) (notes 7.25)
Mat hankel_table(const Basis& B, const std::function<double(double)>& f,
                 const Vec& r) {
  const int l = B.L;
  double nu_max = B.nu.maxCoeff();
  const double pmax = std::sqrt(4.0 * nu_max * 45.0);
  // composite Gauss-Legendre in p, panel width resolving the fastest
  // oscillation sin(p r_max)
  const double rmax = r.maxCoeff();
  const double h = std::min(0.5, 1.0 / rmax * 3.0);
  const int panels = static_cast<int>(std::ceil(pmax / h));
  const QuadRule& g = gauss_legendre01(10);
  std::vector<double> p, wp;
  p.reserve(panels * g.x.size());
  for (int a = 0; a < panels; ++a)
    for (int k = 0; k < g.x.size(); ++k) {
      p.push_back((a + g.x(k)) * h);
      wp.push_back(g.w(k) * h);
    }
  const int np = static_cast<int>(p.size());
  Mat A(B.size(), np);  // A_{n,k} = w p^2 f(p) R~_n(p)
  for (int n = 0; n < B.size(); ++n) {
    const double pref = B.norm(n) * std::pow(2 * B.nu(n), -l - 1.5);
    for (int k = 0; k < np; ++k)
      A(n, k) = wp[k] * p[k] * p[k] * f(p[k]) * pref * std::pow(p[k], l) *
                std::exp(-p[k] * p[k] / (4 * B.nu(n)));
  }
  Mat chi(r.size(), B.size());
  const int chunk = 256;
  Mat J(r.size(), chunk);
  for (int k0 = 0; k0 < np; k0 += chunk) {
    const int kn = std::min(chunk, np - k0);
    for (int q = 0; q < r.size(); ++q)
      for (int k = 0; k < kn; ++k) J(q, k) = sph_bessel_j(l, p[k0 + k] * r(q));
    if (k0 == 0) chi.setZero();
    chi += J.leftCols(kn) * A.middleCols(k0, kn).transpose();
  }
  return chi * std::sqrt(2.0 / std::numbers::pi);
}

}  // namespace

Mat sandwich_two_step(const Basis& bra, const Basis& ket,
                      const std::function<double(double)>& f,
                      const std::function<double(double)>& W, double rmax) {
  // radial nodes: geometric panels toward the origin (narrow functions) up to
  // rmax, 16-point Gauss-Legendre each
  const QuadRule& g = gauss_legendre01(16);
  std::vector<double> r, wr;
  std::vector<double> edges = {0.0};
  for (double e = 0.004; e < rmax; e *= 1.6) edges.push_back(e);
  edges.push_back(rmax);
  for (std::size_t a = 0; a + 1 < edges.size(); ++a)
    for (int k = 0; k < g.x.size(); ++k) {
      r.push_back(edges[a] + g.x(k) * (edges[a + 1] - edges[a]));
      wr.push_back(g.w(k) * (edges[a + 1] - edges[a]));
    }
  Vec rv = Eigen::Map<Vec>(r.data(), r.size());
  const Mat chib = hankel_table(bra, f, rv),
            chik = (ket.L == bra.L && ket.nu.size() == bra.nu.size() &&
                    (ket.nu - bra.nu).norm() == 0)
                       ? chib
                       : hankel_table(ket, f, rv);
  Vec wq(r.size());
  for (std::size_t q = 0; q < r.size(); ++q)
    wq(q) = wr[q] * r[q] * r[q] * W(r[q]);
  return chib.transpose() * wq.asDiagonal() * chik;
}

}  // namespace qgem
