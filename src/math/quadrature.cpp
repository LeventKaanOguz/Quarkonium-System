#include "math/quadrature.hpp"

#include <cmath>
#include <map>
#include <mutex>
#include <stdexcept>
#include <utility>

namespace qgem {
namespace {

QuadRule golub_welsch(const Vec& diag, const Vec& offdiag, double mu0) {
  const int n = static_cast<int>(diag.size());
  Mat T = Mat::Zero(n, n);
  for (int k = 0; k < n; ++k) T(k, k) = diag(k);
  for (int k = 1; k < n; ++k) T(k, k - 1) = T(k - 1, k) = offdiag(k);
  Eigen::SelfAdjointEigenSolver<Mat> es(T);
  QuadRule q{es.eigenvalues(), Vec(n)};
  for (int i = 0; i < n; ++i)
    q.w(i) = mu0 * es.eigenvectors()(0, i) * es.eigenvectors()(0, i);
  return q;
}

QuadRule legendre01(int n) {
  Vec d = Vec::Zero(n), e(n);
  e(0) = 0.0;
  for (int k = 1; k < n; ++k) e(k) = k / std::sqrt(4.0 * k * k - 1.0);
  QuadRule q = golub_welsch(d, e, 2.0);  // on [-1, 1]
  q.x = 0.5 * (q.x.array() + 1.0);
  q.w *= 0.5;
  return q;
}

QuadRule laguerre(int n, double alpha) {
  Vec d(n), e(n);
  e(0) = 0.0;
  for (int k = 0; k < n; ++k) d(k) = 2.0 * k + alpha + 1.0;
  for (int k = 1; k < n; ++k) e(k) = std::sqrt(k * (k + alpha));
  return golub_welsch(d, e, std::tgamma(alpha + 1.0));
}

}  // namespace

const QuadRule& gauss_legendre01(int n) {
  if (n < 1) throw std::invalid_argument("gauss_legendre01: n >= 1");
  static std::mutex m;
  static std::map<int, QuadRule> cache;
  std::lock_guard lock(m);
  auto it = cache.find(n);
  if (it == cache.end()) it = cache.emplace(n, legendre01(n)).first;
  return it->second;
}

const QuadRule& gauss_laguerre(int n, double alpha) {
  if (n < 1 || !(alpha > -1.0))
    throw std::invalid_argument("gauss_laguerre: n >= 1, alpha > -1");
  static std::mutex m;
  static std::map<std::pair<int, double>, QuadRule> cache;
  std::lock_guard lock(m);
  auto key = std::make_pair(n, alpha);
  auto it = cache.find(key);
  if (it == cache.end()) it = cache.emplace(key, laguerre(n, alpha)).first;
  return it->second;
}

}  // namespace qgem
