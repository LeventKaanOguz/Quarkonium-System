#include "models/cache.hpp"

namespace qgem {

const Basis& MatrixCache::basis(int L) {
  auto it = bases_.find(L);
  if (it == bases_.end()) it = bases_.emplace(L, make_basis(L, spec(L))).first;
  return it->second;
}

const Mat& MatrixCache::get(const Key& key, bool param_dependent,
                            const std::function<Mat()>& make) {
  auto& m = param_dependent ? dependent_ : fixed_;
  auto it = m.find(key);
  if (it != m.end()) return it->second;
  if (param_dependent && m.size() > 256) m.clear();
  return m.emplace(key, make()).first->second;
}

const Mat& MatrixCache::overlap(int L) {
  return get({0, L, L, 0, 0.0}, false, [&] { return qgem::overlap(basis(L)); });
}
const Mat& MatrixCache::p2(int L) {
  return get({1, L, L, 0, 0.0}, false, [&] { return qgem::p2(basis(L)); });
}
const Mat& MatrixCache::r_power(int Lb, int Lk, int k) {
  return get({2, Lb, Lk, k, 0.0}, false,
             [&] { return qgem::r_power(basis(Lb), basis(Lk), k); });
}
const Mat& MatrixCache::gaussian(int Lb, int Lk, double beta, int k) {
  return get({3, Lb, Lk, k, beta}, true,
             [&] { return qgem::gaussian(basis(Lb), basis(Lk), beta, k); });
}
const Mat& MatrixCache::exp_power(int Lb, int Lk, double mu, int k) {
  return get({4, Lb, Lk, k, mu}, true,
             [&] { return qgem::exp_power(basis(Lb), basis(Lk), mu, k); });
}
const Mat& MatrixCache::local(int Lb, int Lk, int tag, double param,
                              const std::function<double(double)>& V) {
  return get({5, Lb, Lk, tag, param}, true,
             [&] { return qgem::local(basis(Lb), basis(Lk), V); });
}

}  // namespace qgem
