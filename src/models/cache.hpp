#pragma once
// Matrix cache: parameter-independent matrices are built once per basis; the
// few that depend on a nonlinear parameter (sigma, mu_s) are keyed by its value
// and evicted when the cache grows. One cache per thread (notes 11.12).
#include <map>
#include <tuple>

#include "basis/basis.hpp"
#include "core/types.hpp"
#include "me/radial.hpp"

namespace qgem {

class MatrixCache {
 public:
  explicit MatrixCache(BasisSpec spec = {}) : default_spec_(spec) {}

  /// The basis of partial wave L (one spec for all L).
  const BasisSpec& spec(int /*L*/) const { return default_spec_; }
  const Basis& basis(int L);

  const Mat& overlap(int L);
  const Mat& p2(int L);
  const Mat& r_power(int Lb, int Lk, int k);
  const Mat& gaussian(int Lb, int Lk, double beta, int k);
  const Mat& exp_power(int Lb, int Lk, double mu, int k);
  /// <V> for a named local function; `name` + `param` identify it (the caller
  /// owns the meaning).
  const Mat& local(int Lb, int Lk, int tag, double param,
                   const std::function<double(double)>& V);

 private:
  using Key =
      std::tuple<int, int, int, int, double>;  // kind, Lb, Lk, k/tag, param
  const Mat& get(const Key& key, bool param_dependent,
                 const std::function<Mat()>& make);

  BasisSpec default_spec_;
  std::map<int, Basis> bases_;
  std::map<Key, Mat> fixed_, dependent_;
};

}  // namespace qgem
