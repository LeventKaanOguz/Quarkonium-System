#pragma once
// Joint fit of two systems with shared parameters: the flavour independence of
// the static potential (one string tension for c c-bar and b b-bar). One
// parameter vector, the whitened residuals of both systems stacked.
#include <string>
#include <vector>

#include "fit/chi2_problem.hpp"

namespace qgem {

class JointProblem final : public LeastSquaresProblem {
 public:
  /// `shared`: names free in both problems, with the same transform. The joint
  /// vector is [u_a, the free parameters of b that are not shared].
  JointProblem(Chi2Problem& a, Chi2Problem& b,
               const std::vector<std::string>& shared);
  int n_data() const override { return a_.n_data() + b_.n_data(); }
  int n_free() const override { return a_.n_free() + b_.n_free() - k_; }
  Result<Vec> residuals(const Vec& u) override;
  Result<Mat> jacobian(const Vec& u) override;

  Vec u_a(const Vec& u) const { return u.head(a_.n_free()); }
  Vec u_b(const Vec& u) const;
  /// Joint vector from separate ones (shared values from a, or from b).
  Vec join(const Vec& ua, const Vec& ub, bool shared_from_b = false) const;
  /// Joint parameter names: a's free names, then b's unshared ones as
  /// "name(b)".
  std::vector<std::string> names() const;
  /// d theta / d u of each joint parameter (its transform).
  Vec dtheta_du(const Vec& u) const;
  int k() const { return k_; }  ///< number of shared parameters
  Chi2Problem& a() { return a_; }
  Chi2Problem& b() { return b_; }

 private:
  Chi2Problem& a_;
  Chi2Problem& b_;
  int k_ = 0;
  std::vector<int> b_src_;  ///< joint index of each free parameter of b
};

struct JointFit {
  FitResult fit;
  double sigma_a = 0, sigma_b = 0;  ///< floors per system, GeV
  double nu_a = 0, nu_b = 0;        ///< N_s - sum_{i in s} h_ii
  int iterations = 0;
};

/// Joint fit with a floor per system, each solved from its own share of the
/// degrees of freedom: chi2_s = nu_s = N_s - sum_{i in s} h_ii, with h the
/// leverages of the joint fit (nu_a + nu_b = N - p). Iterated to convergence.
JointFit fit_joint(JointProblem& jp, Vec u, double sigma_a, double sigma_b,
                   const FitOptions& opt = {});

/// Leverages h = diag(J (J^T J)^+ J^T) of a whitened Jacobian.
Vec leverages(const Mat& J);

}  // namespace qgem
