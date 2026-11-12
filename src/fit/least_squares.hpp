#pragma once
// Bounded-by-transform Levenberg-Marquardt on whitened residuals (notes
// Alg. 11.9), the theory floor (Alg. 11.10) and the uncertainty analysis at the
// minimum (Alg. 11.11). Nothing here knows about physics: the problems are
// abstract, so the same code is tested on a straight line (notes Toy 2).
#include <functional>

#include "core/error.hpp"
#include "core/types.hpp"
#include "models/params.hpp"

namespace qgem {

class LeastSquaresProblem {
 public:
  virtual ~LeastSquaresProblem() = default;
  virtual int n_data() const = 0;
  virtual int n_free() const = 0;
  /// Whitened residuals z_i = (y - mu)_i / s_i (or L^-1 (y - mu) for a
  /// covariance V = L L^T).
  virtual Result<Vec> residuals(const Vec& u) = 0;
  /// J = dz/du (N x p): minus the whitened model Jacobian.
  virtual Result<Mat> jacobian(const Vec& u) = 0;
};

struct FitResult {
  Vec u;
  Vec z;  ///< whitened residuals at u
  Mat J;  ///< dz/du at u
  double chi2 = 0;
  int iterations = 0;
  bool converged = false;
  std::string message;
};

struct FitOptions {
  double lambda0 = 1e-3;
  int max_iter = 200;
  double ftol = 1e-12;  ///< relative decrease of chi^2
  double xtol = 1e-10;  ///< relative step
  /// Alternative minimiser of sum z^2 (e.g. Minuit2 from the optional ROOT
  /// bridge). Empty: Levenberg-Marquardt.
  std::function<FitResult(LeastSquaresProblem&, const Vec&, const FitOptions&)>
      engine;
};

FitResult levenberg_marquardt(LeastSquaresProblem& problem, Vec u0,
                              const FitOptions& opt = {});
/// `opt.engine` if set, else Levenberg-Marquardt.
FitResult minimize(LeastSquaresProblem& problem, const Vec& u0,
                   const FitOptions& opt = {});

/// The theory floor at which the chi^2 of the residuals takes the value
/// `target`: chi2(x) = Delta^T (V_exp + x diag(w^2))^{-1} Delta with x =
/// sigma_th^2 (w empty: all 1), which decreases in x; solved by bisection.
/// Returns sigma_th (not squared); 0 if chi2(0) <= target. With target = nu
/// this is the fit's floor (chi2/nu = 1); with chi^2_nu quantiles it gives the
/// ends of the floor's confidence interval.
double solve_theory_floor(const Vec& Delta, const Mat& V_exp, double target,
                          const Vec& w = Vec());

struct Uncertainty {
  Mat C;       ///< covariance of the free parameters in the u variables
  Mat Ctheta;  ///< covariance of the free physical parameters (chain rule
               ///< through the transform)
  Vec sigma;   ///< sqrt(diag(Ctheta))
  Mat rho;     ///< correlation matrix of the physical parameters
  Vec sigma_conditional;   ///< 1/sqrt(F_aa): diagnostic only, never quote
  Vec fisher_eigenvalues;  ///< of the relative (dimensionless) Fisher matrix,
                           ///< ascending
  Mat fisher_eigenvectors;
  double condition = 0;
};

/// F = J^T J in u; C = F^-1 (pseudo-inverse with a warning flag if singular);
/// transformed to theta by C_theta = K C K^T with K = dtheta_free/du.
Uncertainty analyse(const FitResult& fit, const ParamSet& params);

}  // namespace qgem
