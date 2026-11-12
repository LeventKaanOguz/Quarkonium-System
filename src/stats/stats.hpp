#pragma once
// Statistical tools of notes 9-10: leave-one-out, information criteria, profile
// likelihood, parametric bootstrap, Metropolis sampling with convergence
// diagnostics, error propagation, Aitken extrapolation.
#include <functional>
#include <span>
#include <utility>
#include <vector>

#include "fit/least_squares.hpp"

namespace qgem {

// ---- model comparison
// ---------------------------------------------------------------------------------

double aic(double chi2, int p);
double bic(double chi2, int p, int N);
/// P(chi^2_nu >= x): regularized upper incomplete gamma Q(nu/2, x/2).
double chi2_sf(double x, int nu);
/// The q-quantile of chi^2_nu (0 < q < 1), by bisection on chi2_sf.
double chi2_quantile(double q, int nu);
/// P(|T| < k) for T ~ t_nu (Simpson's rule on the density): the coverage of a
/// +-k interval when the scale is estimated with nu degrees of freedom.
double t_coverage(double k, int nu);
/// Likelihood-ratio test p-value for nested models (Wilks): P(chi^2_dp >=
/// chi2_null - chi2_alt). `boundary`: the null value sits on the boundary of
/// the parameter space (Chernoff: half the naive p).
double nested_test_pvalue(double chi2_null, double chi2_alt, int dp,
                          bool boundary = false);

// ---- leave-one-out (notes Result 10.4)
// ------------------------------------------------------------------

struct LooResult {
  Vec leverage;  ///< h_i, diag of the hat matrix J C J^T (whitened); sum = p
  Vec loo_residual;  ///< (y_i - mu_i^{(-i)}) = Delta_i / (1 - h_i), in whitened
                     ///< units
  Vec pull;  ///< z_i / sqrt(1 - h_i), whitened residual divided by sqrt(1 -
             ///< h_i)
};
/// From one fit, linear order. `fit.J` and `fit.z` are the whitened Jacobian
/// and residuals at the minimum.
LooResult leave_one_out(const FitResult& fit);

// ---- error propagation
// ----------------------------------------------------------------------------------

/// sigma^2(P) = g^T C g  (eq. 10.14), and the covariance of two predictions.
double propagate_variance(const Vec& g, const Mat& C);
double propagate_covariance(const Vec& g1, const Vec& g2, const Mat& C);
/// Aitken's delta-squared limit of a geometric sequence P_k, P_{k+1}, P_{k+2}
/// (eq. 10.20).
double aitken_limit(double Pk, double Pk1, double Pk2);

// ---- error budget (notes 10.8)
// ---------------------------------------------------------------------------

/// Every quoted number carries separate contributions. Only independent,
/// Gaussian-like ones add in quadrature
/// (`total()`); the convention ambiguity and the model dependence are reported
/// on their own, never added.
struct Budget {
  double value = 0;
  double sigma_par =
      0;  ///< sqrt(g^T C g): how well the data constrain the parameters
  double sigma_num = 0;  ///< basis/quadrature: |P(2n-1) - P(n)|
  double sigma_th =
      0;  ///< theory floor (masses) -- 0 for quantities without one
  double sigma_amb = 0;     ///< half the envelope over reasonable convention
                            ///< choices (M in widths, alpha_s scale, ...)
  double model_spread = 0;  ///< half the envelope of the same prediction across
                            ///< fitted models: not a statistical error
  double total() const;     ///< sqrt(par^2 + num^2 + th^2)
};

/// g: gradient of the prediction with respect to the fitted parameters, C their
/// covariance (same variables).
Budget error_budget(double value, const Vec& g, const Mat& C, double sigma_num,
                    double sigma_th,
                    const std::vector<double>& convention_variants = {},
                    const std::vector<double>& other_models = {});

// ---- profile likelihood
// ---------------------------------------------------------------------------------

/// chi^2_prof at the given values of free parameter `a` (in u variables),
/// refitting all others from the neighbouring grid point. Returns (value,
/// chi^2_prof - chi^2_min) pairs.
std::vector<std::pair<double, double>> profile(LeastSquaresProblem& problem,
                                               const FitResult& best, int a,
                                               std::span<const double> grid,
                                               const FitOptions& opt = {});

// ---- parametric bootstrap (Alg. 11.12)
// ------------------------------------------------------------------

struct BootstrapResult {
  Mat u;     ///< B x p refitted parameters (u variables)
  Vec chi2;  ///< chi^2 of each replica
  int failed = 0;
};
/// y*(b) = mu_hat + L xi with V = L L^T; `refit(y_star)` returns the fitted u
/// (or an error). Replica b uses its own RNG stream (seed, b), so results do
/// not depend on the number of threads.
BootstrapResult parametric_bootstrap(
    const Vec& mu_hat, const Mat& V, int B, std::uint64_t seed,
    const std::function<Result<std::pair<Vec, double>>(const Vec&)>& refit,
    int threads = 1);

// ---- Metropolis-Hastings (Alg. 11.13)
// --------------------------------------------------------------------

struct Chain {
  Mat samples;  ///< n x p (after burn-in)
  Vec logp;
  double acceptance = 0;
};
struct McmcResult {
  std::vector<Chain> chains;
  Vec rhat;       ///< Gelman-Rubin per parameter (10.17)
  Vec tau_int;    ///< integrated autocorrelation time per parameter (mean over
                  ///< chains)
  Vec n_eff;      ///< m n / tau_int
  Vec mean, std;  ///< pooled
};
/// `log_post(u)` = -chi^2/2 (+ log prior); return -inf outside the support.
/// Proposal covariance 2.38^2/p C.
McmcResult metropolis(const std::function<double(const Vec&)>& log_post,
                      const Vec& u_hat, const Mat& C, int chains, int steps,
                      int burn, std::uint64_t seed,
                      double overdispersion = 3.0);

double gelman_rubin(std::span<const Vec> chains);  ///< one parameter: chains as
                                                   ///< vectors of samples
double integrated_autocorrelation_time(const Vec& x);

}  // namespace qgem
