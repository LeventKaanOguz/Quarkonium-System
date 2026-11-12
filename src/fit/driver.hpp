#pragma once
// High-level fitting workflow (notes 10.11): fit with the theory floor iterated
// to chi^2/nu = 1, exact leave-one-out refits, parametric bootstrap and
// posterior sampling of a Chi2Problem.
#include "fit/chi2_problem.hpp"
#include "stats/stats.hpp"

namespace qgem {

struct FloorFit {
  FitResult fit;
  double sigma_th = 0;  ///< GeV
  int iterations = 0;
};

/// Alg. 11.10: iterate [fit at sigma_th -> residuals -> new sigma_th from
/// chi^2/nu = 1, full PDG covariance] to convergence. With `fixed_floor >= 0`
/// the floor is not iterated (used for model comparison with a common floor).
FloorFit fit_with_floor(Chi2Problem& prob, Vec u0, double sigma_th0 = 0.010,
                        const FitOptions& opt = {}, double fixed_floor = -1.0);

/// One held-out refit: refit without the held-out data (indices into
/// prob.data().d, all in the fit; basis fixed) and predict every datum. With
/// `refloor` the floor is re-solved on the kept data (chi^2/nu = 1), so that
/// the held-out data do not set their own predictive scale; with nu = 0 left
/// it cannot be, and sigma_th is NaN (the refit keeps the full-data weights).
/// Without `refloor` the floor is held at the full fit's value.
struct HeldOutFit {
  Vec mu;           ///< model value of every datum, GeV (empty: refit failed)
  Vec sigma_par;    ///< its parametric error from the refit's covariance, GeV
  Vec theta;        ///< refitted physical parameters
  double sigma_th;  ///< floor of the refit, GeV (see above)
};
HeldOutFit refit_held_out(Chi2Problem& prob, const Vec& u_hat,
                          const std::vector<int>& heldout,
                          const FitOptions& opt = {}, bool refloor = false);

/// Exact leave-one-out: one held-out refit per fitted datum, in fit order.
/// The predictive variance of miss_i is sigma_th_i^2 w_i^2 + s_exp,i^2 +
/// sigma_par_i^2, all from data without datum i.
struct LooRefits {
  Vec miss;       ///< y_i - mu_i^{(-i)}, GeV (NaN if the refit failed)
  Mat theta;      ///< n_fit x n_total refitted parameters (NaN rows: failed)
  Vec sigma_th;   ///< floor of each refit, GeV
  Vec sigma_par;  ///< parametric error of mu_i^{(-i)}, GeV
};
LooRefits exact_leave_one_out_refits(Chi2Problem& prob, const Vec& u_hat,
                                     const FitOptions& opt = {},
                                     bool refloor = false);
/// The misses alone, floor held fixed.
Vec exact_leave_one_out(Chi2Problem& prob, const Vec& u_hat,
                        const FitOptions& opt = {});

/// Refit without the held-out data (floor held at the full fit's value: the
/// random-split RMS scores depend on it only through the relative weights)
/// and return the model value of every datum. Error if the refit fails.
Result<Vec> refit_without(Chi2Problem& prob, const Vec& u_hat,
                          const std::vector<int>& heldout,
                          const FitOptions& opt = {});

/// Leave-group-out: one held-out refit without a whole group; the misses
/// y - mu^{(-group)} of the group's data in the given order, their parametric
/// errors and the refit's floor (NaN where unavailable).
struct GroupRefit {
  Vec miss, sigma_par;  ///< GeV
  double sigma_th;      ///< GeV
};
GroupRefit leave_group_out(Chi2Problem& prob, const Vec& u_hat,
                           const std::vector<int>& group,
                           const FitOptions& opt = {}, bool refloor = false);

/// Closure test of an omitted term. The model's own values at u_hat, shifted
/// by the omitted contribution (`shift`, every datum, GeV: truth = mu -
/// shift), are refitted. R = (mu - shift) - mu' is the part of the term the
/// parameters cannot absorb, as data - model residuals for every datum.
/// sigma_a solves sum R^2/(s^2 + w^2 sigma_a^2) = nu (the fit's floor
/// equation), with the weights made self-consistent by a second refit at
/// sigma_a. Empty R if a refit fails.
struct Truncation {
  Vec R;               ///< GeV, every datum
  double sigma_a = 0;  ///< GeV
  Vec theta;           ///< parameters that absorb the shifted values
};
Truncation absorb_shift(const Chi2Problem& prob, const Vec& u_hat,
                        const Vec& shift, const FitOptions& opt = {});

/// Delete-one jackknife standard error of each column of `samples` (rows with
/// NaN are skipped): sqrt((n-1)/n * sum_i (x_i - mean)^2).
Vec jackknife_sigma(const Mat& samples);

/// Parametric bootstrap of the whole fit. Each worker thread uses its own copy
/// of the problem.
BootstrapResult bootstrap_fit(Chi2Problem& prob, const Vec& u_hat, int B,
                              std::uint64_t seed, int threads = 1,
                              const FitOptions& opt = {});

/// Posterior sampling, flat priors in the u variables (flat in ln theta for
/// log-transformed parameters).
McmcResult sample_posterior(Chi2Problem& prob, const FitResult& fit,
                            const ParamSet& params, int chains, int steps,
                            int burn, std::uint64_t seed);

/// Numerical error of every datum's predicted mass: |M(2n-1) - M(n)| from
/// nested basis doubling (all partial waves), at the physical parameters theta.
/// The criterion of the notes is sigma_num <= 1e-2 min_i s_i.
Vec numerical_errors(Chi2Problem& prob, const Vec& theta,
                     const BasisSpec& spec);

}  // namespace qgem
