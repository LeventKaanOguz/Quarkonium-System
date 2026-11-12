#include "fit/driver.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <vector>

namespace qgem {

namespace {

/// Alternate [minimise at sigma_th -> re-solve sigma_th from chi^2/nu = 1 with
/// the full PDG covariance] to convergence; `min(prob, u0)` minimises.
template <class Minimiser>
FloorFit iterate_floor(Chi2Problem& prob, Vec u, double sth, Minimiser&& min) {
  FloorFit out;
  const double nu = prob.n_data() - prob.n_free();
  for (int it = 0; it < 30; ++it) {
    out.iterations = it + 1;
    prob.set_sigma_th(sth);
    out.fit = min(prob, u);
    u = out.fit.u;
    auto raw = prob.raw_residuals(u);
    if (!raw) break;
    const double next = solve_theory_floor(*raw, prob.exp_covariance(), nu,
                                           prob.floor_weights());
    const bool done = std::fabs(next - sth) < 1e-9;
    sth = next;
    if (done) break;
  }
  prob.set_sigma_th(sth);
  out.fit = min(prob, u);
  out.sigma_th = sth;
  return out;
}

}  // namespace

FloorFit fit_with_floor(Chi2Problem& prob, Vec u, double sth,
                        const FitOptions& opt, double fixed) {
  auto min = [&](Chi2Problem& p, const Vec& u0) {
    return minimize(p, u0, opt);
  };
  if (fixed >= 0) {
    FloorFit out;
    prob.set_sigma_th(fixed);
    out.fit = min(prob, u);
    out.sigma_th = fixed;
    return out;
  }
  return iterate_floor(prob, std::move(u), sth, min);
}

HeldOutFit refit_held_out(Chi2Problem& prob, const Vec& u_hat,
                          const std::vector<int>& heldout,
                          const FitOptions& opt, bool refloor) {
  HeldOutFit out;
  out.sigma_th = std::numeric_limits<double>::quiet_NaN();
  Dataset reduced = prob.data();
  for (int k : heldout) reduced.d.at(k).in_fit = false;
  Chi2Problem sub(prob.model(), reduced, prob.params(), prob.cache(),
                  prob.sigma_th());
  auto lm = [&](Chi2Problem& p, const Vec& u0) {
    return levenberg_marquardt(p, u0, opt);
  };
  FitResult f;
  if (!refloor) {
    f = lm(sub, u_hat);
    out.sigma_th = prob.sigma_th();
  } else if (sub.n_data() - sub.n_free() >= 1) {
    FloorFit ff = iterate_floor(sub, u_hat, prob.sigma_th(), lm);
    f = std::move(ff.fit);
    out.sigma_th = ff.sigma_th;
  } else {
    f = lm(sub, u_hat);  // nu = 0: the kept data cannot fix a floor
  }
  if (!f.J.size()) return out;
  out.theta = prob.params().theta(f.u);
  auto pg = sub.predictions_with_gradients(out.theta);
  if (!pg) return out;
  const Mat K = prob.params().dtheta_du(f.u);
  const Mat C = K * analyse(f, prob.params()).C * K.transpose();
  out.mu = pg->first;
  out.sigma_par.resize(out.mu.size());
  for (int i = 0; i < out.mu.size(); ++i)
    out.sigma_par(i) = std::sqrt(std::max(
        0.0, propagate_variance(Vec(pg->second.row(i).transpose()), C)));
  return out;
}

LooRefits exact_leave_one_out_refits(Chi2Problem& prob, const Vec& u_hat,
                                     const FitOptions& opt, bool refloor) {
  const int N = prob.n_measurements();
  const double nan = std::numeric_limits<double>::quiet_NaN();
  LooRefits out{Vec::Constant(N, nan),
                Mat::Constant(N, prob.params().n_total(), nan),
                Vec::Constant(N, nan), Vec::Constant(N, nan)};
  int i = 0;
  for (std::size_t k = 0; k < prob.data().d.size(); ++k) {
    if (!prob.data().d[k].in_fit) continue;
    const int kk = static_cast<int>(k);
    const HeldOutFit h = refit_held_out(prob, u_hat, {kk}, opt, refloor);
    out.sigma_th(i) = h.sigma_th;
    if (h.mu.size() && std::isfinite(h.mu(kk))) {
      out.miss(i) = prob.data().d[k].y - h.mu(kk);
      out.sigma_par(i) = h.sigma_par(kk);
      out.theta.row(i) = h.theta.transpose();
    }
    ++i;
  }
  return out;
}

Result<Vec> refit_without(Chi2Problem& prob, const Vec& u_hat,
                          const std::vector<int>& heldout,
                          const FitOptions& opt) {
  Dataset reduced = prob.data();
  for (int k : heldout) reduced.d.at(k).in_fit = false;
  Chi2Problem sub(prob.model(), reduced, prob.params(), prob.cache(),
                  prob.sigma_th());
  const FitResult f = levenberg_marquardt(sub, u_hat, opt);
  if (!f.J.size()) return fail("refit failed: " + f.message);
  return sub.predictions(prob.params().theta(f.u));
}

GroupRefit leave_group_out(Chi2Problem& prob, const Vec& u_hat,
                           const std::vector<int>& group, const FitOptions& opt,
                           bool refloor) {
  const int n = static_cast<int>(group.size());
  const double nan = std::numeric_limits<double>::quiet_NaN();
  GroupRefit out{Vec::Constant(n, nan), Vec::Constant(n, nan), nan};
  const HeldOutFit h = refit_held_out(prob, u_hat, group, opt, refloor);
  out.sigma_th = h.sigma_th;
  if (!h.mu.size()) return out;
  for (int g = 0; g < n; ++g) {
    const double m = h.mu(group[g]);
    if (!std::isfinite(m)) continue;
    out.miss(g) = prob.data().d[group[g]].y - m;
    out.sigma_par(g) = h.sigma_par(group[g]);
  }
  return out;
}

Truncation absorb_shift(const Chi2Problem& prob, const Vec& u_hat,
                        const Vec& shift, const FitOptions& opt) {
  Chi2Problem c(prob);
  const auto mu = c.predictions(c.params().theta(u_hat));
  if (!mu) return {};
  std::vector<int> fit;
  for (std::size_t k = 0; k < c.data().d.size(); ++k)
    if (c.data().d[k].in_fit) fit.push_back(static_cast<int>(k));
  const Vec target = *mu - shift;
  Vec y(static_cast<int>(fit.size()));
  for (std::size_t i = 0; i < fit.size(); ++i)
    y(static_cast<int>(i)) = target(fit[i]);
  c.set_measurements(y);
  const int nu = c.n_data() - c.n_free();
  Truncation out;
  Vec u = u_hat;
  for (int pass = 0; pass < 2; ++pass) {  // the fit's weights, then sigma_a
    const FitResult f = levenberg_marquardt(c, u, opt);
    if (!f.J.size()) return {};
    u = f.u;
    const auto mu2 = c.predictions(c.params().theta(u));
    if (!mu2) return {};
    out.R = target - *mu2;
    Vec Rf(y.size());
    for (std::size_t i = 0; i < fit.size(); ++i)
      Rf(static_cast<int>(i)) = out.R(fit[i]);
    out.sigma_a =
        solve_theory_floor(Rf, c.exp_covariance(), nu, c.floor_weights());
    c.set_sigma_th(out.sigma_a);
  }
  out.theta = c.params().theta(u);
  return out;
}

Vec exact_leave_one_out(Chi2Problem& prob, const Vec& u_hat,
                        const FitOptions& opt) {
  return exact_leave_one_out_refits(prob, u_hat, opt).miss;
}

Vec jackknife_sigma(const Mat& samples) {
  Vec out(samples.cols());
  for (int c = 0; c < samples.cols(); ++c) {
    std::vector<double> x;
    for (int r = 0; r < samples.rows(); ++r)
      if (samples.row(r).allFinite()) x.push_back(samples(r, c));
    const double n = static_cast<double>(x.size());
    if (n < 2) {
      out(c) = std::numeric_limits<double>::quiet_NaN();
      continue;
    }
    double mean = 0, ss = 0;
    for (double v : x) mean += v / n;
    for (double v : x) ss += (v - mean) * (v - mean);
    out(c) = std::sqrt((n - 1) / n * ss);
  }
  return out;
}

BootstrapResult bootstrap_fit(Chi2Problem& prob, const Vec& u_hat, int B,
                              std::uint64_t seed, int threads,
                              const FitOptions& opt) {
  const Vec theta_hat = prob.params().theta(u_hat);
  const Vec mu_all = *prob.predictions(theta_hat);
  Vec mu(prob.n_data());
  int i = 0;
  for (std::size_t k = 0; k < prob.data().d.size(); ++k)
    if (prob.data().d[k].in_fit) mu(i++) = mu_all(static_cast<int>(k));
  const Mat V = prob.covariance();
  auto refit = [&](const Vec& ystar) -> Result<std::pair<Vec, double>> {
    thread_local std::unique_ptr<Chi2Problem> local;
    thread_local const Chi2Problem* owner = nullptr;
    if (!local || owner != &prob) {
      local = std::make_unique<Chi2Problem>(prob);
      owner = &prob;
    }
    local->set_measurements(ystar);
    const FitResult f = levenberg_marquardt(*local, u_hat, opt);
    if (!f.converged) return fail("bootstrap refit did not converge");
    return std::make_pair(f.u, f.chi2);
  };
  return parametric_bootstrap(mu, V, B, seed, refit, threads);
}

McmcResult sample_posterior(Chi2Problem& prob, const FitResult& fit,
                            const ParamSet& params, int chains, int steps,
                            int burn, std::uint64_t seed) {
  const Uncertainty U = analyse(fit, params);
  auto log_post = [&](const Vec& u) {
    auto z = prob.residuals(u);
    return z ? -0.5 * z->squaredNorm()
             : -std::numeric_limits<double>::infinity();
  };
  return metropolis(log_post, fit.u, U.C, chains, steps, burn, seed);
}

McmcResult sample_posterior_floor(const Chi2Problem& prob, const FitResult& fit,
                                  double sigma_th, const ParamSet& params,
                                  int chains, int steps, int burn,
                                  std::uint64_t seed) {
  Chi2Problem c(prob);
  const int p = params.n_free();
  const Uncertainty U = analyse(fit, params);
  const double lo = std::log(1e-5), hi = std::log(0.2);
  auto log_post = [&](const Vec& v) {
    const double s = v(p);
    if (s < lo || s > hi) return -std::numeric_limits<double>::infinity();
    c.set_sigma_th(std::exp(s));
    const auto z = c.residuals(v.head(p));
    return z ? -0.5 * (z->squaredNorm() + c.log_det_covariance())
             : -std::numeric_limits<double>::infinity();
  };
  const int nu = std::max(1, prob.n_data() - p);
  Mat C = Mat::Zero(p + 1, p + 1);
  C.topLeftCorner(p, p) = U.C;
  C(p, p) = 1.0 / (2.0 * nu);  // sd of ln sigma ~ 1/sqrt(2 nu)
  Vec v0(p + 1);
  v0 << fit.u, std::log(std::max(sigma_th, 2e-5));
  return metropolis(log_post, v0, C, chains, steps, burn, seed);
}

PosteriorCheck posterior_predictive(const Chi2Problem& prob, const Mat& samples,
                                    int max_draws) {
  PosteriorCheck out;
  Chi2Problem c(prob);
  const int p = static_cast<int>(samples.cols()) - 1;
  const int N = c.n_measurements();
  const int step =
      std::max<int>(1, static_cast<int>(samples.rows()) / max_draws);
  for (int r = 0; r < samples.rows(); r += step) {
    c.set_sigma_th(std::exp(samples(r, p)));
    const auto z = c.residuals(Vec(samples.row(r).head(p).transpose()));
    if (!z) continue;
    const Vec zz = z->head(N);
    out.ppp_chi2 += chi2_sf(zz.squaredNorm(), N);
    const double m = zz.cwiseAbs().maxCoeff();
    out.ppp_max += 1 - std::pow(1 - std::erfc(m / std::sqrt(2.0)), N);
    ++out.draws;
  }
  if (out.draws) out.ppp_chi2 /= out.draws, out.ppp_max /= out.draws;
  return out;
}

Vec numerical_errors(Chi2Problem& prob, const Vec& theta,
                     const BasisSpec& spec) {
  BasisSpec fine = spec;
  fine.n = 2 * spec.n - 1;
  Chi2Problem a(prob.model(), prob.data(), prob.params(), MatrixCache(spec),
                prob.sigma_th());
  Chi2Problem b(prob.model(), prob.data(), prob.params(), MatrixCache(fine),
                prob.sigma_th());
  const auto ma = a.predictions(theta), mb = b.predictions(theta);
  if (!ma || !mb)
    throw std::runtime_error("numerical_errors: prediction failed");
  return (*mb - *ma).cwiseAbs();
}

}  // namespace qgem
