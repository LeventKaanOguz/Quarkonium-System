#include "stats/stats.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <numbers>
#include <thread>

#include "stats/rng.hpp"

namespace qgem {

double aic(double chi2, int p) { return chi2 + 2.0 * p; }
double bic(double chi2, int p, int N) {
  return chi2 + p * std::log(static_cast<double>(N));
}

double chi2_sf(double x, int nu) {
  if (x <= 0) return 1.0;
  const double a = 0.5 * nu, z = 0.5 * x;
  const double gln = std::lgamma(a);
  if (z < a + 1.0) {  // series for the lower function P
    double ap = a, sum = 1.0 / a, del = sum;
    for (int n = 0; n < 1000; ++n) {
      ap += 1.0;
      del *= z / ap;
      sum += del;
      if (std::fabs(del) < std::fabs(sum) * 1e-16) break;
    }
    return 1.0 - sum * std::exp(-z + a * std::log(z) - gln);
  }
  // continued fraction for Q (modified Lentz)
  double b = z + 1.0 - a, c = 1.0 / 1e-300, d = 1.0 / b, h = d;
  for (int i = 1; i < 1000; ++i) {
    const double an = -i * (i - a);
    b += 2.0;
    d = an * d + b;
    if (std::fabs(d) < 1e-300) d = 1e-300;
    c = b + an / c;
    if (std::fabs(c) < 1e-300) c = 1e-300;
    d = 1.0 / d;
    const double del = d * c;
    h *= del;
    if (std::fabs(del - 1.0) < 1e-16) break;
  }
  return std::exp(-z + a * std::log(z) - gln) * h;
}

double chi2_quantile(double q, int nu) {
  double lo = 0.0, hi = nu + 10.0 * std::sqrt(2.0 * nu) + 10.0;
  while (1.0 - chi2_sf(hi, nu) < q) hi *= 2;
  for (int it = 0; it < 200 && hi - lo > 1e-14 * hi; ++it) {
    const double mid = 0.5 * (lo + hi);
    (1.0 - chi2_sf(mid, nu) < q ? lo : hi) = mid;
  }
  return 0.5 * (lo + hi);
}

double t_coverage(double k, int nu) {
  const double v = nu;
  const double c = std::exp(std::lgamma(0.5 * (v + 1)) - std::lgamma(0.5 * v)) /
                   std::sqrt(v * std::numbers::pi);
  auto pdf = [&](double x) {
    return c * std::pow(1 + x * x / v, -0.5 * (v + 1));
  };
  const int n = 2000;  // even
  const double h = k / n;
  double s = pdf(0) + pdf(k);
  for (int i = 1; i < n; ++i) s += (i % 2 ? 4 : 2) * pdf(i * h);
  return 2 * s * h / 3;  // symmetric: twice the integral over [0, k]
}

double nested_test_pvalue(double chi2_null, double chi2_alt, int dp,
                          bool boundary) {
  const double p = chi2_sf(std::max(0.0, chi2_null - chi2_alt), dp);
  return boundary ? 0.5 * p : p;
}

LooResult leave_one_out(const FitResult& fit) {
  const Mat F = fit.J.transpose() * fit.J;
  const Mat C = F.inverse();
  const int N = static_cast<int>(fit.J.rows());
  LooResult r;
  r.leverage = Vec(N);
  r.loo_residual = Vec(N);
  r.pull = Vec(N);
  for (int i = 0; i < N; ++i) {
    const Vec Ji = fit.J.row(i);
    r.leverage(i) = Ji.dot(C * Ji);
    r.loo_residual(i) = fit.z(i) / (1.0 - r.leverage(i));
    r.pull(i) = fit.z(i) / std::sqrt(1.0 - r.leverage(i));
  }
  return r;
}

double propagate_variance(const Vec& g, const Mat& C) { return g.dot(C * g); }
double propagate_covariance(const Vec& g1, const Vec& g2, const Mat& C) {
  return g1.dot(C * g2);
}

double aitken_limit(double Pk, double Pk1, double Pk2) {
  const double den = Pk2 - 2 * Pk1 + Pk;
  if (std::fabs(den) < 1e-300) return Pk2;
  return Pk - (Pk1 - Pk) * (Pk1 - Pk) / den;
}

double Budget::total() const {
  return std::sqrt(sigma_par * sigma_par + sigma_num * sigma_num +
                   sigma_th * sigma_th);
}

Budget error_budget(double value, const Vec& g, const Mat& C, double sigma_num,
                    double sigma_th, const std::vector<double>& variants,
                    const std::vector<double>& models) {
  Budget b;
  b.value = value;
  b.sigma_par = std::sqrt(std::max(0.0, propagate_variance(g, C)));
  b.sigma_num = sigma_num;
  b.sigma_th = sigma_th;
  auto half_envelope = [](const std::vector<double>& v) {
    if (v.size() < 2) return 0.0;
    return 0.5 * (*std::max_element(v.begin(), v.end()) -
                  *std::min_element(v.begin(), v.end()));
  };
  b.sigma_amb = half_envelope(variants);
  b.model_spread = half_envelope(models);
  return b;
}

namespace {
// Fixes one free parameter of a problem at a value.
class FixedParameter final : public LeastSquaresProblem {
 public:
  FixedParameter(LeastSquaresProblem& base, int a, double value)
      : base_(base), a_(a), value_(value) {}
  int n_data() const override { return base_.n_data(); }
  int n_free() const override { return base_.n_free() - 1; }
  Result<Vec> residuals(const Vec& v) override {
    return base_.residuals(expand(v));
  }
  Result<Mat> jacobian(const Vec& v) override {
    auto J = base_.jacobian(expand(v));
    if (!J) return J;
    Mat out(J->rows(), J->cols() - 1);
    for (int c = 0, k = 0; c < J->cols(); ++c)
      if (c != a_) out.col(k++) = J->col(c);
    return out;
  }
  Vec expand(const Vec& v) const {
    Vec u(v.size() + 1);
    for (int c = 0, k = 0; c < u.size(); ++c)
      u(c) = (c == a_) ? value_ : v(k++);
    return u;
  }
  Vec reduce(const Vec& u) const {
    Vec v(u.size() - 1);
    for (int c = 0, k = 0; c < u.size(); ++c)
      if (c != a_) v(k++) = u(c);
    return v;
  }

 private:
  LeastSquaresProblem& base_;
  int a_;
  double value_;
};
}  // namespace

std::vector<std::pair<double, double>> profile(LeastSquaresProblem& problem,
                                               const FitResult& best, int a,
                                               std::span<const double> grid,
                                               const FitOptions& opt) {
  std::vector<std::pair<double, double>> out(grid.size());
  // walk outward from the best fit in both directions, warm-starting from the
  // neighbouring point
  std::vector<std::size_t> order(grid.size());
  for (std::size_t i = 0; i < grid.size(); ++i) order[i] = i;
  std::vector<std::size_t> up, down;
  for (std::size_t i : order) (grid[i] >= best.u(a) ? up : down).push_back(i);
  std::sort(up.begin(), up.end(),
            [&](auto x, auto y) { return grid[x] < grid[y]; });
  std::sort(down.begin(), down.end(),
            [&](auto x, auto y) { return grid[x] > grid[y]; });
  for (const auto* dir : {&up, &down}) {
    Vec start = best.u;
    for (std::size_t i : *dir) {
      FixedParameter fp(problem, a, grid[i]);
      const FitResult f = levenberg_marquardt(fp, fp.reduce(start), opt);
      out[i] = {grid[i], f.chi2 - best.chi2};
      start = fp.expand(f.u);
    }
  }
  return out;
}

BootstrapResult parametric_bootstrap(
    const Vec& mu_hat, const Mat& V, int B, std::uint64_t seed,
    const std::function<Result<std::pair<Vec, double>>(const Vec&)>& refit,
    int threads) {
  Eigen::LLT<Mat> llt(V);
  const Mat L = llt.matrixL();
  const int N = static_cast<int>(mu_hat.size());
  std::vector<Vec> us(B);
  std::vector<double> chi2(B, 0.0);
  std::vector<char> ok(B, 0);
  auto work = [&](int b0, int b1) {
    for (int b = b0; b < b1; ++b) {
      Rng rng(seed, static_cast<std::uint64_t>(b));
      Vec xi(N);
      for (int i = 0; i < N; ++i) xi(i) = rng.normal();
      auto r = refit(mu_hat + L * xi);
      if (r) {
        us[b] = r->first;
        chi2[b] = r->second;
        ok[b] = 1;
      }
    }
  };
  threads = std::max(1, threads);
  if (threads == 1) {
    work(0, B);
  } else {
    std::vector<std::jthread> pool;
    for (int t = 0; t < threads; ++t)
      pool.emplace_back(work, B * t / threads, B * (t + 1) / threads);
  }
  BootstrapResult out;
  int good = 0;
  for (int b = 0; b < B; ++b) good += ok[b];
  out.failed = B - good;
  if (good == 0) return out;
  const int p = static_cast<int>(
      us[std::find(ok.begin(), ok.end(), 1) - ok.begin()].size());
  out.u = Mat(good, p);
  out.chi2 = Vec(good);
  for (int b = 0, k = 0; b < B; ++b)
    if (ok[b]) {
      out.u.row(k) = us[b].transpose();
      out.chi2(k++) = chi2[b];
    }
  return out;
}

double gelman_rubin(std::span<const Vec> chains) {
  const int m = static_cast<int>(chains.size());
  const int n = static_cast<int>(chains[0].size());
  double W = 0, grand = 0;
  Vec means(m);
  for (int j = 0; j < m; ++j) {
    means(j) = chains[j].mean();
    grand += means(j) / m;
    W += (chains[j].array() - means(j)).square().sum() / (n - 1) / m;
  }
  const double Bv = n * (means.array() - grand).square().sum() / (m - 1);
  const double Vhat = (n - 1.0) / n * W + Bv / n;
  return std::sqrt(Vhat / W);
}

double integrated_autocorrelation_time(const Vec& x) {
  const int n = static_cast<int>(x.size());
  const Vec d = x.array() - x.mean();
  const double var = d.squaredNorm() / n;
  if (var <= 0) return 1.0;
  // Sokal's automatic window: smallest M with M >= 5 tau(M)
  double tau = 1.0;
  for (int M = 1; M < n / 2; ++M) {
    double rho = 0;
    for (int t = 0; t + M < n; ++t) rho += d(t) * d(t + M);
    rho /= (n - M) * var;
    tau += 2 * rho;
    if (M >= 5 * tau) break;
  }
  return std::max(tau, 1.0);
}

McmcResult metropolis(const std::function<double(const Vec&)>& log_post,
                      const Vec& u_hat, const Mat& C, int chains, int steps,
                      int burn, std::uint64_t seed, double overdispersion) {
  const int p = static_cast<int>(u_hat.size());
  Eigen::LLT<Mat> llt(C);
  const Mat LC = llt.matrixL();
  const Mat LS = std::sqrt(2.38 * 2.38 / p) * LC;
  McmcResult R;
  R.chains.resize(chains);
  auto run = [&](int j) {
    Rng rng(seed, static_cast<std::uint64_t>(j));
    Vec u = u_hat;
    double lp = -std::numeric_limits<double>::infinity();
    for (int tries = 0; tries < 100 && !std::isfinite(lp);
         ++tries) {  // over-dispersed, finite start
      Vec xi(p);
      for (int a = 0; a < p; ++a) xi(a) = rng.normal();
      u = u_hat + overdispersion * (LC * xi) / (1.0 + tries);
      lp = log_post(u);
    }
    Chain& ch = R.chains[j];
    ch.samples = Mat(steps - burn, p);
    ch.logp = Vec(steps - burn);
    int acc = 0;
    for (int t = 0; t < steps; ++t) {
      Vec xi(p);
      for (int a = 0; a < p; ++a) xi(a) = rng.normal();
      const Vec up = u + LS * xi;
      const double lpp = log_post(up);
      if (std::log(rng.uniform()) < lpp - lp) {
        u = up;
        lp = lpp;
        ++acc;
      }
      if (t >= burn) {
        ch.samples.row(t - burn) = u.transpose();
        ch.logp(t - burn) = lp;
      }
    }
    ch.acceptance = static_cast<double>(acc) / steps;
  };
  for (int j = 0; j < chains; ++j)
    run(j);  // serial: log_post may own mutable state
  R.rhat = Vec(p);
  R.tau_int = Vec(p);
  R.n_eff = Vec(p);
  R.mean = Vec(p);
  R.std = Vec(p);
  const int n = steps - burn;
  for (int a = 0; a < p; ++a) {
    std::vector<Vec> cs;
    double tau = 0, sum = 0, sum2 = 0;
    for (int j = 0; j < chains; ++j) {
      cs.push_back(R.chains[j].samples.col(a));
      tau += integrated_autocorrelation_time(cs.back()) / chains;
      sum += cs.back().sum();
      sum2 += cs.back().squaredNorm();
    }
    R.rhat(a) = chains > 1 ? gelman_rubin(cs) : 1.0;
    R.tau_int(a) = tau;
    R.n_eff(a) = static_cast<double>(chains) * n / tau;
    R.mean(a) = sum / (static_cast<double>(chains) * n);
    R.std(a) = std::sqrt(std::max(
        0.0, sum2 / (static_cast<double>(chains) * n) - R.mean(a) * R.mean(a)));
  }
  return R;
}

}  // namespace qgem
