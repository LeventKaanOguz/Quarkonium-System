// M6: uncertainties and validation, verified exactly on the straight-line toy
// of the notes (9.6, 9.12) and against known distributions.
#include <gtest/gtest.h>

#include <cmath>

#include "fit/least_squares.hpp"
#include "stats/rng.hpp"
#include "stats/stats.hpp"

using namespace qgem;

namespace {
class Line final : public LeastSquaresProblem {
 public:
  explicit Line(Vec y = (Vec(3) << 1, 2, 4).finished()) : y_(std::move(y)) {}
  int n_data() const override { return 3; }
  int n_free() const override { return 2; }
  Result<Vec> residuals(const Vec& u) override {
    Vec z(3);
    for (int i = 0; i < 3; ++i) z(i) = y_(i) - (u(0) + u(1) * i);
    return z;
  }
  Result<Mat> jacobian(const Vec&) override {
    Mat J(3, 2);
    for (int i = 0; i < 3; ++i) J.row(i) << -1.0, -static_cast<double>(i);
    return J;
  }
  Vec y_;
};
}  // namespace

TEST(Rng, ReproducibleAndNormallyDistributed) {
  Rng a(7, 3), b(7, 3), c(7, 4);
  for (int i = 0; i < 10; ++i) EXPECT_EQ(a.normal(), b.normal());
  EXPECT_NE(Rng(7, 3).normal(), c.normal());
  Rng r(123);
  const int n = 200000;
  double s = 0, s2 = 0, s4 = 0;
  for (int i = 0; i < n; ++i) {
    const double x = r.normal();
    s += x;
    s2 += x * x;
    s4 += x * x * x * x;
  }
  EXPECT_NEAR(s / n, 0.0, 0.01);
  EXPECT_NEAR(s2 / n, 1.0, 0.01);
  EXPECT_NEAR(s4 / n, 3.0, 0.06);
}

TEST(Stats, Chi2SurvivalFunctionAndInformationCriteria) {
  // known values: P(chi2_1 >= 3.841459) = 0.05, P(chi2_4 >= 9.487729) = 0.05,
  // P(chi2_2 >= x) = exp(-x/2)
  EXPECT_NEAR(chi2_sf(3.841458820694124, 1), 0.05, 1e-12);
  EXPECT_NEAR(chi2_sf(9.487729036781154, 4), 0.05, 1e-12);
  EXPECT_NEAR(chi2_sf(5.0, 2), std::exp(-2.5), 1e-14);
  EXPECT_NEAR(chi2_sf(50.0, 10), 2.6694e-7, 1e-10);
  EXPECT_NEAR(nested_test_pvalue(10.0, 7.67, 1), chi2_sf(2.33, 1), 1e-15);
  EXPECT_NEAR(nested_test_pvalue(10.0, 7.67, 1, true), 0.5 * chi2_sf(2.33, 1),
              1e-15);
  // the notes' worked example: dchi2 = 2.33 for one extra parameter, Wilks p =
  // 0.13
  EXPECT_NEAR(chi2_sf(2.33, 1), 0.127, 0.001);
  EXPECT_NEAR(aic(4.0, 4), 12.0, 1e-15);
  EXPECT_NEAR(bic(4.0, 4, 8), 4.0 + 4 * std::log(8.0), 1e-15);  // 12.32
}

TEST(Stats, LeaveOneOutOnTheStraightLineIsExact) {
  Line prob;
  const FitResult f = levenberg_marquardt(prob, Vec::Zero(2));
  const LooResult l = leave_one_out(f);
  EXPECT_NEAR(l.leverage(0), 5.0 / 6.0, 1e-9);
  EXPECT_NEAR(l.leverage(1), 1.0 / 3.0, 1e-9);
  EXPECT_NEAR(l.leverage(2), 5.0 / 6.0, 1e-9);
  EXPECT_NEAR(l.leverage.sum(), 2.0, 1e-9);  // sum h_i = p
  // leave-one-out misses y_i - mu_i^(-i) = (+1, -0.5, +1) (notes 9.12); here z
  // = y - mu
  EXPECT_NEAR(l.loo_residual(0), 1.0, 1e-9);
  EXPECT_NEAR(l.loo_residual(1), -0.5, 1e-9);
  EXPECT_NEAR(l.loo_residual(2), 1.0, 1e-9);
  // exact refits: drop point i, fit the line through the other two, predict
  // point i
  const double y[] = {1, 2, 4};
  for (int i = 0; i < 3; ++i) {
    int j = (i + 1) % 3, k = (i + 2) % 3;
    const double slope = (y[k] - y[j]) / (k - j), icpt = y[j] - slope * j;
    const double miss = y[i] - (icpt + slope * i);
    EXPECT_NEAR(l.loo_residual(i), miss, 1e-9) << i;
  }
}

TEST(Stats, ProfileLikelihoodOfALinearModelIsExactlyQuadratic) {
  Line prob;
  const FitResult f = levenberg_marquardt(prob, Vec::Zero(2));
  const double sa = std::sqrt(5.0 / 6.0);
  std::vector<double> grid = {f.u(0) - 2 * sa, f.u(0) - sa, f.u(0), f.u(0) + sa,
                              f.u(0) + 2 * sa};
  const auto pr = profile(prob, f, 0, grid);
  for (std::size_t i = 0; i < grid.size(); ++i) {
    const double t = (grid[i] - f.u(0)) / sa;
    EXPECT_NEAR(pr[i].second, t * t, 1e-8)
        << i;  // Delta chi^2 = 1 at +-sigma_marginal (notes 10.4.3)
  }
}

TEST(Stats, BootstrapRecoversTheCovarianceOfALinearFit) {
  Line prob;
  const FitResult f = levenberg_marquardt(prob, Vec::Zero(2));
  const Vec mu =
      (Vec(3) << f.u(0), f.u(0) + f.u(1), f.u(0) + 2 * f.u(1)).finished();
  auto refit = [&](const Vec& ystar) -> Result<std::pair<Vec, double>> {
    Line p(ystar);
    FitResult r = levenberg_marquardt(p, f.u);
    return std::make_pair(r.u, r.chi2);
  };
  const Mat V = Mat::Identity(3, 3);
  const BootstrapResult b1 = parametric_bootstrap(mu, V, 20000, 99, refit, 1);
  const BootstrapResult b4 = parametric_bootstrap(mu, V, 20000, 99, refit, 4);
  EXPECT_EQ(b1.failed, 0);
  EXPECT_LT((b1.u - b4.u).cwiseAbs().maxCoeff(),
            1e-15);  // independent of the thread count
  const Vec mean = b1.u.colwise().mean();
  EXPECT_NEAR(mean(0), f.u(0), 0.02);
  EXPECT_NEAR(mean(1), f.u(1), 0.015);
  const Mat d = b1.u.rowwise() - mean.transpose();
  const Mat cov = d.transpose() * d / (b1.u.rows() - 1);
  EXPECT_NEAR(std::sqrt(cov(0, 0)), std::sqrt(5.0 / 6.0), 0.01);
  EXPECT_NEAR(std::sqrt(cov(1, 1)), std::sqrt(0.5), 0.008);
  EXPECT_NEAR(cov(0, 1) / std::sqrt(cov(0, 0) * cov(1, 1)),
              -3.0 / std::sqrt(15.0), 0.012);
  // chi^2 of the replicas follows chi^2_1 (nu = 1): mean 1
  EXPECT_NEAR(b1.chi2.mean(), 1.0, 0.04);  // ~3.5 standard errors at B = 20000
}

TEST(Stats, MetropolisSamplesTheGaussianPosteriorWithGoodDiagnostics) {
  Line prob;
  const FitResult f = levenberg_marquardt(prob, Vec::Zero(2));
  Mat C(2, 2);
  C << 5.0 / 6.0, -0.5, -0.5, 0.5;
  auto lp = [&](const Vec& u) {
    const Vec z = *prob.residuals(u);
    return -0.5 * z.squaredNorm();
  };
  const McmcResult m = metropolis(lp, f.u, C, 4, 6000, 1000, 2026);
  EXPECT_NEAR(m.chains[0].acceptance, 0.35, 0.15);  // 2.38^2/p scaling: ~0.3
  for (int a = 0; a < 2; ++a) {
    EXPECT_LT(m.rhat(a), 1.01) << a;
    EXPECT_GT(m.n_eff(a), 1000.0) << a;
  }
  EXPECT_NEAR(m.mean(0), 5.0 / 6.0, 0.06);
  EXPECT_NEAR(m.mean(1), 1.5, 0.05);
  EXPECT_NEAR(m.std(0), std::sqrt(5.0 / 6.0), 0.05);
  EXPECT_NEAR(m.std(1), std::sqrt(0.5), 0.04);
}

TEST(Stats, AutocorrelationTimeOfAnAr1SeriesAndRhatOfDistinctChains) {
  // AR(1) with phi = 0.8: tau_int = (1+phi)/(1-phi) = 9
  Rng r(5);
  const int n = 200000;
  Vec x(n);
  x(0) = 0;
  for (int t = 1; t < n; ++t) x(t) = 0.8 * x(t - 1) + r.normal();
  EXPECT_NEAR(integrated_autocorrelation_time(x), 9.0, 0.9);
  // chains with different means: R-hat >> 1
  std::vector<Vec> cs(2, Vec(1000));
  Rng q(9);
  for (int i = 0; i < 1000; ++i) {
    cs[0](i) = q.normal();
    cs[1](i) = 5.0 + q.normal();
  }
  EXPECT_GT(gelman_rubin(cs), 2.0);
  // identical-distribution chains: ~1
  for (int i = 0; i < 1000; ++i) cs[1](i) = q.normal();
  EXPECT_LT(gelman_rubin(cs), 1.02);
}

TEST(Stats, PropagationAndAitken) {
  Mat C(2, 2);
  C << 5.0 / 6.0, -0.5, -0.5, 0.5;
  EXPECT_NEAR(propagate_variance((Vec(2) << 1, 3).finished(), C), 7.0 / 3.0,
              1e-14);  // toy 2: sigma_P^2
  EXPECT_NEAR(propagate_covariance((Vec(2) << 1, 0).finished(),
                                   (Vec(2) << 0, 1).finished(), C),
              -0.5, 1e-15);
  // geometric sequence P_k = P_inf + A rho^k  ->  exact limit from three terms
  const double Pinf = 1.2345, A = 0.7, rho = 0.31;
  EXPECT_NEAR(aitken_limit(Pinf + A, Pinf + A * rho, Pinf + A * rho * rho),
              Pinf, 1e-13);
}

TEST(Stats, ErrorBudgetKeepsAmbiguityAndModelSpreadSeparate) {
  Mat C(2, 2);
  C << 5.0 / 6.0, -0.5, -0.5, 0.5;
  const Vec g = (Vec(2) << 1, 3).finished();  // sigma_par^2 = 7/3 (toy 2)
  const Budget b =
      error_budget(5.33, g, C, 0.004, 0.0083, {11.6, 12.9}, {3.79, 3.81, 3.83});
  EXPECT_NEAR(b.sigma_par, std::sqrt(7.0 / 3.0), 1e-12);
  EXPECT_NEAR(
      b.sigma_amb, 0.65,
      1e-12);  // half of 12.9 - 11.6 (the notes' Gamma_ee(J/psi) example)
  EXPECT_NEAR(b.model_spread, 0.02, 1e-12);  // not part of total()
  EXPECT_NEAR(b.total(), std::sqrt(7.0 / 3.0 + 0.004 * 0.004 + 0.0083 * 0.0083),
              1e-12);
  // no variants: no ambiguity
  EXPECT_EQ(error_budget(1, g, C, 0, 0).sigma_amb, 0.0);
}
