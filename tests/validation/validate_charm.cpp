// Slow validation (ctest label "slow"): the statistics of the notes' worked
// example (10.10) on the physics fit: leverages, leave-one-out with exact
// refits, profile of m_c, the m_c-c degeneracy, AIC/BIC, bootstrap, MCMC.
#include <gtest/gtest.h>

#include <cmath>
#include <filesystem>

#include "data/dataset.hpp"
#include "fit/driver.hpp"
#include "models/nr_model.hpp"

using namespace qgem;

namespace {
struct Charm {
  PdgDatabase pdg{std::string(QGEM_SOURCE_DIR) + "/data/pdg-2026.0.sqlite"};
  NRModel model{{.screened = false,
                 .spin = SpinMode::FirstOrder,
                 .contact_in_all_L = true}};
  ParamSet ps{std::vector<ParamSpec>{
      {"m1", 1.48, true, Transform::Log, -1e300, 1e300, ""},
      {"m2", 1.48, false, Transform::Identity, -1e300, 1e300, "m1"},
      {"alpha_s", 0.55, true, Transform::Log, -1e300, 1e300, ""},
      {"b", 0.14, true, Transform::Log, -1e300, 1e300, ""},
      {"c", 0.0, false, Transform::Identity, -1e300, 1e300, ""},
      {"sigma", 1.1, true, Transform::Log, -1e300, 1e300, ""},
      {"mu_s", 0.1, false, Transform::Identity, -1e300, 1e300, ""}}};
  Dataset ds = load_dataset(
      std::string(QGEM_SOURCE_DIR) + "/config/ccbar_fit.states", pdg);
  Chi2Problem prob{
      model, ds, ps,
      MatrixCache(BasisSpec{.n = 30, .r1 = 0.02, .rmax = 30.0, .eta = 1e-12}),
      0.01};
  FloorFit ff = fit_with_floor(prob, ps.u0());
};
Charm& charm() {
  static Charm c;
  return c;
}
bool have_data() {
  return std::filesystem::exists(std::string(QGEM_SOURCE_DIR) +
                                 "/data/pdg-2026.0.sqlite");
}
}  // namespace

TEST(Validation, LeverageAndLeaveOneOutAgainstTheNotesTable) {
  if (!have_data()) GTEST_SKIP();
  Charm& c = charm();
  const LooResult l = leave_one_out(c.ff.fit);
  const double h[] = {0.895, 0.666, 0.231, 0.182, 0.259, 0.596, 0.530, 0.643};
  const double pull[] = {0.75, 0.96, -1.89, 0.12, 0.81, -0.54, -0.06, 0.77};
  for (int i = 0; i < 8; ++i) {
    EXPECT_NEAR(l.leverage(i), h[i], 2e-3) << i;
    EXPECT_NEAR(l.pull(i), pull[i], 0.08)
        << i;  // 1-h_i as small as 0.105 amplifies tiny data differences
  }
  EXPECT_NEAR(l.leverage.sum(), 4.0, 1e-8);  // sum h_i = p
  // exact refits: (y_i - mu_i^(-i)) in MeV: +19.4, +14.0, -18.1, +1.1, +7.9,
  // -7.0, -0.8, +10.9
  const Vec exact = exact_leave_one_out(c.prob, c.ff.fit.u);
  const double ex[] = {19.4, 14.0, -18.1, 1.1, 7.9, -7.0, -0.8, 10.9};
  for (int i = 0; i < 8; ++i) EXPECT_NEAR(1e3 * exact(i), ex[i], 0.3) << i;
  // the one-fit formula Delta_i / (1 - h_i) reproduces the refits to < 0.7 MeV
  // (notes)
  const Vec raw = *c.prob.raw_residuals(c.ff.fit.u);
  for (int i = 0; i < 8; ++i)
    EXPECT_NEAR(1e3 * raw(i) / (1 - l.leverage(i)), 1e3 * exact(i), 0.8) << i;
}

TEST(Validation, ProfileOfCharmMassIsNearlyParabolic) {
  if (!have_data()) GTEST_SKIP();
  Charm& c = charm();
  const Uncertainty U = analyse(c.ff.fit, c.ps);
  const double m = c.ps.theta(c.ff.fit.u)(kM1), sm = U.sigma(0);
  std::vector<double> grid;
  for (double t : {-2.0, -1.0, 1.0, 2.0}) grid.push_back(std::log(m + t * sm));
  const auto pr = profile(c.prob, c.ff.fit, 0, grid);
  // notes: 1.086, 1.064 at -+1 sigma and 4.39, 4.22 at -+2 sigma (parabola: 1
  // and 4)
  EXPECT_NEAR(pr[1].second, 1.086, 0.03);
  EXPECT_NEAR(pr[2].second, 1.064, 0.03);
  EXPECT_NEAR(pr[0].second, 4.39, 0.1);
  EXPECT_NEAR(pr[3].second, 4.22, 0.1);
}

TEST(Validation, FreeingTheConstantReproducesTheMcCDegeneracy) {
  if (!have_data()) GTEST_SKIP();
  Charm& c = charm();
  ParamSet ps5 = c.ps.with_values(c.ps.theta(c.ff.fit.u))
                     .with_free({"m1", "alpha_s", "b", "sigma", "c"});
  Chi2Problem p5(
      c.model, c.ds, ps5,
      MatrixCache(BasisSpec{.n = 30, .r1 = 0.02, .rmax = 30.0, .eta = 1e-12}),
      c.ff.sigma_th);
  Vec u0 = ps5.u0();
  const FloorFit f5 = fit_with_floor(
      p5, u0, c.ff.sigma_th, {}, c.ff.sigma_th);  // same floor for both fits
  ASSERT_TRUE(f5.fit.converged);
  const Vec th = ps5.theta(f5.fit.u);
  const Uncertainty U = analyse(f5.fit, ps5);
  // notes: m_c = 1.26 +- 0.13, c = 0.39 +- 0.23, rho(m_c, c) = -0.987, chi^2
  // = 1.67, dchi2 = 2.33 (p = 0.13)
  EXPECT_NEAR(f5.fit.chi2, 1.67, 0.05);
  EXPECT_NEAR(c.ff.fit.chi2 - f5.fit.chi2, 2.33, 0.05);
  EXPECT_NEAR(th(kM1), 1.26, 0.03);
  EXPECT_NEAR(th(kC), 0.39, 0.03);
  EXPECT_NEAR(U.sigma(0), 0.13, 0.015);
  // free parameters are ordered (m1, alpha_s, b, c, sigma): c has index 3
  EXPECT_NEAR(U.sigma(3), 0.23, 0.02);
  EXPECT_NEAR(U.rho(0, 3), -0.987, 0.005);
  EXPECT_GT(U.condition, 1e3);  // near-degenerate direction
  EXPECT_NEAR(nested_test_pvalue(c.ff.fit.chi2, f5.fit.chi2, 1), 0.13, 0.01);
  // AIC 12.00 -> 11.67, BIC 12.32 -> 12.07: the data do not ask for c
  EXPECT_NEAR(aic(c.ff.fit.chi2, 4), 12.00, 1e-4);
  EXPECT_NEAR(aic(f5.fit.chi2, 5), 11.67, 0.05);
  EXPECT_NEAR(bic(c.ff.fit.chi2, 4, 8), 12.32, 0.01);
  EXPECT_NEAR(bic(f5.fit.chi2, 5, 8), 12.07, 0.05);
}

TEST(Validation, BootstrapAndPosteriorAgreeWithTheLinearErrors) {
  if (!have_data()) GTEST_SKIP();
  Charm& c = charm();
  const Uncertainty U = analyse(c.ff.fit, c.ps);
  // bootstrap (B = 300 in the notes: std 0.0225, 0.0266, 0.0076, 0.118 -- 7-10%
  // above sqrt(C_aa))
  const BootstrapResult b = bootstrap_fit(c.prob, c.ff.fit.u, 400, 20261011, 4);
  EXPECT_LE(b.failed, 4);
  Mat th(b.u.rows(), 4);
  for (int i = 0; i < b.u.rows(); ++i) {
    const Vec t = c.ps.theta(b.u.row(i).transpose());
    th.row(i) << t(kM1), t(kAlphaS), t(kB), t(kSigma);
  }
  const Vec mean = th.colwise().mean();
  const Mat d = th.rowwise() - mean.transpose();
  const Vec sd = (d.array().square().colwise().sum() / (th.rows() - 1)).sqrt();
  for (int a = 0; a < 4; ++a) {
    EXPECT_GT(sd(a) / U.sigma(a), 0.95) << a;
    EXPECT_LT(sd(a) / U.sigma(a), 1.25)
        << a;  // non-linearity widens the sampling distribution slightly
  }
  // MCMC: notes (4 chains x 6000, flat priors, burn-in 1000) mean (1.4648,
  // 0.5308, 0.1492, 1.148), std (0.0203, 0.0238, 0.0070, 0.112)
  const McmcResult m =
      sample_posterior(c.prob, c.ff.fit, c.ps, 4, 4000, 800, 77);
  for (int a = 0; a < 4; ++a) {
    EXPECT_LT(m.rhat(a), 1.02) << a;
    EXPECT_GT(m.n_eff(a), 400) << a;
  }
  EXPECT_NEAR(m.chains[0].acceptance, 0.3, 0.1);
  // posterior moments in theta from the u-samples
  Vec tmean = Vec::Zero(4), tsq = Vec::Zero(4);
  long n = 0;
  for (const Chain& ch : m.chains)
    for (int i = 0; i < ch.samples.rows(); ++i) {
      const Vec t = c.ps.theta(ch.samples.row(i).transpose());
      Vec v(4);
      v << t(kM1), t(kAlphaS), t(kB), t(kSigma);
      tmean += v;
      tsq += v.cwiseProduct(v);
      ++n;
    }
  tmean /= n;
  const Vec tstd = (tsq / n - tmean.cwiseProduct(tmean)).cwiseSqrt();
  const double nm[] = {1.4648, 0.5308, 0.1492, 1.148},
               ns[] = {0.0203, 0.0238, 0.0070, 0.112};
  for (int a = 0; a < 4; ++a) {
    EXPECT_NEAR(tmean(a), nm[a], 0.25 * ns[a] + 1e-9)
        << a;  // posterior mean is shifted by < 0.25 sigma (notes: ~0.1 sigma)
    EXPECT_NEAR(tstd(a), ns[a], 0.15 * ns[a]) << a;
  }
}

TEST(Validation, PriorOnTheConstantTamesTheMcCDegeneracy) {
  if (!have_data()) GTEST_SKIP();
  Charm& c = charm();
  // c free with a Gaussian prior c ~ N(0, 0.1 GeV): notes 10.9.4, the posterior
  // covariance is (F + Lambda)^-1
  ParamSet ps5 = c.ps.with_values(c.ps.theta(c.ff.fit.u))
                     .with_free({"m1", "alpha_s", "b", "sigma", "c"});
  Chi2Problem p5(
      c.model, c.ds, ps5,
      MatrixCache(BasisSpec{.n = 30, .r1 = 0.02, .rmax = 30.0, .eta = 1e-12}),
      c.ff.sigma_th);
  p5.add_prior(kC, 0.0, 0.1);
  EXPECT_EQ(p5.n_data(), 9);
  const FloorFit f5 =
      fit_with_floor(p5, ps5.u0(), c.ff.sigma_th, {}, c.ff.sigma_th);
  ASSERT_TRUE(f5.fit.converged);
  const Vec th = ps5.theta(f5.fit.u);
  const Uncertainty U = analyse(f5.fit, ps5);
  EXPECT_LT(
      std::fabs(th(kC)),
      0.39);  // pulled from the free-fit value 0.39 toward the prior mean 0
  EXPECT_GT(std::fabs(th(kC)), 0.0);
  EXPECT_LT(U.sigma(3), 0.1);  // tighter than the prior itself, far tighter
                               // than the free fit (0.23)
  EXPECT_NEAR(U.rho(0, 3), -0.987,
              0.2);  // the degeneracy is broken: |rho| drops
  EXPECT_LT(std::fabs(U.rho(0, 3)), 0.987);
  // chi^2 includes the prior term exactly
  const Vec z = *p5.residuals(f5.fit.u);
  EXPECT_NEAR(z.squaredNorm(), f5.fit.chi2, 1e-10);
  EXPECT_NEAR(z(8), (0.0 - th(kC)) / 0.1, 1e-12);
}
