// M5: least squares, theory floor, PDG data, and the full worked example of the
// notes (Test 6).
#include <gtest/gtest.h>

#include <cmath>
#include <filesystem>

#include "data/dataset.hpp"
#include "fit/chi2_problem.hpp"
#include "fit/driver.hpp"
#include "fit/joint.hpp"
#include "fit/least_squares.hpp"
#include "models/nr_model.hpp"

using namespace qgem;

namespace {
// Straight line mu = a + b x through (0,1), (1,2), (2,4), sigma = 1 (notes Toy
// 2, 9.6).
class StraightLine final : public LeastSquaresProblem {
 public:
  int n_data() const override { return 3; }
  int n_free() const override { return 2; }
  Result<Vec> residuals(const Vec& u) override {
    Vec z(3);
    const double x[] = {0, 1, 2}, y[] = {1, 2, 4};
    for (int i = 0; i < 3; ++i) z(i) = y[i] - (u(0) + u(1) * x[i]);
    return z;
  }
  Result<Mat> jacobian(const Vec&) override {
    Mat J(3, 2);
    for (int i = 0; i < 3; ++i) J.row(i) << -1.0, -static_cast<double>(i);
    return J;
  }
};
std::string data_path() {
  return std::string(QGEM_SOURCE_DIR) + "/data/pdg-2026.0.sqlite";
}
}  // namespace

TEST(LevenbergMarquardt, StraightLineToyOfTheNotes) {
  StraightLine prob;
  ParamSet ps({{"a", 0.0, true}, {"b", 0.0, true}});
  const FitResult f = levenberg_marquardt(prob, Vec::Zero(2));
  ASSERT_TRUE(f.converged);
  EXPECT_NEAR(f.u(0), 5.0 / 6.0, 1e-9);
  EXPECT_NEAR(f.u(1), 3.0 / 2.0, 1e-9);
  EXPECT_NEAR(f.chi2, 1.0 / 6.0, 1e-12);
  const Uncertainty U = analyse(f, ps);
  EXPECT_NEAR(U.sigma(0), std::sqrt(5.0 / 6.0), 1e-9);            // 0.912871
  EXPECT_NEAR(U.sigma(1), std::sqrt(0.5), 1e-9);                  // 0.707107
  EXPECT_NEAR(U.rho(0, 1), -3.0 / std::sqrt(15.0), 1e-9);         // -0.774597
  EXPECT_NEAR(U.sigma_conditional(0), 1 / std::sqrt(3.0), 1e-9);  // 0.577
  EXPECT_NEAR(U.sigma_conditional(1), 1 / std::sqrt(5.0), 1e-9);  // 0.447
  // prediction at x = 3: P = a + 3 b = 16/3, sigma_P^2 = g^T C g = 7/3
  const Vec g = (Vec(2) << 1.0, 3.0).finished();
  EXPECT_NEAR(f.u(0) + 3 * f.u(1), 16.0 / 3.0, 1e-9);
  EXPECT_NEAR(g.dot(U.Ctheta * g), 7.0 / 3.0, 1e-9);
}

TEST(LevenbergMarquardt, PositivityByLogTransformNeverLeavesTheBound) {
  // minimize (x - 2)^2 + (y - 0.01)^2 style problem with u = log(theta):
  // solution stays > 0.
  class Prob final : public LeastSquaresProblem {
   public:
    int n_data() const override { return 2; }
    int n_free() const override { return 1; }
    Result<Vec> residuals(const Vec& u) override {
      return (Vec(2) << std::exp(u(0)) - 3.0, 0.1 * (std::exp(u(0)) - 5.0))
          .finished();
    }
    Result<Mat> jacobian(const Vec& u) override {
      return (Mat(2, 1) << std::exp(u(0)), 0.1 * std::exp(u(0))).finished();
    }
  } prob;
  const FitResult f =
      levenberg_marquardt(prob, Vec::Constant(1, std::log(100.0)));
  ASSERT_TRUE(f.converged);
  EXPECT_NEAR(std::exp(f.u(0)), (3.0 + 0.01 * 5.0) / (1.0 + 0.01), 1e-9);
}

TEST(TheoryFloor, SolvesChiSquarePerDofEqualsOne) {
  Vec Delta(6), s(6);
  Delta << 0.004, -0.009, 0.017, 0.002, -0.011, 0.006;
  s << 0.0003, 0.00001, 0.0002, 0.0001, 0.0006, 0.00005;
  const int nu = 3;
  const Mat V = s.array().square().matrix().asDiagonal();
  const double sth = solve_theory_floor(Delta, V, nu);
  double chi2 = 0;
  for (int i = 0; i < 6; ++i)
    chi2 += Delta(i) * Delta(i) / (s(i) * s(i) + sth * sth);
  EXPECT_NEAR(chi2, nu, 1e-9);
  // no floor needed when the experimental errors already explain the residuals
  EXPECT_EQ(solve_theory_floor(Delta, Mat(1e6 * V), nu), 0.0);
}

TEST(TheoryFloor, UsesTheFullExperimentalCovariance) {
  // a correlated pair (common systematic): chi2 = Delta^T (V + x I)^-1 Delta
  Vec D(3);
  D << 0.003, 0.004, -0.002;
  Mat V = Mat::Zero(3, 3);
  V.diagonal() << 1e-6, 1e-6, 4e-6;
  V(0, 1) = V(1, 0) = 0.8e-6;
  const double sth = solve_theory_floor(D, V, 2.0);
  Mat C = V;
  C.diagonal().array() += sth * sth;
  EXPECT_NEAR(D.dot(C.llt().solve(D)), 2.0, 1e-9);
  // the diagonal approximation gives a different floor
  const Mat Vd = V.diagonal().asDiagonal();
  EXPECT_GT(std::fabs(solve_theory_floor(D, Vd, 2.0) - sth), 1e-5);
}

TEST(TheoryFloor, ConfidenceIntervalInvertsTheChiSquareTest) {
  for (int nu : {1, 4, 9}) {
    EXPECT_NEAR(1 - chi2_sf(chi2_quantile(0.16, nu), nu), 0.16, 1e-10);
    EXPECT_NEAR(1 - chi2_sf(chi2_quantile(0.975, nu), nu), 0.975, 1e-10);
  }
  // t coverage: P(|T| < 2) = 0.860674 for nu = 3; -> 0.9545 as nu grows
  EXPECT_NEAR(t_coverage(2.0, 3), 0.860674031441, 1e-9);
  EXPECT_NEAR(t_coverage(1.0, 8), 0.653406492913, 1e-9);
  EXPECT_NEAR(t_coverage(2.0, 100000), 0.954500, 1e-5);
  // the floor's interval contains the estimate and widens with the level
  Vec D(6);
  D << 0.004, -0.009, 0.017, 0.002, -0.011, 0.006;
  const Mat V = Mat(Vec::Constant(6, 1e-8).asDiagonal());
  const int nu = 3;
  const double est = solve_theory_floor(D, V, nu);
  const double lo68 = solve_theory_floor(D, V, chi2_quantile(0.84135, nu));
  const double hi68 = solve_theory_floor(D, V, chi2_quantile(0.15865, nu));
  const double lo95 = solve_theory_floor(D, V, chi2_quantile(0.97725, nu));
  const double hi95 = solve_theory_floor(D, V, chi2_quantile(0.02275, nu));
  EXPECT_LT(lo95, lo68);
  EXPECT_LT(lo68, est);
  EXPECT_LT(est, hi68);
  EXPECT_LT(hi68, hi95);
  // with V = s^2 I, chi2 = |D|^2 / (s^2 + sigma^2) exactly
  EXPECT_NEAR(lo68,
              std::sqrt(D.squaredNorm() / chi2_quantile(0.84135, nu) - 1e-8),
              1e-9 * lo68);
}

TEST(ParamSet, TransformsTiesAndChainRule) {
  ParamSet ps({{"m1", 1.5, true, Transform::Log},
               {"m2", 0.0, false, Transform::Identity, -1e300, 1e300, "m1"},
               {"b", 0.18, true, Transform::Log},
               {"c", 0.3, false}});
  const Vec u = ps.u0();
  ASSERT_EQ(u.size(), 2);
  const Vec th = ps.theta(u);
  EXPECT_NEAR(th(0), 1.5, 1e-15);
  EXPECT_NEAR(th(1), 1.5, 1e-15);  // tied
  EXPECT_NEAR(th(3), 0.3, 1e-15);
  const Mat K = ps.dtheta_du(u);
  EXPECT_NEAR(K(0, 0), 1.5, 1e-15);
  EXPECT_NEAR(K(1, 0), 1.5, 1e-15);  // the tied parameter moves with its source
  EXPECT_NEAR(K(2, 1), 0.18, 1e-15);
  EXPECT_NEAR(K(3, 0), 0.0, 1e-15);
}

TEST(Pdg, ReaderGivesTheNotesCharmoniumTable) {
  if (!std::filesystem::exists(data_path()))
    GTEST_SKIP() << "PDG snapshot not present (data/ is git-ignored)";
  PdgDatabase pdg(data_path());
  auto eta = pdg.mass(441), jpsi = pdg.mass(443), chi0 = pdg.mass(10441),
       psi2s = pdg.mass(100443);
  ASSERT_TRUE(eta && jpsi && chi0 && psi2s);
  EXPECT_NEAR(eta->value, 2.98409, 5e-5);  // 2984.09 +- 0.30 MeV
  EXPECT_NEAR(eta->err, 0.0003, 2e-5);
  EXPECT_NEAR(jpsi->value, 3.0969, 1e-4);
  EXPECT_NEAR(jpsi->err, 6.2e-6, 1e-6);
  EXPECT_NEAR(chi0->value, 3.4155, 1e-4);
  EXPECT_NEAR(psi2s->value, 3.686097,
              1e-6);  // the fit value (in the summary table), not the average
  EXPECT_NEAR(psi2s->err, 1.055e-5, 1e-7);
  // chi_b(1P): the second (common) error is only in the display text
  auto chib1 = pdg.mass(20553);
  ASSERT_TRUE(chib1);
  EXPECT_NEAR(chib1->err, 0.00026, 1e-8);
  EXPECT_NEAR(chib1->common_err, 0.00031, 1e-8);
  EXPECT_FALSE(pdg.mass(999999999).has_value());
  // a partial width the PDG gives only as BR x total width (chi_c0 -> gamma
  // gamma: 1.96e-4 x 12.40 MeV)
  const auto gg0 = pdg.partial_width(10441, "gamma gamma");
  ASSERT_TRUE(gg0);
  EXPECT_TRUE(gg0->derived);
  EXPECT_NEAR(gg0->value_keV, 2.433, 0.005);
  EXPECT_FALSE(pdg.partial_width(443, "e+ e-")->derived);  // listed directly
  EXPECT_FALSE(pdg.partial_width(551, "gamma gamma").has_value());  // not seen
  EXPECT_DOUBLE_EQ(second_error_from_display("9859.44+-0.42+-0.31"), 0.31);
  EXPECT_DOUBLE_EQ(second_error_from_display("3096.900+-0.006"), 0.0);
}

TEST(Dataset, CommonSystematicsCorrelateTheGroup) {
  if (!std::filesystem::exists(data_path()))
    GTEST_SKIP() << "PDG snapshot not present";
  PdgDatabase pdg(data_path());
  const std::string f =
      std::string(QGEM_SOURCE_DIR) + "/config/bbar_fit.states";
  const Dataset ds = load_dataset(f, pdg);
  const Mat V = ds.V_exp();
  ASSERT_GT(V.rows(), 3);
  // chi_bJ(1P) share a 0.31 MeV common error: off-diagonal V = 0.31^2 MeV^2,
  // diagonal err^2 + 0.31^2
  int i0 = -1, i1 = -1, k = 0;
  for (const Datum* d : ds.fitted()) {
    if (d->mcid == 10551) i0 = k;
    if (d->mcid == 20553) i1 = k;
    ++k;
  }
  ASSERT_GE(i0, 0);
  ASSERT_GE(i1, 0);
  EXPECT_NEAR(V(i0, i1), 0.00031 * 0.00031, 1e-14);
  EXPECT_NEAR(V(i0, i0), 0.00042 * 0.00042 + 0.00031 * 0.00031, 1e-14);
}

// ---- notes Test 6: Model A fitted to the eight charmonium states below
// threshold (10.10) ---------------------
TEST(Test6, CharmoniumWorkedExampleOfTheNotes) {
  if (!std::filesystem::exists(data_path()))
    GTEST_SKIP() << "PDG snapshot not present (data/ is git-ignored)";
  PdgDatabase pdg(data_path());
  Dataset ds = load_dataset(
      std::string(QGEM_SOURCE_DIR) + "/config/ccbar_fit.states", pdg);
  NRModel model({.screened = false,
                 .spin = SpinMode::FirstOrder,
                 .contact_in_all_L = true});
  ParamSet ps({{"m1", 1.48, true, Transform::Log},
               {"m2", 1.48, false, Transform::Identity, -1e300, 1e300, "m1"},
               {"alpha_s", 0.55, true, Transform::Log},
               {"b", 0.14, true, Transform::Log},
               {"c", 0.0, false},
               {"sigma", 1.1, true, Transform::Log},
               {"mu_s", 0.1, false}});
  Chi2Problem prob(
      model, ds, ps,
      MatrixCache(BasisSpec{.n = 30, .r1 = 0.02, .rmax = 30.0, .eta = 1e-12}),
      0.010);
  Vec u = ps.u0();
  double sth = 0.010;
  FitResult f;
  for (int it = 0; it < 20; ++it) {
    prob.set_sigma_th(sth);
    f = levenberg_marquardt(prob, u);
    u = f.u;
    const double next =
        solve_theory_floor(*prob.raw_residuals(u), prob.exp_covariance(),
                           prob.n_data() - prob.n_free());
    if (std::fabs(next - sth) < 1e-8) break;
    sth = next;
  }
  prob.set_sigma_th(sth);
  f = levenberg_marquardt(prob, u);
  ASSERT_TRUE(f.converged);
  EXPECT_NEAR(1e3 * sth, 8.34, 0.01);  // notes: sigma_th = 8.34 MeV
  EXPECT_NEAR(f.chi2, 4.0, 1e-6);  // chi^2 = nu by construction of the floor
  const Vec th = ps.theta(f.u);
  EXPECT_NEAR(th(kM1), 1.4670, 2e-4);
  EXPECT_NEAR(th(kAlphaS), 0.5339, 2e-4);
  EXPECT_NEAR(th(kB), 0.1486, 2e-4);
  EXPECT_NEAR(th(kSigma), 1.130, 2e-3);
  const Uncertainty U = analyse(f, ps);
  EXPECT_NEAR(U.sigma(0), 0.0209, 2e-4);
  EXPECT_NEAR(U.sigma(1), 0.0242, 3e-4);
  EXPECT_NEAR(U.sigma(2), 0.0071, 2e-4);
  EXPECT_NEAR(U.sigma(3), 0.107, 2e-3);
  // correlation matrix (m_c, alpha_s, b, sigma)
  EXPECT_NEAR(U.rho(0, 1), 0.939, 2e-3);
  EXPECT_NEAR(U.rho(0, 2), -0.969, 2e-3);
  EXPECT_NEAR(U.rho(0, 3), -0.466, 2e-3);
  EXPECT_NEAR(U.rho(1, 2), -0.839, 2e-3);
  EXPECT_NEAR(U.rho(1, 3), -0.626, 2e-3);
  EXPECT_NEAR(U.rho(2, 3), 0.363, 2e-3);
  // relative Fisher eigenvalues (94.5, 364, 4.19e3, 7.88e5)
  EXPECT_NEAR(U.fisher_eigenvalues(0), 94.5, 0.5);
  EXPECT_NEAR(U.fisher_eigenvalues(1), 364.0, 2.0);
  EXPECT_NEAR(U.fisher_eigenvalues(3) / 7.88e5, 1.0, 5e-3);
  // residuals (MeV): +2.1, +4.5, -13.8, +0.9, +5.9, -3.0, -0.4, +4.1
  const Vec r = *prob.raw_residuals(f.u);
  const double expect[] = {2.1, 4.5, -13.8, 0.9, 5.9, -3.0, -0.4, 4.1};
  for (int i = 0; i < 8; ++i) EXPECT_NEAR(1e3 * r(i), expect[i], 0.15) << i;
  // numerical error by nested doubling n = 30 -> 59: masses move by <= 6e-5 MeV
  // (notes 10.10), far below 1e-2 min s_i (= 1e-2 * 8.3 MeV) required of the
  // numerics
  const Vec snum = numerical_errors(
      prob, th, BasisSpec{.n = 30, .r1 = 0.02, .rmax = 30.0, .eta = 1e-12});
  for (int i = 0; i < 8; ++i) {  // the eight fitted states; high-n predictions
                                 // need a larger basis (see qgem report)
    EXPECT_LT(1e3 * snum(i), 1e-3) << i;
    EXPECT_LT(snum(i), 1e-2 * 0.5 * 8.3e-3);
  }
  // Jacobian by Feynman-Hellmann against central differences of the whitened
  // residuals
  const Mat J = *prob.jacobian(f.u);
  for (int a = 0; a < 4; ++a) {
    Vec up = f.u, um = f.u;
    const double h = 1e-5;
    up(a) += h;
    um(a) -= h;
    const Vec fd = (*prob.residuals(up) - *prob.residuals(um)) / (2 * h);
    EXPECT_LT((J.col(a) - fd).cwiseAbs().maxCoeff(),
              5e-5 * std::max(1.0, fd.cwiseAbs().maxCoeff()))
        << a;
  }
  // Held-out refit with the floor re-solved: chi2/nu = 1 on the kept data, so
  // the held-out state does not set its own predictive scale. Removing the
  // largest residual (chi_c0, -13.8 MeV) lowers the floor.
  const HeldOutFit h = refit_held_out(prob, f.u, {2}, {}, true);
  ASSERT_EQ(h.mu.size(), ds.d.size());
  EXPECT_LT(h.sigma_th, 0.8 * sth);
  EXPECT_GT(h.sigma_par(2), 0.0);
  Dataset reduced = ds;
  reduced.d[2].in_fit = false;
  Chi2Problem sub(
      model, reduced, ps,
      MatrixCache(BasisSpec{.n = 30, .r1 = 0.02, .rmax = 30.0, .eta = 1e-12}),
      h.sigma_th);
  EXPECT_NEAR(sub.residuals(ps.u(h.theta))->squaredNorm(), 3.0, 1e-6);
  // with nu = 0 left the floor cannot be solved
  EXPECT_TRUE(
      std::isnan(refit_held_out(prob, f.u, {0, 1, 2, 3}, {}, true).sigma_th));
  // Closure: a shift along dmu/db (1 MeV at most) is absorbed by b alone, so
  // nothing is left (R = O(shift^2)) and the a-priori floor vanishes.
  const auto pgr = prob.predictions_with_gradients(th);
  ASSERT_TRUE(pgr.has_value());
  const Vec g = pgr->second.col(kB);
  const double eps = 1e-3 / g.cwiseAbs().maxCoeff();
  const Truncation tc = absorb_shift(prob, f.u, Vec(eps * g));
  ASSERT_EQ(tc.R.size(), static_cast<int>(ds.d.size()));
  EXPECT_LT(tc.R.cwiseAbs().maxCoeff(), 1e-5);  // < 0.01 MeV
  EXPECT_NEAR(tc.theta(kB), th(kB) - eps, 1e-3 * eps);
  EXPECT_LT(tc.sigma_a, 1e-5);
}

TEST(Chi2Problem, CollapsedSpectrumIsOutsideTheModel) {
  // A point-like contact term (sigma -> large) makes H unbounded below in the
  // 1S0 channel: the lowest level is set by the basis, not the physics. The
  // chi^2 must refuse such points (a wall for the minimiser), while level
  // tracking would otherwise follow the "physical" level past the collapse.
  if (!std::filesystem::exists(data_path())) GTEST_SKIP();
  PdgDatabase pdg(data_path());
  Dataset ds = load_dataset(
      std::string(QGEM_SOURCE_DIR) + "/config/ccbar_fit.states", pdg);
  NRModel m;
  ParamSet ps = m.default_params(1.48, 1.48, 0.55, 0.14, 0.0, 1.1);
  auto sp = ps.specs();
  sp[1].tied_to = "m1";
  ps = ParamSet(sp).with_free({"m1", "alpha_s", "b", "sigma"});
  Chi2Problem prob(
      m, ds, ps, MatrixCache({.n = 36, .r1 = 0.02, .rmax = 60.0, .eta = 1e-12}),
      0.01);
  Vec th = ps.theta0();
  EXPECT_TRUE(prob.residuals(ps.u(th)).has_value());
  th(ps.index("sigma")) = 1000.0;
  const auto z = prob.residuals(ps.u(th));
  ASSERT_FALSE(z.has_value());
  EXPECT_NE(z.error().message.find("collapsed"), std::string::npos);
}

TEST(Jackknife, DeleteOneFormulaAndFailedRows) {
  // x = 1,2,3,4: mean 2.5, sum of squares 5, sigma = sqrt(3/4 * 5)
  Mat s(5, 2);
  s << 1, 10, 2, 10, 3, 10, 4, 10, std::nan(""), 0;  // last row: failed refit
  const Vec j = jackknife_sigma(s);
  EXPECT_NEAR(j(0), std::sqrt(0.75 * 5.0), 1e-14);
  EXPECT_NEAR(j(1), 0.0, 1e-14);
  EXPECT_TRUE(std::isnan(jackknife_sigma(Mat::Constant(1, 1, 2.0))(0)));
}

TEST(Chi2Problem, UnresolvedSmearedContactIsOutsideTheModel) {
  // A smeared contact inside H narrower than the narrowest Gaussian
  // (sigma * r1 > 1) has an energy set by the basis: refused.
  if (!std::filesystem::exists(data_path())) GTEST_SKIP();
  PdgDatabase pdg(data_path());
  Dataset ds = load_dataset(
      std::string(QGEM_SOURCE_DIR) + "/config/ccbar_fit.states", pdg);
  NRModel m({.screened = false, .spin = SpinMode::Smeared});
  ParamSet ps = m.default_params(1.6, 1.6, 0.9, 0.12, 0.0, 1.5);
  auto sp = ps.specs();
  sp[1].tied_to = "m1";
  ps = ParamSet(sp).with_free({"m1", "alpha_s", "b", "sigma"});
  Chi2Problem prob(
      m, ds, ps, MatrixCache({.n = 36, .r1 = 0.05, .rmax = 15.0, .eta = 1e-12}),
      0.01);
  Vec th = ps.theta0();
  EXPECT_TRUE(prob.residuals(ps.u(th)).has_value());  // sigma r1 = 0.075
  th(ps.index("sigma")) = 30.0;                       // sigma r1 = 1.5
  const auto z = prob.residuals(ps.u(th));
  ASSERT_FALSE(z.has_value());
  EXPECT_NE(z.error().message.find("resolution"), std::string::npos);
}

TEST(TheoryFloor, StateDependentWeightsSolveTheWeightedEquation) {
  // sum_i Delta_i^2 / (s_i^2 + w_i^2 sigma^2) = nu
  const Vec D = (Vec(4) << 0.010, -0.004, 0.006, 0.020).finished();
  const Vec s = Vec::Constant(4, 0.001);
  const Mat V = s.array().square().matrix().asDiagonal();
  const Vec w = (Vec(4) << 1.0, 0.5, 1.0, 2.0).finished();
  const int nu = 2;
  const double sig = solve_theory_floor(D, V, nu, w);
  double chi2 = 0;
  for (int i = 0; i < 4; ++i)
    chi2 += D(i) * D(i) / (s(i) * s(i) + w(i) * w(i) * sig * sig);
  EXPECT_NEAR(chi2, nu, 1e-9);
  // uniform weights reproduce the original floor
  EXPECT_NEAR(solve_theory_floor(D, V, nu, Vec::Ones(4)),
              solve_theory_floor(D, V, nu), 1e-15);
}

TEST(Joint, SharedParametersJacobianAndSeparateLimit) {
  if (!std::filesystem::exists(data_path()))
    GTEST_SKIP() << "PDG snapshot not present (data/ is git-ignored)";
  PdgDatabase pdg(data_path());
  Dataset ds = load_dataset(
      std::string(QGEM_SOURCE_DIR) + "/config/ccbar_fit.states", pdg);
  NRModel model({.screened = false,
                 .spin = SpinMode::FirstOrder,
                 .contact_in_all_L = true});
  ParamSet ps({{"m1", 1.48, true, Transform::Log},
               {"m2", 1.48, false, Transform::Identity, -1e300, 1e300, "m1"},
               {"alpha_s", 0.55, true, Transform::Log},
               {"b", 0.14, true, Transform::Log},
               {"c", 0.0, false},
               {"sigma", 1.1, true, Transform::Log},
               {"mu_s", 0.1, false}});
  const BasisSpec bs{.n = 30, .r1 = 0.02, .rmax = 30.0, .eta = 1e-12};
  Chi2Problem one(model, ds, ps, MatrixCache(bs), 0.010);
  const FloorFit sep = fit_with_floor(one, ps.u0());
  ASSERT_TRUE(sep.fit.converged);
  // two copies of the same system: nothing shared reproduces the separate
  // fit in each, with nu_s = N - p
  Chi2Problem a(model, ds, ps, MatrixCache(bs), 0.010),
      b(model, ds, ps, MatrixCache(bs), 0.010);
  JointProblem none(a, b, {});
  EXPECT_EQ(none.n_free(), 8);
  const JointFit j0 =
      fit_joint(none, none.join(sep.fit.u, sep.fit.u), 0.010, 0.010);
  EXPECT_NEAR(j0.sigma_a, sep.sigma_th, 1e-7);
  EXPECT_NEAR(j0.sigma_b, sep.sigma_th, 1e-7);
  EXPECT_NEAR(j0.nu_a, 4.0, 1e-6);
  EXPECT_LT((none.u_a(j0.fit.u) - sep.fit.u).cwiseAbs().maxCoeff(), 1e-5);
  // sharing b: one parameter fewer; the shared column is the sum of both
  // systems' columns (finite differences)
  JointProblem jb(a, b, {"b"});
  EXPECT_EQ(jb.n_free(), 7);
  EXPECT_EQ(jb.k(), 1);
  EXPECT_EQ(jb.names().back(), "sigma(b)");
  a.set_sigma_th(sep.sigma_th), b.set_sigma_th(sep.sigma_th);
  const Vec u = jb.join(sep.fit.u, sep.fit.u);
  const Mat J = *jb.jacobian(u);
  for (int c = 0; c < jb.n_free(); ++c) {
    Vec up = u, um = u;
    up(c) += 1e-5, um(c) -= 1e-5;
    const Vec fd = (*jb.residuals(up) - *jb.residuals(um)) / 2e-5;
    EXPECT_LT((J.col(c) - fd).cwiseAbs().maxCoeff(),
              5e-5 * std::max(1.0, fd.cwiseAbs().maxCoeff()))
        << c;
  }
  // identical systems agree on b: the constraint costs nothing
  const FitResult fs = levenberg_marquardt(jb, u);
  EXPECT_NEAR(fs.chi2, 2 * 4.0, 1e-6);
  EXPECT_THROW(JointProblem(a, b, {"c"}), std::invalid_argument);  // not free
}
