// The optional ROOT bridge: Minuit2 must find the same minimum as
// Levenberg-Marquardt, and its Hesse errors must match (J^T J)^-1 where the
// problem is linear (notes Toy 2) and the worked charmonium fit (notes 10.10).
#include <gtest/gtest.h>

#include <cmath>
#include <filesystem>

#include "data/dataset.hpp"
#include "fit/driver.hpp"
#include "models/nr_model.hpp"
#include "root/minuit2_engine.hpp"

using namespace qgem;

namespace {
// y = (1, 2, 4) at x = (0, 1, 2), unit errors: a + b x.
class StraightLine final : public LeastSquaresProblem {
 public:
  int n_data() const override { return 3; }
  int n_free() const override { return 2; }
  Result<Vec> residuals(const Vec& u) override {
    Vec z(3);
    const double y[3] = {1, 2, 4};
    for (int i = 0; i < 3; ++i) z(i) = y[i] - (u(0) + u(1) * i);
    return z;
  }
  Result<Mat> jacobian(const Vec&) override {
    Mat J(3, 2);
    for (int i = 0; i < 3; ++i) J.row(i) << -1.0, -static_cast<double>(i);
    return J;
  }
};
}  // namespace

TEST(Minuit2, StraightLineMatchesTheExactSolution) {
  StraightLine prob;
  const FitResult f = minuit2_minimize(prob, Vec::Zero(2), {});
  ASSERT_TRUE(f.converged);
  EXPECT_NEAR(f.u(0), 5.0 / 6.0, 1e-6);
  EXPECT_NEAR(f.u(1), 3.0 / 2.0, 1e-6);
  const MinuitErrors e = minuit2_errors(prob, f.u, true);
  ASSERT_TRUE(e.valid);
  EXPECT_NEAR(e.cov_u(0, 0), 5.0 / 6.0, 1e-6);  // Hesse = (J^T J)^-1, exact
  EXPECT_NEAR(e.cov_u(1, 1), 0.5, 1e-6);
  // a linear model: Minos is symmetric and equals Hesse
  EXPECT_NEAR(e.minos_hi(0), std::sqrt(5.0 / 6.0), 1e-4);
  EXPECT_NEAR(e.minos_lo(0), -std::sqrt(5.0 / 6.0), 1e-4);
}

TEST(Minuit2, WorkedCharmoniumFitAgreesWithLevenbergMarquardt) {
  const std::string db =
      std::string(QGEM_SOURCE_DIR) + "/data/pdg-2026.0.sqlite";
  if (!std::filesystem::exists(db)) GTEST_SKIP();
  PdgDatabase pdg(db);
  Dataset ds = load_dataset(
      std::string(QGEM_SOURCE_DIR) + "/config/ccbar_fit.states", pdg);
  NRModel m;
  ParamSet ps = m.default_params(1.48, 1.48, 0.55, 0.14, 0.0, 1.1);
  auto sp = ps.specs();
  sp[1].tied_to = "m1";
  ps = ParamSet(sp).with_free({"m1", "alpha_s", "b", "sigma"});
  const BasisSpec spec{.n = 36, .r1 = 0.02, .rmax = 60.0, .eta = 1e-12};
  Chi2Problem lm(m, ds, ps, MatrixCache(spec), 0.00834);
  const FitResult a = levenberg_marquardt(lm, ps.u0());
  Chi2Problem mn(m, ds, ps, MatrixCache(spec), 0.00834);
  const FitResult b = minuit2_minimize(mn, ps.u0(), {});
  ASSERT_TRUE(a.converged);
  ASSERT_TRUE(b.converged);
  EXPECT_NEAR(a.chi2, b.chi2, 1e-6);
  for (int i = 0; i < a.u.size(); ++i) EXPECT_NEAR(a.u(i), b.u(i), 1e-4) << i;
}
