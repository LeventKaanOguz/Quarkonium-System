// M7: basis optimization (Ky Fan), nested doubling, Nelder-Mead.
#include <gtest/gtest.h>

#include <algorithm>
#include <cmath>

#include "data/dataset.hpp"
#include "hyper/basis_opt.hpp"
#include "math/optimize.hpp"
#include "models/nr_model.hpp"
#include "obs/observables.hpp"

using namespace qgem;

namespace {
Vec theta_of(double m1, double m2, double als, double b, double c, double sig,
             double mu = 0.1) {
  Vec t(kNRParams);
  t << m1, m2, als, b, c, sig, mu;
  return t;
}
Channel single(int L, int S, int J) {
  Channel c;
  c.J = J;
  c.comps = {{L, S}};
  return c;
}
}  // namespace

TEST(NelderMead, MinimizesRosenbrockAndAQuadratic) {
  auto rosen = [](const Vec& x) {
    return 100 * std::pow(x(1) - x(0) * x(0), 2) + std::pow(1 - x(0), 2);
  };
  NelderMeadOptions o;
  o.max_iter = 2000;
  o.xtol = 1e-10;
  const auto r = nelder_mead(rosen, (Vec(2) << -1.2, 1.0).finished(), o);
  EXPECT_NEAR(r.x(0), 1.0, 1e-5);
  EXPECT_NEAR(r.x(1), 1.0, 1e-5);
  auto quad = [](const Vec& x) {
    return (x(0) - 3) * (x(0) - 3) + 2 * (x(1) + 1) * (x(1) + 1);
  };
  const auto q = nelder_mead(quad, Vec::Zero(2));
  EXPECT_NEAR(q.x(0), 3.0, 1e-5);
  EXPECT_NEAR(q.x(1), -1.0, 1e-5);
}

TEST(NestedConvergence, MassesDecreaseMonotonicallyAndAitkenExtrapolates) {
  NRModel model({.screened = false, .spin = SpinMode::None});
  const Vec th = theta_of(1.5, 1.5, 0.5, 0.18, 0.0, 1.1);
  const BasisSpec start{.n = 8, .r1 = 0.03, .rmax = 40.0, .eta = 1e-13};
  const Observable mass3S = [](const SolvedChannel& sc, MatrixCache&) {
    return sc.level(2).M;
  };
  const NestedSequence s =
      nested_convergence(model, th, single(0, 1, 1), start, mass3S, 3);
  ASSERT_EQ(s.n.size(), 4u);
  EXPECT_EQ(s.n[1], 15);
  EXPECT_EQ(s.n[2], 29);
  EXPECT_EQ(s.n[3], 57);
  EXPECT_TRUE(s.monotone);
  // converged value from a very large reference basis; the last doubling is
  // within sigma_num-ish of it
  MatrixCache big(BasisSpec{.n = 113, .r1 = 0.03, .rmax = 40.0, .eta = 1e-13});
  const double exact =
      solve_channel(model, single(0, 1, 1), th, big).level(2).M;
  EXPECT_LE(std::fabs(s.value.back() - exact), s.sigma_num + 1e-9);
  EXPECT_GE(s.value.back(), exact - 1e-9);
}

TEST(NestedConvergence, WaveFunctionObservablesConvergeMoreSlowlyThanEnergies) {
  // notes 10.8: the energy error is quadratic in the wave-function error,
  // |R(0)|^2 linear.
  NRModel model({.screened = false, .spin = SpinMode::None});
  const Vec th = theta_of(1.5, 1.5, 0.5, 0.18, 0.0, 1.1);
  const BasisSpec start{.n = 10, .r1 = 0.02, .rmax = 40.0, .eta = 1e-13};
  const Observable e1 = [](const SolvedChannel& sc, MatrixCache&) {
    return sc.level(0).M;
  };
  const Observable r0 = [](const SolvedChannel& sc, MatrixCache& c) {
    return R0_squared(sc, 0, c);
  };
  const NestedSequence se =
      nested_convergence(model, th, single(0, 1, 1), start, e1, 3);
  const NestedSequence sr =
      nested_convergence(model, th, single(0, 1, 1), start, r0, 3);
  EXPECT_LT(se.sigma_num / std::fabs(se.value.back()),
            sr.sigma_num / std::fabs(sr.value.back()));
}

TEST(BasisOpt, DatasetBasisLowersTheLevelSumAndStaysAnUpperBound) {
  // the variational objective of a dataset: sum over its channels of the
  // lowest levels the data use (here 1S0, 3S1 up to 2S, 3P2)
  NRModel model({.screened = false, .spin = SpinMode::None});
  const Vec th = theta_of(1.5, 1.5, 0.5, 0.18, 0.0, 1.1);
  Dataset ds;
  for (auto [ch, n] :
       {std::pair{"1S0", 0}, {"3S1", 0}, {"3S1", 1}, {"3P2", 0}}) {
    Datum x;
    x.label = std::string(ch) + "#" + std::to_string(n);
    x.state.channel = ch;
    x.state.n = n;
    ds.d.push_back(x);
  }
  const BasisSpec poor{.n = 12, .r1 = 0.02, .rmax = 30.0, .eta = 1e-13};
  const BasisOptResult r = optimize_dataset_basis(model, th, ds, poor);
  EXPECT_LT(r.objective, r.start_objective - 0.01);  // 14.016 -> 14.000 GeV
  const double exact =
      dataset_trace(model, th, ds,
                    BasisSpec{.n = 80, .r1 = 0.01, .rmax = 60.0, .eta = 1e-13});
  EXPECT_GE(r.objective, exact - 1e-9);
  EXPECT_LT(r.objective - exact, r.start_objective - exact);
  // outside the allowed ranges the objective refuses the point
  EXPECT_TRUE(std::isinf(dataset_trace(
      model, th, ds,
      BasisSpec{.n = 12, .r1 = 5.0, .rmax = 30.0, .eta = 1e-13})));
}

TEST(BasisSize, NestedSequenceConvergesMonotonicallyAndChoosesTheSmallestSize) {
  // n = 10, 19, 37, 73 at fixed ranges: nested bases, so every level can only
  // decrease; the chosen size is the first whose change to the next is below
  // the target
  NRModel model({.screened = false, .spin = SpinMode::None});
  const Vec th = theta_of(1.5, 1.5, 0.5, 0.18, 0.0, 1.1);
  Dataset ds;
  for (auto [ch, n] : {std::pair{"3S1", 0}, {"3S1", 2}, {"3P2", 1}}) {
    Datum x;
    x.label = std::string(ch) + "#" + std::to_string(n);
    x.state.channel = ch;
    x.state.n = n;
    ds.d.push_back(x);
  }
  const BasisSpec spec{.n = 10, .r1 = 0.05, .rmax = 15.0, .eta = 1e-12};
  const BasisSizeChoice c = choose_basis_size(model, th, ds, spec, 10, 3, 1e-4);
  ASSERT_EQ(c.n, (std::vector<int>{10, 19, 37, 73}));
  ASSERT_EQ(c.max_change.size(), 3u);
  EXPECT_GT(c.max_change[0], c.max_change[1]);  // converging
  EXPECT_GT(c.max_change[1], c.max_change[2]);
  EXPECT_LT(c.max_increase, 1e-6);  // variational: no level rises
  ASSERT_GT(c.chosen, 0);
  const auto k = std::find(c.n.begin(), c.n.end(), c.chosen) - c.n.begin();
  EXPECT_LT(c.max_change[k], 1e-4);
  if (k > 0) {
    EXPECT_GE(c.max_change[k - 1], 1e-4);  // the smallest such size
  }
}
