// Running coupling (notes 5): forms, limits, the static potential and the
// smeared Coulomb functions.
#include <gtest/gtest.h>

#include <cmath>
#include <numbers>

#include "me/coulomb.hpp"
#include "models/alpha_s.hpp"

using namespace qgem;
using std::numbers::pi;

TEST(AlphaS, FormsAndTheirLimits) {
  EXPECT_NEAR(beta0(3), 9.0, 1e-15);
  EXPECT_NEAR(beta0(4), 25.0 / 3.0, 1e-15);
  // Richardson: alpha -> 4 pi Lambda^2/(beta0 Q^2) as Q -> 0 ("confining"),
  // one-loop at large Q
  const auto rich = make_richardson_coupling(0.4, 3);
  EXPECT_NEAR(rich->alpha(1e-7) * 1e-7 / (4 * pi * 0.16 / 9.0), 1.0,
              1e-6);  // ln(1+x) = x - x^2/2: corrections ~ Q^2/2 Lambda^2
  EXPECT_NEAR(rich->alpha(1e6), 4 * pi / (9.0 * std::log(1e6 / 0.16)), 1e-4);
  // APT: finite at Q = 0 with the value 4 pi / beta0; finite at Q^2 = Lambda^2
  // (bracket -> 1/2)
  const auto apt = make_apt_coupling(0.35, 3);
  EXPECT_NEAR(apt->alpha(0.0), 4 * pi / 9.0, 1e-12);
  EXPECT_NEAR(apt->alpha(0.35 * 0.35), 4 * pi / 9.0 * 0.5, 1e-9);
  EXPECT_NEAR(apt->alpha(0.35 * 0.35 * (1 + 1e-4)) / apt->alpha(0.35 * 0.35),
              1.0, 1e-3);  // no pole
  // frozen (gluon mass)
  const auto fr = make_frozen_coupling(0.3, 0.5, 3);
  EXPECT_NEAR(fr->alpha(0.0), 4 * pi / (9.0 * std::log(4 * 0.25 / 0.09)),
              1e-12);
  // Godfrey-Isgur Gaussian sum saturates at 0.6 and falls off
  const auto gi = make_gaussian_sum_coupling(gi_alpha(), gi_gamma());
  EXPECT_NEAR(gi->alpha(0.0), 0.6, 1e-15);
  EXPECT_LT(gi->alpha(1e5), 1e-3);  // only the widest Gaussian (4 gamma_3^2 =
                                    // 1000) survives to Q^2 ~ 100
  EXPECT_THROW(make_frozen_coupling(0.6, 0.2, 3), std::invalid_argument);
}

TEST(AlphaS, StaticPotentialOfAConstantAndOfTheRichardsonForm) {
  // constant alpha: V = -4 alpha / (3 r)
  const auto c = make_constant_coupling(0.4);
  for (double r : {0.2, 1.0, 5.0})
    EXPECT_NEAR(static_potential(*c, r) * r / (-4 * 0.4 / 3), 1.0, 1e-6) << r;
  // Richardson: V(r) -> b r + const with b = 8 pi Lambda^2/(33 - 2 nf)
  // (Exercise 5.3)
  const double L = 0.4;
  const int nf = 3;
  const auto rich = make_richardson_coupling(L, nf);
  const double b = 8 * pi * L * L / (33 - 2 * nf);
  const double r = 25.0, h = 1.0;
  const double slope =
      (static_potential(*rich, r + h) - static_potential(*rich, r - h)) /
      (2 * h);
  EXPECT_NEAR(slope / b, 1.0,
              0.02);  // ~1% finite-r corrections at r = 25 GeV^-1
  // the Richardson V(r) carries an infrared-divergent additive constant ("up to
  // an infinite constant", Ex. 5.3): only slopes are meaningful short distance:
  // V r -> -(4/3) alpha(1/r^2)-like, i.e. asymptotically free: |V r| decreases
  // as r -> 0
  EXPECT_LT(std::fabs(static_potential(*rich, 0.02) * 0.02),
            std::fabs(static_potential(*rich, 0.5) * 0.5));
  // APT is finite at all r and approaches -(4/3) (4 pi/beta0) ... never exceeds
  // its IR value 4 pi/beta0 in |V r|
  const auto apt = make_apt_coupling(0.35, 3);
  for (double rr : {0.05, 0.5, 3.0})
    EXPECT_LT(std::fabs(static_potential(*apt, rr) * rr),
              4.0 / 3.0 * 4 * pi / 9.0 + 1e-9)
        << rr;
  const auto g = make_gaussian_sum_coupling(gi_alpha(), gi_gamma());
  EXPECT_THROW(static_potential(*g, 1.0), std::invalid_argument);
}

TEST(AlphaS, SmearedCoulombTransformReproducesTheAnalyticGaussianSum) {
  // Smearing the GI Gaussian sum with width sigma composes the widths (6.7):
  // analytic erf forms. The numerical transform of eq. (6.x) must reproduce G,
  // G'/r, lap G, T for every r.
  const auto g = make_gaussian_sum_coupling(gi_alpha(), gi_gamma());
  const double sigma = 1.9;
  const auto al = gi_alpha();
  const auto ga = gi_gamma();
  for (double r : {0.01, 0.05, 0.3, 1.0, 2.5, 6.0}) {
    const SmearedCoulombValues v = smeared_coulomb_transform(*g, sigma, r);
    double G = 0, Gp = 0, lap = 0, T = 0;
    for (std::size_t k = 0; k < al.size(); ++k) {
      const double tau =
          sigma * ga[k] / std::sqrt(ga[k] * ga[k] + sigma * sigma);
      G += -4.0 / 3 * al[k] * std::erf(tau * r) / r;
      Gp += -4.0 / 3 * al[k] * coulomb_gp_over_r(tau, r);
      T += -4.0 / 3 * al[k] * coulomb_tensor(tau, r);
      lap += 16 * pi / 3 * al[k] * std::pow(tau, 3) / std::pow(pi, 1.5) *
             std::exp(-tau * tau * r * r);
    }
    EXPECT_NEAR(v.G, G, 1e-9 * (1 + std::fabs(G))) << r;
    EXPECT_NEAR(v.Gp_over_r, Gp, 1e-7 * (1 + std::fabs(Gp))) << r;
    EXPECT_NEAR(v.lap, lap, 1e-8 * (1 + std::fabs(lap))) << r;
    EXPECT_NEAR(v.T, T, 1e-6 * (1 + std::fabs(T))) << r;
  }
}

TEST(AlphaS, SplineTableReproducesTheTransform) {
  const auto rich = make_frozen_coupling(
      0.3, 0.5, 3);  // finite at Q = 0: no infrared-divergent constant
  const double sigma = 1.5;
  const SmearedCoulombTable tab(*rich, sigma);
  for (double r : {0.003, 0.02, 0.2, 0.9, 3.0, 12.0, 40.0}) {
    const SmearedCoulombValues v = smeared_coulomb_transform(*rich, sigma, r);
    EXPECT_NEAR(tab.G(r), v.G, 2e-6 * (1 + std::fabs(v.G))) << r;
    EXPECT_NEAR(tab.Gp_over_r(r), v.Gp_over_r,
                2e-6 * (1 + std::fabs(v.Gp_over_r)))
        << r;
    EXPECT_NEAR(tab.lap(r), v.lap, 2e-6 * (1 + std::fabs(v.lap))) << r;
    EXPECT_NEAR(tab.T(r), v.T, 2e-6 * (1 + std::fabs(v.T))) << r;
  }
  // smearing cures the Coulomb singularity: G is finite at the origin
  EXPECT_TRUE(std::isfinite(tab.G(1e-9)));
}

TEST(AlphaS, GaussianSumShapesOfFiniteCouplings) {
  // frozen and APT couplings are finite at Q = 0: their shapes alpha/alpha(0)
  // are fitted by non-negative Gaussian sums with sum w = 1, accurate to a few
  // per cent over 0 < Q < 50 GeV (the shapes fall like 1/ln Q^2, Gaussians do
  // not, so the tail sets the error).
  for (const auto& a :
       {make_frozen_coupling(0.30, 0.50, 4), make_apt_coupling(0.30, 4)}) {
    const GaussianShape g = fit_gaussian_shape(*a, default_shape_gamma());
    double tot = 0;
    for (double w : g.w) {
      EXPECT_GE(w, 0.0);
      tot += w;
    }
    EXPECT_NEAR(tot, 1.0, 1e-12) << a->name();
    EXPECT_LT(g.max_rel_err, 0.05) << a->name();
  }
  // a confining (divergent at Q = 0) coupling has no such shape
  EXPECT_THROW(fit_gaussian_shape(*make_richardson_coupling(0.4, 3),
                                  default_shape_gamma()),
               std::invalid_argument);
}
