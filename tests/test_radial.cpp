// M1a: closed-form radial matrix elements against independent quadrature.
#include <gtest/gtest.h>

#include <cmath>
#include <numbers>

#include "basis/basis.hpp"
#include "me/radial.hpp"
#include "reference_quadrature.hpp"

using namespace qgem;

namespace {
// R_nl(r) with the normalization recomputed here by quadrature of r^2 R^2 (not
// from the closed form).
struct RefFn {
  int l;
  double nu, N;
  double operator()(double r) const {
    return N * std::pow(r, l) * std::exp(-nu * r * r);
  }
};
RefFn ref_fn(int l, double nu) {
  const double rc = std::sqrt(45.0 / (2 * nu));
  const double n2 = qgem::test::integrate(
      [&](double r) {
        return std::pow(r, 2 * l + 2) * std::exp(-2 * nu * r * r);
      },
      0, rc, 200);
  return {l, nu, 1.0 / std::sqrt(n2)};
}
double ref_me(const RefFn& f, const RefFn& g,
              const std::function<double(double)>& O) {
  const double a = f.nu + g.nu, rc = std::sqrt(45.0 / a);
  return qgem::test::integrate(
      [&](double r) { return r * r * f(r) * O(r) * g(r); }, 0, rc, 400);
}
}  // namespace

TEST(Radial, BasisNormalizationAndGeometricRanges) {
  const Basis b = make_basis(1, {.n = 10, .r1 = 0.05, .rmax = 30.0});
  EXPECT_NEAR(1 / std::sqrt(b.nu(0)), 0.05, 1e-15);
  EXPECT_NEAR(1 / std::sqrt(b.nu(9)), 30.0, 1e-12);
  const double q = 1 / std::sqrt(b.nu(1)) / (1 / std::sqrt(b.nu(0)));
  EXPECT_NEAR(1 / std::sqrt(b.nu(5)), 0.05 * std::pow(q, 5), 1e-12);
  const Mat S = overlap(b);
  for (int i = 0; i < 10; ++i) EXPECT_NEAR(S(i, i), 1.0, 1e-15);
}

TEST(Radial, OverlapKineticAndPowersMatchQuadrature) {
  for (int l : {0, 1, 2}) {
    const Basis b = make_basis(l, {.n = 6, .r1 = 0.1, .rmax = 8.0});
    const Mat S = overlap(b), P2 = p2(b), r1 = r_power(b, b, 1),
              rm1 = r_power(b, b, -1), r2 = r_power(b, b, 2);
    Mat rm3;
    if (l >= 1) {
      rm3 = r_power(b, b, -3);
    }
    for (int i = 0; i < 6; ++i)
      for (int j = 0; j < 6; ++j) {
        const RefFn f = ref_fn(l, b.nu(i)), g = ref_fn(l, b.nu(j));
        EXPECT_NEAR(S(i, j), ref_me(f, g, [](double) { return 1.0; }), 1e-12);
        EXPECT_NEAR(r1(i, j) / ref_me(f, g, [](double r) { return r; }), 1.0,
                    1e-12);
        EXPECT_NEAR(rm1(i, j) / ref_me(f, g, [](double r) { return 1 / r; }),
                    1.0, 1e-12);
        EXPECT_NEAR(r2(i, j) / ref_me(f, g, [](double r) { return r * r; }),
                    1.0, 1e-12);
        if (l >= 1) {
          EXPECT_NEAR(
              rm3(i, j) /
                  ref_me(f, g, [](double r) { return 1 / (r * r * r); }),
              1.0, 1e-11);
        }
        // <p^2> = Int r^2 [R_f' R_g' + l(l+1)/r^2 R_f R_g] dr (integration by
        // parts; boundary terms vanish) with the analytic derivative R' = (l/r
        // - 2 nu r) R: independent of Lemma 7.6.
        const double a = f.nu + g.nu, rc = std::sqrt(45.0 / a);
        const double refp2 = qgem::test::integrate(
            [&](double r) {
              const double df = (l / r - 2 * f.nu * r) * f(r),
                           dg = (l / r - 2 * g.nu * r) * g(r);
              return r * r * df * dg + l * (l + 1) * f(r) * g(r);
            },
            0, rc, 400);
        EXPECT_NEAR(P2(i, j) / refp2, 1.0, 1e-11)
            << "l=" << l << " " << i << "," << j;
      }
    // symmetry
    EXPECT_LT((P2 - P2.transpose()).cwiseAbs().maxCoeff(),
              1e-12 * P2.cwiseAbs().maxCoeff());
  }
}

TEST(Radial, GaussianExpPowerAndErfMatchQuadrature) {
  const Basis b = make_basis(0, {.n = 5, .r1 = 0.2, .rmax = 6.0});
  const Basis bp = make_basis(1, {.n = 5, .r1 = 0.2, .rmax = 6.0});
  const double beta = 1.21, mu = 0.35, tau = 2.0;
  const Mat G = gaussian(b, b, beta), E = erf_over_r(b, b, tau),
            Ex = exp_power(b, b, mu, 0), Ex1 = exp_power(b, b, mu, -1);
  const Mat Gx = gaussian(b, bp, beta), E01 = erf_over_r(b, bp, tau),
            Rx = r_power(b, bp, -1);
  for (int i = 0; i < 5; ++i)
    for (int j = 0; j < 5; ++j) {
      const RefFn f = ref_fn(0, b.nu(i)), g = ref_fn(0, b.nu(j)),
                  g1 = ref_fn(1, bp.nu(j));
      EXPECT_NEAR(
          G(i, j) /
              ref_me(f, g, [&](double r) { return std::exp(-beta * r * r); }),
          1.0, 1e-12);
      EXPECT_NEAR(
          E(i, j) /
              ref_me(f, g, [&](double r) { return std::erf(tau * r) / r; }),
          1.0, 1e-12);
      EXPECT_NEAR(
          Ex(i, j) / ref_me(f, g, [&](double r) { return std::exp(-mu * r); }),
          1.0, 1e-12);
      EXPECT_NEAR(
          Ex1(i, j) /
              ref_me(f, g, [&](double r) { return std::exp(-mu * r) / r; }),
          1.0, 1e-12);
      // different L in bra and ket (tensor coupling, transitions): l_bar = 1/2
      EXPECT_NEAR(
          Gx(i, j) /
              ref_me(f, g1, [&](double r) { return std::exp(-beta * r * r); }),
          1.0, 1e-12);
      EXPECT_NEAR(
          E01(i, j) /
              ref_me(f, g1, [&](double r) { return std::erf(tau * r) / r; }),
          1.0, 1e-12);
      EXPECT_NEAR(Rx(i, j) / ref_me(f, g1, [](double r) { return 1 / r; }), 1.0,
                  1e-12);
    }
  // closed form for l_bar = 0: <erf(tau r)/r> = N'N tau / (2 a sqrt(a + tau^2))
  for (int i = 0; i < 5; ++i) {
    const double a = 2 * b.nu(i);
    EXPECT_NEAR(
        E(i, i),
        b.norm(i) * b.norm(i) * tau / (2 * a * std::sqrt(a + tau * tau)),
        1e-14);
  }
  // tau -> infinity recovers <1/r>
  const Mat Einf = erf_over_r(b, b, 1e6), R1 = r_power(b, b, -1);
  EXPECT_LT((Einf - R1).cwiseAbs().maxCoeff(), 1e-9);
}

TEST(Radial, LocalQuadratureReproducesClosedForms) {
  const Basis b = make_basis(1, {.n = 8, .r1 = 0.05, .rmax = 20.0});
  const Mat Lr = local(b, b, [](double r) { return r; }), Cr = r_power(b, b, 1);
  EXPECT_LT((Lr - Cr).cwiseAbs().maxCoeff(), 1e-12);
  const Mat Lm = local(b, b, [](double r) { return std::exp(-0.4 * r) / r; }),
            Cm = exp_power(b, b, 0.4, -1);
  EXPECT_LT((Lm - Cm).cwiseAbs().maxCoeff(), 1e-12);
}
