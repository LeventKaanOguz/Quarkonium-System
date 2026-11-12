// M1c: independent cross-checks of the GEM solver (notes roadmap V5, V6): the
// linear potential against Airy zeros, and Cornell/screened against a separate
// finite-difference solver.
#include <gtest/gtest.h>

#include <cmath>

#include "basis/basis.hpp"
#include "linalg/gep.hpp"
#include "me/radial.hpp"
#include "reference_fd.hpp"

using namespace qgem;

TEST(Independent, LinearPotentialAgainstAiryZeros) {
  // -u''/(2 mr) + b r u = E u, u(0)=0:  E_n = (b^2 / 2 mr)^{1/3} |a_n|, a_n the
  // zeros of Ai.
  const double mr = 0.75, b = 0.18;
  const double airy[] = {2.338107410459767, 4.087949444130970,
                         5.520559828095551, 6.786708090071759};
  const BasisSpec spec{.n = 40, .r1 = 0.02, .rmax = 60.0};
  const Basis B = make_basis(0, spec);
  const Vec E =
      solve_gep(kinetic_nr(B, mr) + b * r_power(B, B, 1), overlap(B)).E;
  const double scale = std::cbrt(b * b / (2 * mr));
  const double tol[] = {
      5e-8, 5e-7, 3e-6,
      1e-5};  // excited levels need more of the basis (cf. notes Test 1)
  for (int n = 0; n < 4; ++n) {
    EXPECT_GE(E(n), scale * airy[n] - 1e-11) << n;  // upper bound
    EXPECT_NEAR(E(n), scale * airy[n], tol[n]) << n;
  }
  // l = 0 states of the linear potential are l = 0 only: the l = 1 ladder obeys
  // no Airy formula, but the finite-difference check below covers it.
}

TEST(Independent, CornellAndScreenedAgainstFiniteDifference) {
  const double mr = 0.75, kappa = 2.0 / 3.0, b = 0.18, mu = 0.10;
  auto vA = [&](double r) { return -kappa / r + b * r; };
  auto vB = [&](double r) {
    return -kappa / r + b / mu * (1 - std::exp(-mu * r));
  };
  const BasisSpec spec{.n = 40, .r1 = 0.02, .rmax = 50.0, .eta = 1e-13};
  for (int l : {0, 1, 2}) {
    const Basis B = make_basis(l, spec);
    const Mat S = overlap(B), T = kinetic_nr(B, mr),
              C = -kappa * r_power(B, B, -1);
    const Vec EA = solve_gep(T + C + b * r_power(B, B, 1), S, 1e-13).E;
    const Vec EB =
        solve_gep(T + C + b / mu * (S - exp_power(B, B, mu, 0)), S, 1e-13).E;
    for (int k = 0; k < 2; ++k) {
      // FD box 40 GeV^-1 is far beyond the 1S/2S extent; N = 8000 -> 32000 grid
      // points.
      const double fa = test::fd_level(vA, l, mr, 40.0, 8000, k),
                   fb = test::fd_level(vB, l, mr, 40.0, 8000, k);
      EXPECT_NEAR(EA(k), fa, 1e-6) << "A l=" << l << " k=" << k;
      EXPECT_NEAR(EB(k), fb, 1e-6) << "B l=" << l << " k=" << k;
    }
  }
}
