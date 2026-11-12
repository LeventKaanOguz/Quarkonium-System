// M1b/M1c: eigensolver and spin-independent Models A/B against analytic results
// and the independently computed reference numbers of the notes (12.2, Test 1).
#include <gtest/gtest.h>

#include <cmath>

#include "basis/basis.hpp"
#include "linalg/gep.hpp"
#include "me/radial.hpp"

using namespace qgem;

namespace {
Vec cornell(const BasisSpec& spec, int l, double mr, double kappa, double b,
            double mu = 0.0) {
  const Basis B = make_basis(l, spec);
  Mat V = (mu > 0) ? Mat(b / mu * (overlap(B) - exp_power(B, B, mu, 0)))
                   : Mat(b * r_power(B, B, 1));
  Mat H = kinetic_nr(B, mr) - kappa * r_power(B, B, -1) + V;
  return solve_gep(H, overlap(B), spec.eta).E;
}
}  // namespace

TEST(Gep, CanonicalAndCholeskyAgreeOnWellConditionedProblem) {
  const Basis B = make_basis(0, {.n = 8, .r1 = 0.1, .rmax = 10.0});
  Mat H = kinetic_nr(B, 0.75) - (2.0 / 3.0) * r_power(B, B, -1) +
          0.18 * r_power(B, B, 1);
  const Mat S = overlap(B);
  const GepResult a = solve_gep(H, S), b = solve_gep_cholesky(H, S);
  EXPECT_LT((a.E - b.E).cwiseAbs().maxCoeff(), 1e-10);
  // S-orthonormal eigenvectors
  const Mat G = a.C.transpose() * S * a.C;
  EXPECT_LT((G - Mat::Identity(a.n_kept, a.n_kept)).cwiseAbs().maxCoeff(),
            1e-10);
  // X^T S X = 1
  EXPECT_LT((a.X.transpose() * S * a.X - Mat::Identity(a.n_kept, a.n_kept))
                .cwiseAbs()
                .maxCoeff(),
            1e-10);
}

TEST(Gep, DropsNearNullDirectionsOfAnIllConditionedOverlap) {
  // q = 1.05: neighbouring functions are almost parallel; the Cholesky route is
  // not trusted here.
  const Basis B = make_basis(0, {.n = 60, .r1 = 0.1, .rmax = 1.1});
  const Mat S = overlap(B);
  const GepResult r = solve_gep(Mat::Zero(60, 60), S, 1e-12);
  EXPECT_LT(r.n_kept, 60);
  EXPECT_GT(r.n_kept, 5);
}

TEST(Solver, HarmonicOscillatorLadder) {
  // E = omega (2 n_r + l + 3/2) with V = mr omega^2 r^2 / 2. A geometric
  // Gaussian basis represents the polynomial-times-Gaussian oscillator states
  // only through cancellations, so the error grows with the radial excitation;
  // it must always be >= 0 (HUM bound) and fall when the basis is enlarged.
  const double mr = 0.75, omega = 1.3;
  for (int l : {0, 1, 2}) {
    const BasisSpec spec{.n = 40, .r1 = 0.1, .rmax = 10.0};
    const Basis B = make_basis(l, spec);
    const Mat H =
        kinetic_nr(B, mr) + 0.5 * mr * omega * omega * r_power(B, B, 2);
    const Vec E = solve_gep(H, overlap(B)).E;
    const double tol[] = {1e-9, 5e-8, 1e-6, 5e-5};
    for (int nr = 0; nr < 4; ++nr) {
      const double ex = omega * (2 * nr + l + 1.5);
      EXPECT_GE(E(nr), ex - 1e-11) << "l=" << l << " nr=" << nr;
      EXPECT_NEAR(E(nr), ex, tol[nr]) << "l=" << l << " nr=" << nr;
    }
  }
}

TEST(Solver, OscillatorErrorFallsUnderNestedDoubling) {
  const double mr = 0.75, omega = 1.3;
  BasisSpec spec{.n = 11, .r1 = 0.1, .rmax = 10.0};
  double prev_err[3] = {1e9, 1e9, 1e9};
  for (int step = 0; step < 3; ++step, spec.n = 2 * spec.n - 1) {
    const Basis B = make_basis(0, spec);
    const Vec E = solve_gep(kinetic_nr(B, mr) +
                                0.5 * mr * omega * omega * r_power(B, B, 2),
                            overlap(B))
                      .E;
    for (int nr = 0; nr < 3; ++nr) {
      const double err = E(nr) - omega * (2 * nr + 1.5);
      EXPECT_GE(err, -1e-11);
      EXPECT_LE(err, prev_err[nr] + 1e-11) << "n=" << spec.n << " nr=" << nr;
      prev_err[nr] = err;
    }
  }
}

TEST(Solver, SingleGaussianOscillatorIsExactToRoundoff) {
  // One Gaussian with nu = mr omega / 2 is the exact oscillator ground state.
  const double mr = 0.75, omega = 1.3, nu = mr * omega / 2;
  const double r0 = 1 / std::sqrt(nu);
  const Basis B = make_basis(0, {.n = 1, .r1 = r0, .rmax = r0});
  const Mat H = kinetic_nr(B, mr) + 0.5 * mr * omega * omega * r_power(B, B, 2);
  EXPECT_NEAR(solve_gep(H, overlap(B)).E(0), 1.5 * omega, 1e-14);
}

TEST(Solver, HydrogenLevels) {
  // E_n = -mr kappa^2 / (2 n^2), n = n_r + l + 1.  (The notes' list "-0.22,
  // -0.05, -0.02469" rounds loosely: the exact values for mr = 1, kappa = 2/3
  // are -2/9, -1/18, -2/81.)
  const double mr = 1.0, kappa = 2.0 / 3.0;
  const BasisSpec spec{.n = 60, .r1 = 0.005, .rmax = 200.0};
  for (int l : {0, 1, 2}) {
    const Basis B = make_basis(l, spec);
    const Vec E =
        solve_gep(kinetic_nr(B, mr) - kappa * r_power(B, B, -1), overlap(B)).E;
    for (int nr = 0; nr < 3; ++nr) {
      const double n = nr + l + 1, ex = -mr * kappa * kappa / (2 * n * n);
      EXPECT_NEAR(E(nr), ex, 2e-6) << "l=" << l << " nr=" << nr;
      EXPECT_GE(E(nr), ex - 1e-10);  // round-off level
    }
  }
}

TEST(Solver,
     ConstantShiftMovesEveryLevelByExactlyTheShift) {  // dM/dc = 1 because
                                                       // dH/dc = S
  const Basis B = make_basis(0, {.n = 20, .r1 = 0.05, .rmax = 30.0});
  const Mat S = overlap(B);
  const Mat H = kinetic_nr(B, 0.75) - (2.0 / 3.0) * r_power(B, B, -1) +
                0.18 * r_power(B, B, 1);
  const Vec E0 = solve_gep(H, S).E, E1 = solve_gep(H + 0.1 * S, S).E;
  EXPECT_NEAR(((E1 - E0).array() - 0.1).abs().maxCoeff(), 0.0, 1e-10);
}

TEST(Solver, NestedDoublingDecreasesEveryLevelMonotonically) {
  // n -> 2n-1 keeps all old functions: Ritz values can only go down (Cauchy
  // interlacing).
  BasisSpec spec{.n = 6, .r1 = 0.02, .rmax = 50.0};
  Vec prev = cornell(spec, 0, 0.75, 2.0 / 3.0, 0.18);
  for (int step = 0; step < 3; ++step) {
    spec.n = 2 * spec.n - 1;
    const Vec cur = cornell(spec, 0, 0.75, 2.0 / 3.0, 0.18);
    for (int k = 0; k < 4; ++k)
      EXPECT_LE(cur(k), prev(k) + 1e-10) << "n=" << spec.n << " k=" << k;
    prev = cur;
  }
}

// ---- notes 12.2, Test 1: m1 = m2 = 1.5 (mr = 0.75), alpha_s = 0.5 (kappa =
// 2/3), b = 0.18, c = 0. Entries are E = M - 2m in GeV. The FD column is the
// independent finite-difference reference.
struct T1Row {
  const char* name;
  int l, nr;
  double a20, a40, aFD, b20, b40, bFD;
};
static const T1Row kTest1[] = {
    {"1S", 0, 0, 0.2164465, 0.2164462, 0.2164461, 0.1823122, 0.1823119,
     0.1823119},
    {"2S", 0, 1, 0.8673229, 0.8672945, 0.8672945, 0.7227575, 0.7227531,
     0.7227531},
    {"3S", 0, 2, 1.3260325, 1.3240175, 1.3240174, 1.0349623, 1.0349507,
     1.0349507},
    {"1P", 1, 0, 0.6993398, 0.6993328, 0.6993328, 0.6062354, 0.6062335,
     0.6062335},
    {"2P", 1, 1, 1.1753880, 1.1753865, 1.1753865, 0.9492210, 0.9491315,
     0.9491315},
    {"3P", 1, 2, 1.5780260, 1.5719843, 1.5719842, 1.1895326, 1.1887053,
     1.1887053},
    {"1D", 2, 0, 1.0093770, 1.0093090, 1.0093090, 0.8498083, 0.8498066,
     0.8498066},
    {"2D", 2, 1, 1.4216081, 1.4212850, 1.4212847, 1.1135569, 1.1134266,
     1.1134266},
    {"3D", 2, 2, 1.7927322, 1.7838549, 1.7838539, 1.3109787, 1.3085364,
     1.3085363},
};

TEST(Test1, ModelAAndBReproduceTheNotesReferenceTable) {
  const double mr = 0.75, kappa = 2.0 / 3.0, b = 0.18, muS = 0.10;
  const BasisSpec s20{.n = 20, .r1 = 0.05, .rmax = 30.0, .eta = 1e-13};
  const BasisSpec s40{.n = 40, .r1 = 0.02, .rmax = 50.0, .eta = 1e-13};
  for (const T1Row& r : kTest1) {
    EXPECT_NEAR(cornell(s20, r.l, mr, kappa, b)(r.nr), r.a20, 6e-8)
        << r.name << " A n=20";
    EXPECT_NEAR(cornell(s40, r.l, mr, kappa, b)(r.nr), r.a40, 6e-8)
        << r.name << " A n=40";
    EXPECT_NEAR(cornell(s20, r.l, mr, kappa, b, muS)(r.nr), r.b20, 6e-8)
        << r.name << " B n=20";
    EXPECT_NEAR(cornell(s40, r.l, mr, kappa, b, muS)(r.nr), r.b40, 6e-8)
        << r.name << " B n=40";
    // the n=40 GEM level is within its documented numerical error of the
    // finite-difference value
    EXPECT_NEAR(cornell(s40, r.l, mr, kappa, b)(r.nr), r.aFD, 1.5e-6)
        << r.name << " A vs FD";
    EXPECT_NEAR(cornell(s40, r.l, mr, kappa, b, muS)(r.nr), r.bFD, 1.5e-6)
        << r.name << " B vs FD";
  }
}
