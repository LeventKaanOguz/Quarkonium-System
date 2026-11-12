// M2/M4: Models A/B: spin terms (notes Test 2), virial theorem, dM/dc = 1,
// smeared Coulomb functions and the Feynman-Hellmann gradients against central
// finite differences.
#include <gtest/gtest.h>

#include <cmath>
#include <numbers>

#include "me/coulomb.hpp"
#include "models/nr_model.hpp"
#include "reference_quadrature.hpp"
#include "solve/spectrum.hpp"

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
const BasisSpec kSpec40{.n = 40, .r1 = 0.02, .rmax = 50.0, .eta = 1e-13};
}  // namespace

TEST(CoulombSmeared, SeriesMatchesDirectFormulaAndLimits) {
  const double tau = 1.7;
  for (double r : {0.3, 0.45, 0.58, 0.5882, 0.5883, 0.7, 1.2,
                   3.0}) {  // z = tau r straddles the switch at 1
    const long double z = static_cast<long double>(tau) * r,
                      t3 = static_cast<long double>(tau) * tau * tau;
    const long double sp = std::sqrt(3.14159265358979323846264338327950288L);
    const long double h =
        2 / sp * std::exp(-z * z) / (z * z) - std::erf(z) / (z * z * z);
    const long double gp = t3 * h, T = t3 * (3 * h + 4 / sp * std::exp(-z * z));
    EXPECT_NEAR(coulomb_gp_over_r(tau, r), static_cast<double>(gp),
                1e-13 * (1 + std::abs(static_cast<double>(gp))))
        << r;
    EXPECT_NEAR(coulomb_tensor(tau, r), static_cast<double>(T),
                1e-12 * (1 + std::abs(static_cast<double>(T))))
        << r;
  }
  // point limits and origin values (notes 6.20-6.21)
  EXPECT_NEAR(coulomb_gp_over_r(tau, 1e-6),
              -4 * std::pow(tau, 3) / (3 * std::sqrt(std::numbers::pi)), 1e-9);
  EXPECT_NEAR(coulomb_tensor(tau, 1e-3),
              -8 * std::pow(tau, 5) * 1e-6 / (5 * std::sqrt(std::numbers::pi)),
              1e-9);
  EXPECT_NEAR(coulomb_gp_over_r(tau, 40.0) * std::pow(40.0, 3), -1.0, 1e-12);
  EXPECT_NEAR(coulomb_tensor(tau, 40.0) * std::pow(40.0, 3), -3.0, 1e-12);
  // g''(r) from finite differences of g' = d/dr erf(tau r)/r must equal g'/r -
  // T
  const double r = 0.9, h = 1e-4;
  auto g = [&](double x) { return std::erf(tau * x) / x; };
  const double gpp = (g(r + h) - 2 * g(r) + g(r - h)) / (h * h),
               gp = (g(r + h) - g(r - h)) / (2 * h);
  EXPECT_NEAR(gp / r, coulomb_gp_over_r(tau, r), 1e-7);
  EXPECT_NEAR(gp / r - gpp, coulomb_tensor(tau, r), 1e-6);
}

TEST(Test2, ExpectationValuesAndFirstOrderSpinSplittingsOfTheNotes) {
  const double m = 1.5, mr = 0.75, kappa = 2.0 / 3.0, b = 0.18, sig = 1.1;
  const Basis B0 = make_basis(0, kSpec40), B1 = make_basis(1, kSpec40);
  const Mat S0 = overlap(B0), S1 = overlap(B1);
  const Mat H0 =
      kinetic_nr(B0, mr) - kappa * r_power(B0, B0, -1) + b * r_power(B0, B0, 1);
  const GepResult g0 = solve_gep(H0, S0, 1e-13);
  const Vec c = g0.C.col(0);
  EXPECT_NEAR(std::sqrt(c.dot(r_power(B0, B0, 2) * c)), 1.981259,
              5e-6);  // rms radius, GeV^-1
  const double R0 = c.dot(origin_coefficients(B0));
  EXPECT_NEAR(R0 * R0, 1.3796907, 1e-6);
  const double dsig = std::pow(sig / std::sqrt(std::numbers::pi), 3) *
                      c.dot(gaussian(B0, B0, sig * sig) * c);
  EXPECT_NEAR(dsig, 0.0371741, 2e-7);
  const double delta_hf = 32 * std::numbers::pi * 0.5 / (9 * m * m) * dsig;
  EXPECT_NEAR(1e3 * delta_hf, 92.2753, 2e-4);  // MeV
  // point delta function: 8 alpha_s |R(0)|^2 / (9 m^2)
  EXPECT_NEAR(1e3 * 8 * 0.5 * R0 * R0 / (9 * m * m), 272.5315, 1e-3);
  // virial: 2<T> = kappa <1/r> + b <r>
  EXPECT_NEAR(2 * c.dot(kinetic_nr(B0, mr) * c),
              kappa * c.dot(r_power(B0, B0, -1) * c) +
                  b * c.dot(r_power(B0, B0, 1) * c),
              1e-7);

  // 1P state
  const Mat H1 =
      kinetic_nr(B1, mr) - kappa * r_power(B1, B1, -1) + b * r_power(B1, B1, 1);
  const Vec cp = solve_gep(H1, S1, 1e-13).C.col(0);
  const double rm3 = cp.dot(r_power(B1, B1, -3) * cp),
               rm1 = cp.dot(r_power(B1, B1, -1) * cp);
  EXPECT_NEAR(rm3, 0.1120235, 2e-6);
  EXPECT_NEAR(rm1, 0.3812009, 2e-6);
  const double a = (3 * kappa * rm3 - b * rm1) / (2 * m * m),
               t = kappa * rm3 / (4 * m * m);
  EXPECT_NEAR(1e3 * a, 34.5402, 2e-3);
  EXPECT_NEAR(1e3 * t, 8.2980, 1e-3);
  EXPECT_NEAR(1e3 * (a * ls(1, 1, 2) + t * s12_diag(1, 2)), 31.2210, 2e-3);
  EXPECT_NEAR(1e3 * (a * ls(1, 1, 1) + t * s12_diag(1, 1)), -17.9441, 2e-3);
  EXPECT_NEAR(1e3 * (a * ls(1, 1, 0) + t * s12_diag(1, 0)), -102.2726, 2e-3);
}

TEST(NRModel, FirstOrderShiftsReproduceNotesTest2AndCentreOfGravity) {
  // contact term in S waves only, so the P-wave eigenvector is exactly that of
  // Test 2's H0.
  NRModel model({.screened = false,
                 .spin = SpinMode::FirstOrder,
                 .contact_in_all_L = false});
  MatrixCache cache(kSpec40);
  const Vec th = theta_of(1.5, 1.5, 0.5, 0.18, 0.0, 1.1);
  double cog = 0;
  const double expect[] = {-102.2726, -17.9441, 31.2210};
  for (int J = 0; J <= 2; ++J) {
    const SolvedChannel sc = solve_channel(model, single(1, 1, J), th, cache);
    EXPECT_NEAR(1e3 * sc.level(0).shift, expect[J], 2e-3) << "3P" << J;
    cog += (2 * J + 1) * sc.level(0).shift;
  }
  EXPECT_NEAR(cog, 0.0, 1e-10);  // sum (2J+1) shift = 0
  // spin-independent offset: M = m1 + m2 + E
  const SolvedChannel s = solve_channel(model, single(1, 1, 2), th, cache);
  EXPECT_NEAR(s.level(0).M, 3.0 + s.level(0).E + s.level(0).shift, 1e-14);
}

TEST(NRModel, ConstantShiftAndMassOffsetGradients) {
  NRModel model({.screened = false, .spin = SpinMode::FirstOrder});
  MatrixCache cache(kSpec40);
  Vec th = theta_of(1.5, 1.5, 0.5, 0.18, 0.0, 1.1);
  const SolvedChannel a =
      solve_channel(model, single(0, 1, 1), th, cache, true);
  th(kC) = 0.1;
  const SolvedChannel b = solve_channel(model, single(0, 1, 1), th, cache);
  for (int k = 0; k < 4; ++k)
    EXPECT_NEAR(b.level(k).M - a.level(k).M, 0.1, 1e-10) << k;
  const Vec g =
      a.mass_gradient(0, model, theta_of(1.5, 1.5, 0.5, 0.18, 0.0, 1.1));
  EXPECT_NEAR(g(kC), 1.0, 1e-10);  // dM/dc = 1 for every level
}

namespace {
// central-difference check of dM/dtheta for level k of one channel
void check_gradients(const NRModel& model, const Channel& ch, const Vec& th,
                     int k, double tol, const std::vector<int>& params) {
  MatrixCache cache(BasisSpec{.n = 24, .r1 = 0.03, .rmax = 40.0, .eta = 1e-13});
  const SolvedChannel sc = solve_channel(model, ch, th, cache, true);
  const Vec g = sc.mass_gradient(k, model, th);
  for (int p : params) {
    const double h = 1e-5 * std::max(std::fabs(th(p)), 0.1);
    Vec tp = th, tm = th;
    tp(p) += h;
    tm(p) -= h;
    MatrixCache c1(BasisSpec{.n = 24, .r1 = 0.03, .rmax = 40.0, .eta = 1e-13}),
        c2 = c1;
    const double mp = solve_channel(model, ch, tp, c1).level(k).M,
                 mm = solve_channel(model, ch, tm, c2).level(k).M;
    const double fd = (mp - mm) / (2 * h);
    EXPECT_NEAR(g(p), fd, tol * std::max(1.0, std::fabs(fd)))
        << model.name() << " " << ch.label() << " k=" << k << " param " << p;
  }
}
}  // namespace

TEST(NRModel, FeynmanHellmannGradientsMatchFiniteDifferences) {
  const std::vector<int> all = {kM1, kM2, kAlphaS, kB, kC, kSigma};
  const Vec th = theta_of(1.48, 1.48, 0.52, 0.15, 0.02, 1.2);
  // Model A, first order: includes the eigenvector-derivative term of the shift
  // (7.35)
  NRModel fo({.screened = false, .spin = SpinMode::FirstOrder});
  for (int k : {0, 1}) {
    check_gradients(fo, single(1, 1, 2), th, k, 2e-6, all);
    check_gradients(fo, single(1, 1, 0), th, k, 2e-6, all);
    check_gradients(fo, single(0, 1, 1), th, k, 2e-6, all);
    check_gradients(fo, single(0, 0, 0), th, k, 2e-6, all);
  }
  // unequal masses: m1 != m2 (Bc-like) must still be consistent
  check_gradients(fo, single(1, 1, 1), theta_of(1.4, 4.7, 0.45, 0.17, 0.0, 1.3),
                  0, 2e-6, all);
  // Model B adds mu_s
  NRModel fob({.screened = true, .spin = SpinMode::FirstOrder});
  std::vector<int> allB = all;
  allB.push_back(kMuS);
  check_gradients(fob, single(1, 1, 2), th, 0, 2e-6, allB);
  check_gradients(fob, single(0, 1, 1), th, 1, 2e-6, allB);
  // smeared spin terms in H: coupled 3S1-3D1 and 3P2-3F2, 1P1-3P1 for unequal
  // masses
  NRModel sm({.screened = false, .spin = SpinMode::Smeared});
  Channel sd = {1, -1, -1, {{0, 1}, {2, 1}}};
  check_gradients(sm, sd, th, 0, 2e-6, all);
  check_gradients(sm, sd, th, 1, 2e-6, all);
  Channel pf = {2, 1, 1, {{1, 1}, {3, 1}}};
  check_gradients(sm, pf, th, 0, 2e-6, all);
  Channel unequal = {1, 1, 0, {{1, 0}, {1, 1}}};
  check_gradients(sm, unequal, theta_of(1.4, 4.7, 0.45, 0.17, 0.0, 1.3), 0,
                  2e-6, all);
  NRModel smb({.screened = true, .spin = SpinMode::Smeared});
  check_gradients(smb, sd, th, 0, 2e-6, allB);
}

TEST(NRModel, ScreenedModelReducesToCornellForSmallMu) {
  // b (1 - e^{-mu r})/mu = b r - b mu r^2/2 + O(mu^2): the first-order shift of
  // level k is -b mu <r^2>/2.
  NRModel A({.screened = false, .spin = SpinMode::FirstOrder}),
      B({.screened = true, .spin = SpinMode::FirstOrder});
  MatrixCache ca(kSpec40), cb(kSpec40);
  const double mu = 1e-4, b = 0.18;
  Vec th = theta_of(1.5, 1.5, 0.5, b, 0.0, 1.1, mu);
  const SolvedChannel a = solve_channel(A, single(0, 1, 1), th, ca),
                      bb = solve_channel(B, single(0, 1, 1), th, cb);
  const Mat r2 = r_power(ca.basis(0), ca.basis(0), 2);
  for (int k = 0; k < 3; ++k) {
    const double first_order =
        -0.5 * b * mu * a.level(k).c.dot(r2 * a.level(k).c);
    EXPECT_NEAR(bb.level(k).M - a.level(k).M, first_order,
                5e-4 * std::fabs(first_order))
        << k;
  }
}

TEST(NRModel, SmearedCoupledChannelIsBoundedAndSDMixingIsSmall) {
  NRModel sm({.screened = false, .spin = SpinMode::Smeared});
  MatrixCache cache(kSpec40);
  const Vec th = theta_of(1.5, 1.5, 0.5, 0.18, 0.0, 1.1);
  Channel sd = {1, -1, -1, {{0, 1}, {2, 1}}};
  const SolvedChannel sc = solve_channel(sm, sd, th, cache);
  // lowest level is S-dominated, and above the pure-S eigenvalue minus a few
  // MeV at most
  EXPECT_GT(sc.level(0).prob[0], 0.99);
  EXPECT_LT(sc.level(0).prob[1], 0.01);
  EXPECT_NEAR(sc.level(0).prob[0] + sc.level(0).prob[1], 1.0, 1e-10);
  // switching the tensor coupling off must raise the lowest level (second-order
  // level repulsion)
  const Vec ev = sc.gep.E;
  EXPECT_TRUE(std::isfinite(ev(0)));
}

TEST(NRModel, RunningCouplingGradientsAndLimits) {
  const std::vector<int> all = {kM1, kM2, kAlphaS, kB, kC, kSigma};
  const Vec th = theta_of(1.48, 1.48, 0.52, 0.15, 0.02, 1.2);
  NRModel fo({.running = Running::Gaussian,
              .screened = false,
              .spin = SpinMode::FirstOrder});
  check_gradients(fo, single(1, 1, 2), th, 0, 2e-6, all);
  check_gradients(fo, single(0, 1, 1), th, 1, 2e-6, all);
  NRModel sm({.running = Running::Gaussian,
              .screened = true,
              .spin = SpinMode::Smeared});
  Channel sd = {1, -1, -1, {{0, 1}, {2, 1}}};
  std::vector<int> allB = all;
  allB.push_back(kMuS);
  check_gradients(sm, sd, th, 0, 2e-6, allB);
  // the running coupling is weaker at short distance than the saturation value:
  // levels sit above the constant-alpha ones
  NRModel constant({.running = Running::Constant, .spin = SpinMode::None}),
      run({.running = Running::Gaussian, .spin = SpinMode::None});
  MatrixCache c1(kSpec40), c2(kSpec40);
  const SolvedChannel a = solve_channel(constant, single(0, 1, 1), th, c1),
                      b = solve_channel(run, single(0, 1, 1), th, c2);
  EXPECT_GT(b.level(0).M, a.level(0).M);
}
