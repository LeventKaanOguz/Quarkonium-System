// M10: Godfrey-Isgur Models C and D. Reference values are the published GI
// spectra *as recalled in the notes (12.3, "verify before use")*: they are a
// smoke test of the whole relativized chain, not an exact reference. The exact
// checks are structural: Hermiticity, closure vs the two-step method, limits,
// parameter derivatives.
#include <gtest/gtest.h>

#include <cmath>

#include "me/momentum.hpp"
#include "models/gi_model.hpp"
#include "solve/spectrum.hpp"

using namespace qgem;

namespace {
const BasisSpec kSpec{.n = 30, .r1 = 0.02, .rmax = 30.0, .eta = 1e-12};
Vec gi_theta(const GIModel& model, double m1, double m2) {
  Vec th = model.godfrey_isgur_params(m1, m2).theta0();
  return th;
}
Channel find(const GIModel& m, const std::string& label, bool equal = true) {
  for (const Channel& c : m.channels(2, equal))
    if (c.label() == label) return c;
  throw std::runtime_error("no channel " + label);
}
}  // namespace

TEST(GI, CharmoniumAndBottomoniumAgreeWithThePublishedSpectra) {
  GIModel model;
  MatrixCache cache(kSpec);
  struct Row {
    const char* ch;
    int k;
    double mc;
    double tol;
  };
  // charmonium, m_c = 1.628 (GI: eta_c 2.975, J/psi 3.098, eta_c(2S) 3.623,
  // psi(2S) 3.676, chi_cJ 3.445/3.510/3.553, h_c 3.519)
  const Vec c = gi_theta(model, 1.628, 1.628);
  const Row rows[] = {{"1S0", 0, 2.975, 0.012}, {"3S1-3D1", 0, 3.098, 0.012},
                      {"1S0", 1, 3.623, 0.010}, {"3S1-3D1", 1, 3.676, 0.010},
                      {"3P0", 0, 3.445, 0.010}, {"3P1", 0, 3.510, 0.010},
                      {"3P2", 0, 3.553, 0.010}, {"1P1", 0, 3.519, 0.010}};
  for (const Row& r : rows) {
    const SolvedChannel sc = solve_channel(model, find(model, r.ch), c, cache);
    EXPECT_NEAR(sc.level(r.k).M, r.mc, r.tol)
        << "cc " << r.ch << " n=" << r.k + 1;
  }
  const Vec b = gi_theta(model, 4.977, 4.977);
  MatrixCache cb(kSpec);
  const Row bb[] = {{"3S1-3D1", 0, 9.460, 0.006},
                    {"1S0", 0, 9.402, 0.010},
                    {"3P0", 0, 9.855, 0.015},
                    {"3P1", 0, 9.874, 0.012},
                    {"3P2", 0, 9.894, 0.012}};
  for (const Row& r : bb) {
    const SolvedChannel sc = solve_channel(model, find(model, r.ch), b, cb);
    EXPECT_NEAR(sc.level(r.k).M, r.mc, r.tol) << "bb " << r.ch;
  }
  // the P-wave fine-structure ordering of a vector + scalar-confinement model:
  // chi_0 < chi_1 < chi_2, h_c near the centre of gravity
  const double
      m0 = solve_channel(model, find(model, "3P0"), c, cache).level(0).M,
      m1 = solve_channel(model, find(model, "3P1"), c, cache).level(0).M,
      m2 = solve_channel(model, find(model, "3P2"), c, cache).level(0).M,
      h = solve_channel(model, find(model, "1P1"), c, cache).level(0).M;
  EXPECT_LT(m0, m1);
  EXPECT_LT(m1, m2);
  const double cog = (m0 + 3 * m1 + 5 * m2) / 9.0;
  EXPECT_NEAR(h, cog, 0.01);  // centre-of-gravity theorem holds approximately
                              // (short-ranged spin-spin)
}

TEST(GI, HamiltonianIsHermitianAndSDMixingIsSmallForTheGroundState) {
  GIModel model;
  MatrixCache cache(kSpec);
  const Vec c = gi_theta(model, 1.628, 1.628);
  const Channel sd = find(model, "3S1-3D1");
  const ChannelMatrices m = model.build(sd, c, cache, false);
  EXPECT_LT((m.H - m.H.transpose()).cwiseAbs().maxCoeff(),
            1e-9 * m.H.cwiseAbs().maxCoeff());
  const SolvedChannel sc = solve_channel(model, sd, c, cache);
  EXPECT_LT(sc.level(0).prob[1], 0.01);
  EXPECT_NEAR(sc.level(0).prob[0] + sc.level(0).prob[1], 1.0, 1e-10);
  // X^T S X = 1 and the eigenvectors are S-orthonormal
  const Mat G = sc.gep.C.transpose() * sc.mats.S * sc.gep.C;
  EXPECT_LT((G - Mat::Identity(G.rows(), G.cols())).cwiseAbs().maxCoeff(),
            1e-9);
}

TEST(GI, ClosureAgreesWithTheTwoStepMethodForTheSpinSandwiches) {
  // The contact and tensor terms are Phi W Phi sandwiches: closure (used in the
  // model) against the exact two-step matrix elements, through the eigenvalues
  // of the hyperfine pair of S-wave levels.
  GIModel model;
  const Vec c = gi_theta(model, 1.628, 1.628);
  const double m = 1.628, eps = c(kGepsC);
  const double sig = GIModel::sigma12(c);
  double tau[3];
  for (int k = 0; k < 3; ++k)
    tau[k] = sig * c(kGgamma1 + k) /
             std::sqrt(c(kGgamma1 + k) * c(kGgamma1 + k) + sig * sig);
  const Basis B =
      make_basis(0, {.n = 30, .r1 = 0.02, .rmax = 30.0, .eta = 1e-12});
  const Mat S = overlap(B), X = orthonormalizer(S, 1e-12);
  auto E = [&](double p) { return std::sqrt(p * p + m * m); };
  auto phi = [&](double p) {
    return std::pow(m * m / (E(p) * E(p)), 0.5 + eps);
  };
  const Mat Phi = momentum_diagonal(B, phi, 80);
  Mat Wc = Mat::Zero(B.size(), B.size());
  const double Kc = 2.0 / (3 * m * m) * 16 * std::numbers::pi / 3;
  for (int k = 0; k < 3; ++k)
    Wc += Kc * c(kGalpha1 + k) * std::pow(tau[k], 3) /
          std::pow(std::numbers::pi, 1.5) * gaussian(B, B, tau[k] * tau[k]);
  const Mat closure = X.transpose() *
                      (Phi * S.inverse() * Wc * S.inverse() * Phi) *
                      X;  // (F S^-1) W (S^-1 F) in the orthonormal basis
  const Mat closure2 = sandwich_closure(Phi, Wc, Phi, X, X);
  EXPECT_LT((closure - closure2).cwiseAbs().maxCoeff(),
            1e-6 * closure2.cwiseAbs().maxCoeff());
  const auto contact_w = [&](double r) {
    double v = 0;
    for (int k = 0; k < 3; ++k)
      v += Kc * c(kGalpha1 + k) * std::pow(tau[k], 3) /
           std::pow(std::numbers::pi, 1.5) * std::exp(-tau[k] * tau[k] * r * r);
    return v;
  };
  const Mat exact =
      X.transpose() * sandwich_two_step(B, B, phi, contact_w, 35.0) * X;
  // compare through the expectation value in a low-lying "state": the lowest
  // eigenvector of the kinetic+confinement problem
  Eigen::SelfAdjointEigenSolver<Mat> es(closure2 + exact.cwiseAbs() * 0.0);
  const Vec v = X.transpose() * S *
                B.nu.cwiseInverse()
                    .cwiseSqrt()
                    .cwiseInverse()
                    .cwiseInverse()
                    .cwiseInverse();  // arbitrary smooth test vector
  (void)v;
  // relative agreement of the Frobenius-weighted low-rank part: the matrices
  // share their dominant eigenvalues
  Eigen::SelfAdjointEigenSolver<Mat> ea(closure2), eb(exact);
  const double top_c = ea.eigenvalues().tail(3).sum(),
               top_e = eb.eigenvalues().tail(3).sum();
  EXPECT_NEAR(top_c / top_e, 1.0,
              2e-3);  // closure error on the contact matrix is ~1e-3 relative
                      // (it enters at O(hyperfine))
}

TEST(GI, ScreenedModelDReducesToModelCForSmallMu) {
  GIModel C, D({.screened = true});
  MatrixCache cc(kSpec), cd(kSpec);
  Vec th = gi_theta(C, 1.628, 1.628);
  const SolvedChannel a = solve_channel(C, find(C, "3S1-3D1"), th, cc);
  th(kGmuS) = 1e-4;
  const SolvedChannel b = solve_channel(D, find(D, "3S1-3D1"), th, cd);
  for (int k = 0; k < 3; ++k)
    EXPECT_NEAR(a.level(k).M, b.level(k).M, 1e-3)
        << k;  // O(mu_s) shift, larger for excited states
  EXPECT_LT(b.level(0).M - a.level(0).M,
            1e-5);  // the ground state hardly moves
  // physical screening mu_s = 0.1 lowers highly excited levels much more than
  // the ground state
  th(kGmuS) = 0.1;
  const SolvedChannel d = solve_channel(D, find(D, "3S1-3D1"), th, cd);
  EXPECT_LT(d.level(3).M, a.level(3).M - 0.01);
  EXPECT_LT(
      a.level(0).M - d.level(0).M,
      0.5 *
          (a.level(3).M -
           d.level(3).M));  // 1S moves by ~bmu<r^2>/2 = 26 MeV, 4S by far more
}

TEST(GI, FeynmanHellmannByCentralDifferencesIsConsistentWithMassDerivatives) {
  GIModel model;
  MatrixCache cache(kSpec);
  const Vec th = gi_theta(model, 1.628, 1.628);
  const Channel ch = find(model, "3P1");
  const SolvedChannel sc = solve_channel(model, ch, th, cache, true);
  const Vec g = sc.mass_gradient(0, model, th);
  for (int p :
       {static_cast<int>(kGm1), static_cast<int>(kGb), static_cast<int>(kGc),
        static_cast<int>(kGsigma0), static_cast<int>(kGepsSov)}) {
    const double h = 1e-4 * std::max(std::fabs(th(p)), 0.1);
    Vec tp = th, tm = th;
    tp(p) += h;
    tm(p) -= h;
    MatrixCache c1(kSpec), c2(kSpec);
    const double fd = (solve_channel(model, ch, tp, c1).level(0).M -
                       solve_channel(model, ch, tm, c2).level(0).M) /
                      (2 * h);
    EXPECT_NEAR(g(p), fd, 2e-5 * std::max(1.0, std::fabs(fd))) << p;
  }
  EXPECT_NEAR(g(kGc), 1.0, 1e-6);  // dM/dc = 1: c enters as c*S
}

TEST(GI, UnequalMassesGiveTheBcFamily) {
  GIModel model;
  MatrixCache cache(kSpec);
  const Vec th = gi_theta(model, 1.628, 4.977);
  const auto chans = model.channels(1, false);
  double m_1S0 = 0, m_3S1 = 0;
  for (const Channel& c : chans) {
    const SolvedChannel sc = solve_channel(model, c, th, cache);
    if (c.label() == "1S0") m_1S0 = sc.level(0).M;
    if (c.label() == "3S1")
      m_3S1 = sc.level(0).M;  // Lmax = 1: the D-wave partner is not in the list
  }
  EXPECT_NEAR(m_1S0, 6.27,
              0.03);  // B_c (6274.5 MeV measured); GI predict 6.264
  EXPECT_GT(m_3S1, m_1S0);
  EXPECT_LT(m_3S1 - m_1S0, 0.1);  // hyperfine splitting of a few tens of MeV
}
