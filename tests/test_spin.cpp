// M2: spin-angle factors. Independent check: the operators are built explicitly
// as matrices in the |L M_L; S M_S> product basis, with Clebsch-Gordan
// coefficients from the Racah formula, then the coupled states |(L S) J M> are
// formed and <S1.S2>, <L.S>, <S12> evaluated numerically.
#include <gtest/gtest.h>

#include <cmath>
#include <complex>
#include <vector>

#include "core/types.hpp"
#include "spin/spin.hpp"

using namespace qgem;

namespace {
double fact(int n) { return std::tgamma(n + 1.0); }

// <j1 m1 j2 m2 | J M>, Racah formula, integer arguments.
double cg(int j1, int m1, int j2, int m2, int J, int M) {
  if (m1 + m2 != M || J < std::abs(j1 - j2) || J > j1 + j2 ||
      std::abs(m1) > j1 || std::abs(m2) > j2 || std::abs(M) > J)
    return 0.0;
  const double pre =
      std::sqrt((2.0 * J + 1) * fact(J + j1 - j2) * fact(J - j1 + j2) *
                fact(j1 + j2 - J) / fact(j1 + j2 + J + 1)) *
      std::sqrt(fact(J + M) * fact(J - M) * fact(j1 - m1) * fact(j1 + m1) *
                fact(j2 - m2) * fact(j2 + m2));
  double sum = 0;
  for (int k = 0; k <= j1 + j2 + J; ++k) {
    const int a = j1 + j2 - J - k, b = j1 - m1 - k, c = j2 + m2 - k,
              d = J - j2 + m1 + k, e = J - j1 - m2 + k;
    if (a < 0 || b < 0 || c < 0 || d < 0 || e < 0) continue;
    sum += ((k % 2) ? -1.0 : 1.0) /
           (fact(k) * fact(a) * fact(b) * fact(c) * fact(d) * fact(e));
  }
  return pre * sum;
}

using CM = Eigen::MatrixXcd;
using C = std::complex<double>;

// Spin-1/2 operators; two-spin space ordered (m1, m2) with m = +1/2, -1/2.
CM sigma(int axis) {
  CM s(2, 2);
  if (axis == 0) s << 0, 1, 1, 0;
  if (axis == 1) s << 0, C(0, -1), C(0, 1), 0;
  if (axis == 2) s << 1, 0, 0, -1;
  return s;
}
CM kron(const CM& a, const CM& b) {
  CM k(a.rows() * b.rows(), a.cols() * b.cols());
  for (int i = 0; i < a.rows(); ++i)
    for (int j = 0; j < a.cols(); ++j)
      k.block(i * b.rows(), j * b.cols(), b.rows(), b.cols()) = a(i, j) * b;
  return k;
}

// Orbital angular momentum matrices for a given L in the |L m> basis, m =
// L..-L.
std::vector<CM> lmat(int L) {
  const int n = 2 * L + 1;
  CM lz = CM::Zero(n, n), lp = CM::Zero(n, n);
  for (int i = 0; i < n; ++i) {
    const int m = L - i;
    lz(i, i) = m;
    if (i > 0)
      lp(i - 1, i) = std::sqrt(
          static_cast<double>(L * (L + 1) - m * (m + 1)));  // |m+1><m|
  }
  CM lm = lp.adjoint();
  return {0.5 * (lp + lm), (lp - lm) / C(0, 2), lz};
}

// r-hat components projected on fixed L are zero (parity); S12 needs the full
// space {L, L'}. We build S12 from the angular matrix elements <L' m'| r_i r_j
// |L m> by Gaunt integration done via numerical quadrature over the sphere with
// explicit spherical harmonics (Condon-Shortley).
C ylm(int l, int m, double th, double ph) {
  // associated Legendre via recurrence (no CS phase in P_l^m for m>=0 handled
  // below)
  const int am = std::abs(m);
  double pmm = 1.0;
  const double x = std::cos(th), s = std::sin(th);
  for (int i = 1; i <= am; ++i)
    pmm *= -(2.0 * i - 1) * s;  // includes the Condon-Shortley phase
  double p = pmm;
  if (l > am) {
    double pm1 = x * (2 * am + 1) * pmm;
    if (l == am + 1)
      p = pm1;
    else {
      double pl2 = pmm, pl1 = pm1;
      for (int ll = am + 2; ll <= l; ++ll) {
        p = (x * (2 * ll - 1) * pl1 - (ll + am - 1) * pl2) / (ll - am);
        pl2 = pl1;
        pl1 = p;
      }
    }
  }
  const double norm =
      std::sqrt((2.0 * l + 1) / (4 * M_PI) * fact(l - am) / fact(l + am));
  C y = norm * p * std::polar(1.0, am * ph);
  if (m < 0)
    y = ((am % 2) ? -1.0 : 1.0) * std::conj(y);  // Y_{l,-m} = (-1)^m conj(Y_lm)
  return y;
}

// <L' M'| rhat_i rhat_j |L M> by Gauss-Legendre (theta) x uniform (phi)
// quadrature
double nhat(int i, double th, double ph) {
  return i == 0 ? std::sin(th) * std::cos(ph)
                : (i == 1 ? std::sin(th) * std::sin(ph) : std::cos(th));
}

}  // namespace

TEST(Spin, ClosedFormFactors) {
  EXPECT_DOUBLE_EQ(s1s2(0), -0.75);
  EXPECT_DOUBLE_EQ(s1s2(1), 0.25);
  // 3P2, 3P1, 3P0: <L.S> = 1, -1, -2 ; <S12> = -2/5, 2, -4
  EXPECT_DOUBLE_EQ(ls(1, 1, 2), 1.0);
  EXPECT_DOUBLE_EQ(ls(1, 1, 1), -1.0);
  EXPECT_DOUBLE_EQ(ls(1, 1, 0), -2.0);
  EXPECT_NEAR(s12_diag(1, 2), -0.4, 1e-15);
  EXPECT_DOUBLE_EQ(s12_diag(1, 1), 2.0);
  EXPECT_DOUBLE_EQ(s12_diag(1, 0), -4.0);
  EXPECT_DOUBLE_EQ(s12_diag(0, 1), 0.0);
  EXPECT_NEAR(s12_offdiag(1), 2 * std::sqrt(2.0), 1e-15);
  EXPECT_NEAR(ls_minus(1), std::sqrt(2.0), 1e-15);
}

TEST(Spin, CentreOfGravityTheorem) {
  // sum_J (2J+1) <L.S> = 0 and sum_J (2J+1) <S12> = 0 for fixed L >= 1, S = 1.
  for (int L = 1; L <= 6; ++L) {
    double a = 0, t = 0;
    for (int J = L - 1; J <= L + 1; ++J) {
      a += (2 * J + 1) * ls(L, 1, J);
      t += (2 * J + 1) * s12_diag(L, J);
    }
    EXPECT_NEAR(a, 0.0, 1e-12) << L;
    EXPECT_NEAR(t, 0.0, 1e-12) << L;
  }
}

TEST(Spin, OperatorsBuiltFromClebschGordanMatchClosedForms) {
  // Spin space (4-dim), S1.S2 and S12 acting on spin only need r-hat: S12 = 2[3
  // (S.n)^2 - S^2] (notes 3.24).
  CM s1[3], s2[3];
  const CM id = CM::Identity(2, 2);
  for (int a = 0; a < 3; ++a) {
    s1[a] = 0.5 * kron(sigma(a), id);
    s2[a] = 0.5 * kron(id, sigma(a));
  }
  CM S[3];
  for (int a = 0; a < 3; ++a) S[a] = s1[a] + s2[a];
  CM S1S2 = CM::Zero(4, 4);
  for (int a = 0; a < 3; ++a) S1S2 += s1[a] * s2[a];

  const int Ls[] = {0, 1, 2, 3};
  for (int L : Ls) {
    const int nL = 2 * L + 1;
    const std::vector<CM> Lop = lmat(L);
    // coupled states |(L,S=1) J M> in the product basis (index = iL * 4 + iS);
    // spin triplet states in the (m1,m2) basis: |1,1>=uu, |1,0>=(ud+du)/sqrt2,
    // |1,-1>=dd.
    auto triplet = [&](int mS) {
      Eigen::VectorXcd v =
          Eigen::VectorXcd::Zero(4);  // basis order: uu, ud, du, dd
      if (mS == 1) v(0) = 1;
      if (mS == 0) v(1) = v(2) = 1 / std::sqrt(2.0);
      if (mS == -1) v(3) = 1;
      return v;
    };
    // L.S in the product space
    CM LS = CM::Zero(nL * 4, nL * 4);
    for (int a = 0; a < 3; ++a) LS += kron(Lop[a], S[a]);
    for (int J = std::abs(L - 1); J <= L + 1; ++J) {
      if (L == 0 && J != 1) continue;
      const int M = J;  // any M; take the stretched one
      Eigen::VectorXcd psi = Eigen::VectorXcd::Zero(nL * 4);
      for (int mS = -1; mS <= 1; ++mS) {
        const int mL = M - mS;
        if (std::abs(mL) > L) continue;
        const double c = cg(L, mL, 1, mS, J, M);
        const int iL = L - mL;
        psi.segment(iL * 4, 4) += c * triplet(mS);
      }
      ASSERT_NEAR(psi.norm(), 1.0, 1e-12) << "L=" << L << " J=" << J;
      EXPECT_NEAR((psi.adjoint() * LS * psi)(0).real(), ls(L, 1, J), 1e-12)
          << "L=" << L << " J=" << J;
      // <S1.S2> = 1/4 for any triplet state
      CM full = kron(CM::Identity(nL, nL), S1S2);
      EXPECT_NEAR((psi.adjoint() * full * psi)(0).real(), s1s2(1), 1e-12);
    }
  }

  // Tensor operator: S12 = 2[3 (S.n)^2 - S^2] with n = r-hat, matrix elements
  // between the coupled states
  // |(L S=1) J M> and |(L' S=1) J M>, angular integrals by quadrature over the
  // sphere.
  const int nth = 40, nph = 64;
  std::vector<long double> gx, gw;  // Gauss-Legendre in cos(theta)
  {
    gx.assign(nth, 0);
    gw.assign(nth, 0);
    for (int i = 0; i < nth; ++i) {
      long double z = std::cos(M_PI * (i + 0.75L) / (nth + 0.5L)), dp = 1;
      for (int it = 0; it < 100; ++it) {
        long double p0 = 1, p1 = z;
        for (int k = 2; k <= nth; ++k) {
          long double p2 = ((2 * k - 1) * z * p1 - (k - 1) * p0) / k;
          p0 = p1;
          p1 = p2;
        }
        dp = nth * (z * p1 - p0) / (z * z - 1);
        const long double dz = p1 / dp;
        z -= dz;
        if (std::fabs(dz) < 1e-18L) break;
      }
      gx[i] = z;
      gw[i] = 2 / ((1 - z * z) * dp * dp);
    }
  }
  auto s12_elem = [&](int L1, int L2, int J) {
    const int M = J;
    auto spinor = [&](int L, double th, double ph) {
      Eigen::VectorXcd v = Eigen::VectorXcd::Zero(4);
      for (int mS = -1; mS <= 1; ++mS) {
        const int mL = M - mS;
        if (std::abs(mL) > L) continue;
        const double c = cg(L, mL, 1, mS, J, M);
        Eigen::VectorXcd t = Eigen::VectorXcd::Zero(4);
        if (mS == 1) t(0) = 1;
        if (mS == 0) t(1) = t(2) = 1 / std::sqrt(2.0);
        if (mS == -1) t(3) = 1;
        v += c * ylm(L, mL, th, ph) * t;
      }
      return v;
    };
    C sum = 0;
    for (int i = 0; i < nth; ++i)
      for (int k = 0; k < nph; ++k) {
        const double th = std::acos(static_cast<double>(gx[i])),
                     ph = 2 * M_PI * (k + 0.5) / nph;
        CM Sn = CM::Zero(4, 4);
        for (int a = 0; a < 3; ++a) Sn += nhat(a, th, ph) * S[a];
        CM S2 = CM::Zero(4, 4);
        for (int a = 0; a < 3; ++a) S2 += S[a] * S[a];
        const CM S12 = 2.0 * (3.0 * Sn * Sn - S2);
        sum += static_cast<double>(gw[i]) * (2 * M_PI / nph) *
               (spinor(L1, th, ph).adjoint() * S12 * spinor(L2, th, ph))(0);
      }
    return sum.real();
  };
  for (int L = 1; L <= 3; ++L)
    for (int J = L - 1; J <= L + 1; ++J)
      EXPECT_NEAR(s12_elem(L, L, J), s12_diag(L, J), 1e-10)
          << "L=" << L << " J=" << J;
  EXPECT_NEAR(s12_elem(0, 2, 1), s12_offdiag(1), 1e-10);  // 3S1-3D1: 2 sqrt 2
  EXPECT_NEAR(s12_elem(1, 3, 2), s12_offdiag(2), 1e-10);  // 3P2-3F2
  EXPECT_NEAR(s12_elem(2, 4, 3), s12_offdiag(3), 1e-10);  // 3D3-3G3
}

TEST(Spin, ChannelEnumeration) {
  const auto eq = enumerate_channels(2, true);
  std::vector<std::string> names;
  for (const auto& c : eq) names.push_back(c.label());
  auto has = [&](const std::string& s) {
    return std::find(names.begin(), names.end(), s) != names.end();
  };
  EXPECT_TRUE(has("1S0"));
  EXPECT_TRUE(has("3S1-3D1"));
  EXPECT_TRUE(has("1P1"));
  EXPECT_TRUE(has("3P0"));
  EXPECT_TRUE(has("3P1"));
  EXPECT_TRUE(has("3P2"));  // L = 3 > Lmax: 3P2 stays a single component
  EXPECT_TRUE(has("1D2"));
  // quantum numbers: J/psi 1--, chi_c1 1++, h_c 1+-, eta_c 0-+
  for (const auto& c : eq) {
    if (c.label() == "3S1-3D1") {
      EXPECT_EQ(c.P, -1);
      EXPECT_EQ(c.C, -1);
    }
    if (c.label() == "3P1") {
      EXPECT_EQ(c.P, 1);
      EXPECT_EQ(c.C, 1);
    }
    if (c.label() == "1P1") {
      EXPECT_EQ(c.P, 1);
      EXPECT_EQ(c.C, -1);
    }
    if (c.label() == "1S0") {
      EXPECT_EQ(c.P, -1);
      EXPECT_EQ(c.C, 1);
    }
  }
  // no invalid (L = 0, S = 1, J = 0) channel, and J = 0 is the singlet only for
  // unequal masses
  for (const auto& c : eq) EXPECT_NE(c.label(), "3S0");
  const auto uneq = enumerate_channels(1, false);
  for (const auto& c : uneq) EXPECT_NE(c.label(), "1S0-3S0");
  bool found = false;
  for (const auto& c : uneq)
    if (c.label() == "1P1-3P1") found = true;
  EXPECT_TRUE(found);
}
