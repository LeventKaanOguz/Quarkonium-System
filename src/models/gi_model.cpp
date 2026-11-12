#include "models/gi_model.hpp"

#include <algorithm>
#include <cmath>
#include <map>
#include <numbers>
#include <stdexcept>
#include <tuple>

#include "linalg/gep.hpp"
#include "me/analytic.hpp"
#include "me/coulomb.hpp"
#include "me/momentum.hpp"
#include "me/radial.hpp"
#include "me/smeared.hpp"
#include "models/alpha_s.hpp"

namespace qgem {
namespace {
constexpr double kPi = std::numbers::pi;
}

std::vector<std::string> GIModel::param_names() const {
  return {"m1",     "m2",     "b",       "c",       "sigma0", "s",
          "eps_c",  "eps_t",  "eps_sov", "eps_sos", "alpha1", "alpha2",
          "alpha3", "gamma1", "gamma2",  "gamma3",  "mu_s"};
}

std::vector<Channel> GIModel::channels(int Lmax, bool equal_masses) const {
  return enumerate_channels(Lmax, equal_masses);
}

double GIModel::sigma12(const Vec& th) {
  const double m1 = th(kGm1), m2 = th(kGm2), s0 = th(kGsigma0), s = th(kGs);
  const double r = 4 * m1 * m2 / ((m1 + m2) * (m1 + m2));
  return std::sqrt(s0 * s0 * (0.5 + 0.5 * std::pow(r, 4)) +
                   s * s * std::pow(2 * m1 * m2 / (m1 + m2), 2));
}

ParamSet GIModel::godfrey_isgur_params(double m1, double m2, GiRunning running,
                                       bool with_c) const {
  const auto a = gi_alpha(), g = gi_gamma();
  auto spec = [](const char* name, double v, bool free,
                 Transform tr = Transform::Identity) {
    ParamSpec p;
    p.name = name;
    p.value = v;
    p.free = free;
    p.tr = tr;
    return p;
  };
  ParamSpec mass2 = spec("m2", m2, false);
  // physical bounds: a smearing width below ~0.1 GeV makes the contact/tensor
  // terms singular again
  auto bounded = [](ParamSpec p, double lo) {
    p.lo = lo;
    return p;
  };
  if (m1 == m2) mass2.tied_to = "m1";
  std::vector<ParamSpec> p = {
      spec("m1", m1, true, Transform::Log),
      mass2,
      spec("b", 0.18, true, Transform::Log),
      spec("c", with_c ? -0.253 : 0.0, with_c),
      bounded(spec("sigma0", 1.80, true, Transform::Log), 0.1),
      spec("s", 1.55, false),
      spec("eps_c", -0.168, false),
      spec("eps_t", 0.025, false),
      spec("eps_sov", -0.035, false),
      spec("eps_sos", 0.055, false)};
  const bool free_alpha = running == GiRunning::Free;
  if (running == GiRunning::Constant) {
    // one coupling: alpha_1 (free) with gamma_1 -> infinity; the other two
    // Gaussians switched off
    p.push_back(spec("alpha1", 0.5, true, Transform::Log));
    p.push_back(spec("alpha2", 0.0, false));
    p.push_back(spec("alpha3", 0.0, false));
    p.push_back(spec("gamma1", 1e4, false));
    p.push_back(spec("gamma2", g[1], false));
    p.push_back(spec("gamma3", g[2], false));
  } else {
    for (int k = 0; k < 3; ++k)
      p.push_back(spec(("alpha" + std::to_string(k + 1)).c_str(), a[k],
                       free_alpha, Transform::Log));
    for (int k = 0; k < 3; ++k)
      p.push_back(spec(("gamma" + std::to_string(k + 1)).c_str(), g[k], false));
  }
  // screening mass: bounded at 1e-4 GeV (range 2000 fm: no screening, model
  // C) so finite-difference steps stay above the 1e-5 limit of build()
  p.push_back(bounded(spec("mu_s", 0.1, opt_.screened, Transform::Log), 1e-4));
  return ParamSet(std::move(p));
}

ChannelMatrices GIModel::assemble(const Channel& ch, const Vec& th,
                                  MatrixCache& cache,
                                  MomentumMemo& memo) const {
  if (th.size() != kGIParams)
    throw std::invalid_argument("GIModel: theta has the wrong size");
  const double m1 = th(kGm1), m2 = th(kGm2), b = th(kGb), c = th(kGc);
  if (!(m1 > 0 && m2 > 0 && th(kGsigma0) > 0))
    throw std::invalid_argument("GIModel: masses and sigma0 must be > 0");
  if (opt_.screened && !(th(kGmuS) >= 1e-5))
    throw std::invalid_argument("GIModel (D): mu_s >= 1e-5");
  const double sig = sigma12(th), mu = th(kGmuS);
  const double alpha[3] = {th(kGalpha1), th(kGalpha2), th(kGalpha3)},
               gamma[3] = {th(kGgamma1), th(kGgamma2), th(kGgamma3)};
  double tau[3];
  for (int k = 0; k < 3; ++k)
    tau[k] = sig * gamma[k] / std::sqrt(gamma[k] * gamma[k] + sig * sig);
  const double eps[4] = {th(kGepsC), th(kGepsT), th(kGepsSov),
                         th(kGepsSos)};  // c, t, sov, sos

  auto Ea = [](double m, double p) { return std::sqrt(p * p + m * m); };
  const int nq = opt_.nq;

  // exact analytic matrices of the smeared radial functions (G~'/r, tensor,
  // S~'/r)
  auto Gp_matrix = [&](const Basis& bb, const Basis& kk) {
    Mat M = Mat::Zero(bb.size(), kk.size());
    for (int k = 0; k < 3; ++k)
      M += -4.0 / 3.0 * alpha[k] * gprime_over_r_matrix(bb, kk, tau[k]);
    return M;
  };
  auto T_matrix = [&](const Basis& bb, const Basis& kk) {
    Mat M = Mat::Zero(bb.size(), kk.size());
    for (int k = 0; k < 3; ++k)
      M += -4.0 / 3.0 * alpha[k] * tensor_matrix(bb, kk, tau[k]);
    return M;
  };
  auto Sp_matrix = [&](const Basis& bb) -> Mat {
    if (!opt_.screened) return b * smeared_r_prime_over_r_matrix(bb, bb, sig);
    return -b / mu * smeared_exp_prime_over_r_matrix(bb, sig, mu);
  };

  ChannelMatrices out;
  std::vector<int> off;
  std::vector<Mat> Xs;
  int N = 0, Np = 0;
  for (const Component& cp : ch.comps) {
    const int n = cache.basis(cp.L).size();
    out.L.push_back(cp.L);
    out.block_size.push_back(n);
    off.push_back(N);
    N += n;
    Xs.push_back(orthonormalizer(cache.overlap(cp.L), cache.spec(cp.L).eta));
    Np += static_cast<int>(Xs.back().cols());
  }
  std::vector<int> offp(ch.comps.size(), 0);
  for (std::size_t i = 1; i < ch.comps.size(); ++i)
    offp[i] = offp[i - 1] + static_cast<int>(Xs[i - 1].cols());
  Mat Hp = Mat::Zero(Np, Np);
  out.S = Mat::Zero(N, N);
  out.X = Mat::Zero(N, Np);
  for (std::size_t i = 0; i < ch.comps.size(); ++i) {
    out.S.block(off[i], off[i], out.block_size[i], out.block_size[i]) =
        cache.overlap(ch.comps[i].L);
    out.X.block(off[i], offp[i], Xs[i].rows(), Xs[i].cols()) = Xs[i];
  }

  // per-component primed matrices
  struct Prim {
    Mat K, FG, G, Sconf, Phi_c, Phi_t, Phi_v11, Phi_v22, Phi_v12, Phi_s11,
        Phi_s22, Gp, Sp;
  };
  // momentum-diagonal matrices, memoised (all are symmetric in ma <-> mb)
  enum Kind { kPhi, kKinetic, kFG };
  auto momentum = [&](Kind kind, const Basis& B, double ma, double mb,
                      double e) -> const Mat& {
    const auto key = std::make_tuple(static_cast<int>(kind), B.L,
                                     std::min(ma, mb), std::max(ma, mb), e);
    auto it = memo.find(key);
    if (it != memo.end()) return it->second;
    Mat M;
    switch (kind) {
      case kPhi:
        M = momentum_diagonal(
            B,
            [=](double p) {
              return std::pow(ma * mb / (Ea(ma, p) * Ea(mb, p)), 0.5 + e);
            },
            nq);
        break;
      case kKinetic:
        M = momentum_diagonal(
            B, [=](double p) { return Ea(ma, p) + Ea(mb, p); }, nq);
        break;
      case kFG:
        M = momentum_diagonal(
            B,
            [=](double p) {
              return std::sqrt(1.0 + p * p / (Ea(ma, p) * Ea(mb, p)));
            },
            nq);
        break;
    }
    return memo.emplace(key, std::move(M)).first->second;
  };
  auto phi = [&](const Basis& B, double ma, double mb, double e) -> const Mat& {
    return momentum(kPhi, B, ma, mb, e);
  };
  std::vector<Prim> P(ch.comps.size());
  for (std::size_t i = 0; i < ch.comps.size(); ++i) {
    const int L = ch.comps[i].L;
    const Basis& B = cache.basis(L);
    const Mat& X = Xs[i];
    auto pr = [&](const Mat& A) { return Mat(X.transpose() * A * X); };
    P[i].K = pr(momentum(kKinetic, B, m1, m2, 0.0));
    P[i].FG = pr(momentum(kFG, B, m1, m2, 0.0));
    Mat G = Mat::Zero(B.size(), B.size());
    for (int k = 0; k < 3; ++k)
      G += -4.0 / 3.0 * alpha[k] * erf_over_r(B, B, tau[k]);
    P[i].G = pr(G);
    Mat Sc = c * cache.overlap(L);
    if (!opt_.screened) {
      Sc += b * smeared_r_matrix(B, B, sig);
    } else {
      Sc += b / mu * (cache.overlap(L) - smeared_exp_matrix(B, sig, mu));
    }
    P[i].Sconf = pr(Sc);
    P[i].Phi_c = pr(phi(B, m1, m2, eps[0]));
    P[i].Phi_t = pr(phi(B, m1, m2, eps[1]));
    P[i].Phi_v11 = pr(phi(B, m1, m1, eps[2]));
    P[i].Phi_v22 = pr(phi(B, m2, m2, eps[2]));
    P[i].Phi_v12 = pr(phi(B, m1, m2, eps[2]));
    P[i].Phi_s11 = pr(phi(B, m1, m1, eps[3]));
    P[i].Phi_s22 = pr(phi(B, m2, m2, eps[3]));
    if (L >= 1) {
      P[i].Gp = pr(Gp_matrix(B, B));
      P[i].Sp = pr(Sp_matrix(B));
    }
    Hp.block(offp[i], offp[i], X.cols(), X.cols()) +=
        P[i].K + P[i].FG * P[i].G * P[i].FG + P[i].Sconf;
  }

  // spin-orbit coefficient matrices A (coefficient of L.S) and B (coefficient
  // of L.S_-), in the primed basis of block i
  auto so_AB = [&](std::size_t i, Mat& A, Mat& Bm) {
    const Prim& q = P[i];
    // a_k = Phi_kk^{sov} Gp/(2 m_k^2) Phi_kk^{sov} + Phi_12^{sov} Gp/(m1 m2)
    // Phi_12^{sov} - Phi_kk^{sos} Sp/(2 m_k^2) Phi_kk^{sos}
    const Mat a1 = q.Phi_v11 * q.Gp * q.Phi_v11 / (2 * m1 * m1) +
                   q.Phi_v12 * q.Gp * q.Phi_v12 / (m1 * m2) -
                   q.Phi_s11 * q.Sp * q.Phi_s11 / (2 * m1 * m1);
    const Mat a2 = q.Phi_v22 * q.Gp * q.Phi_v22 / (2 * m2 * m2) +
                   q.Phi_v12 * q.Gp * q.Phi_v12 / (m1 * m2) -
                   q.Phi_s22 * q.Sp * q.Phi_s22 / (2 * m2 * m2);
    A = 0.5 * (a1 + a2);
    Bm = 0.5 * (a1 - a2);
  };

  const double Kc = 2.0 / (3.0 * m1 * m2) * 16.0 * kPi / 3.0;
  for (std::size_t i = 0; i < ch.comps.size(); ++i) {
    const int L = ch.comps[i].L, S = ch.comps[i].S, J = ch.J;
    const Basis& B = cache.basis(L);
    const Mat& X = Xs[i];
    const int o = offp[i], nk = static_cast<int>(X.cols());
    // contact: Phi_c [ (2/3 m1 m2) lap G~ ] Phi_c S1.S2,  lap G~ = (16 pi/3)
    // sum a_k rho_tau_k
    Mat Wc = Mat::Zero(B.size(), B.size());
    for (int k = 0; k < 3; ++k)
      Wc += Kc * alpha[k] * std::pow(tau[k], 3) / std::pow(kPi, 1.5) *
            gaussian(B, B, tau[k] * tau[k]);
    Hp.block(o, o, nk, nk) +=
        s1s2(S) * P[i].Phi_c * (X.transpose() * Wc * X) * P[i].Phi_c;
    if (L >= 1 && S == 1) {
      Mat Wt = T_matrix(B, B);
      Hp.block(o, o, nk, nk) += s12_diag(L, J) * (1.0 / (12 * m1 * m2)) *
                                P[i].Phi_t * (X.transpose() * Wt * X) *
                                P[i].Phi_t;
      Mat A, Bm;
      so_AB(i, A, Bm);
      Hp.block(o, o, nk, nk) += ls(L, S, J) * A;
    }
  }
  if (ch.comps.size() == 2) {
    const Component &c1 = ch.comps[0], &c2 = ch.comps[1];
    const Basis &B1 = cache.basis(c1.L), &B2 = cache.basis(c2.L);
    Mat off12 = Mat::Zero(Xs[0].cols(), Xs[1].cols());
    if (c1.S == 1 && c2.S == 1 && c1.L != c2.L) {
      const Mat Wt12 = T_matrix(B1, B2);
      off12 = s12_offdiag(ch.J) * (1.0 / (12 * m1 * m2)) * P[0].Phi_t *
              (Xs[0].transpose() * Wt12 * Xs[1]) * P[1].Phi_t;
    } else if (c1.L == c2.L && c1.S == 0 && c2.S == 1 && c1.L >= 1) {
      Mat A, Bm;
      so_AB(0, A, Bm);
      off12 = ls_minus(c1.L) * Bm;
    }
    Hp.block(offp[0], offp[1], Xs[0].cols(), Xs[1].cols()) += off12;
    Hp.block(offp[1], offp[0], Xs[1].cols(), Xs[0].cols()) += off12.transpose();
  }
  Hp = 0.5 * (Hp + Hp.transpose());
  // H in the original basis such that X^T H X = H' and H = S X H' X^T S
  out.H = out.S * out.X * Hp * out.X.transpose() * out.S;
  return out;
}

ChannelMatrices GIModel::build(const Channel& ch, const Vec& theta,
                               MatrixCache& cache, bool derivatives) const {
  // Resolution: the smeared terms must be wider than the narrowest Gaussian
  // (see NRModel::build).
  for (const auto& comp : ch.comps)
    if (sigma12(theta) * cache.spec(comp.L).r1 > 1.0)
      throw std::invalid_argument(
          "GIModel::build: smearing 1/sigma12 below the basis resolution r1");
  MomentumMemo memo;
  ChannelMatrices out = assemble(ch, theta, cache, memo);
  if (!derivatives) return out;
  out.n_dp = kGIParams;
  for (int a = 0; a < kGIParams; ++a) {
    if (!active_.empty() &&
        std::find(active_.begin(), active_.end(), a) == active_.end()) {
      out.dH[a] = Mat::Zero(out.H.rows(), out.H.cols());
      continue;
    }
    const double h = 1e-5 * std::max(std::fabs(theta(a)), 0.1);
    Vec tp = theta, tm = theta;
    tp(a) += h;
    tm(a) -= h;
    out.dH[a] =
        (assemble(ch, tp, cache, memo).H - assemble(ch, tm, cache, memo).H) /
        (2 * h);
  }
  return out;
}

}  // namespace qgem
