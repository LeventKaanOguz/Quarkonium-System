#include "models/nr_model.hpp"

#include <cmath>
#include <numbers>
#include <stdexcept>

#include "me/analytic.hpp"
#include "me/momentum.hpp"
#include "models/alpha_s.hpp"

namespace qgem {
namespace {

constexpr double kSqrtPi = 1.7724538509055160273;

/// Accumulates  H += coef * M  together with dH_a += dcoef_a * M  (+ coef *
/// dscale * dM when M itself depends on parameter p), for one block position.
struct Block {
  Mat* H;
  std::array<Mat, kMaxP>* dH;
  int r0, c0;
  int nd;
  void add(const Dual& coef, const Mat& M, int p = -1, const Mat* dM = nullptr,
           double dscale = 1.0) const {
    H->block(r0, c0, M.rows(), M.cols()) += coef.v * M;
    for (int a = 0; a < nd; ++a)
      if (coef.d[a] != 0.0)
        (*dH)[a].block(r0, c0, M.rows(), M.cols()) += coef.d[a] * M;
    if (p >= 0 && dM != nullptr && nd > p)
      (*dH)[p].block(r0, c0, M.rows(), M.cols()) += (coef.v * dscale) * (*dM);
  }
};

/// One term of the coupling: weight w, central range gamma (infinite for a
/// constant coupling) and the spin-term width tau(sigma) = sigma gamma /
/// sqrt(gamma^2 + sigma^2) (= sigma for a constant coupling).
struct Term {
  double w, gamma;
  bool point;
};

}  // namespace

std::string NRModel::name() const {
  std::string s = opt_.screened ? "B" : "A";
  if (opt_.running == Running::Gaussian)
    s += " (running alpha_s, " +
         (opt_.shape_label.empty() ? std::string("Godfrey-Isgur")
                                   : opt_.shape_label) +
         " shape)";
  else
    s += " (constant alpha_s)";
  switch (opt_.spin) {
    case SpinMode::None:
      s += "/spin-independent";
      break;
    case SpinMode::FirstOrder:
      s += "/first-order LS+T";
      break;
    case SpinMode::Smeared:
      s += "/smeared LS+T in H";
      break;
  }
  if (opt_.relativistic_kinetic) s += "/relativistic kinetic";
  return s;
}

std::vector<std::string> NRModel::param_names() const {
  return {"m1", "m2", "alpha_s", "b", "c", "sigma", "mu_s"};
}

Vec NRModel::mass_offset_gradient(const Vec&) const {
  Vec g = Vec::Zero(kNRParams);
  g(kM1) = 1.0;
  g(kM2) = 1.0;
  return g;
}

ParamSet NRModel::default_params(double m1, double m2, double alpha_s, double b,
                                 double c, double sigma, double mu_s) const {
  auto spec = [](const char* name, double v, bool free, Transform tr) {
    ParamSpec p;
    p.name = name;
    p.value = v;
    p.free = free;
    p.tr = tr;
    return p;
  };
  // screening mass: bounded at 1e-4 GeV (range 2000 fm: no screening, model
  // A), above the 1e-5 limit of build()
  ParamSpec screening = spec("mu_s", mu_s, opt_.screened, Transform::Log);
  screening.lo = 1e-4;
  return ParamSet({spec("m1", m1, true, Transform::Log),
                   spec("m2", m2, true, Transform::Log),
                   spec("alpha_s", alpha_s, true, Transform::Log),
                   spec("b", b, true, Transform::Log),
                   spec("c", c, false, Transform::Identity),
                   spec("sigma", sigma, true, Transform::Log), screening});
}

std::vector<Channel> NRModel::channels(int Lmax, bool equal_masses) const {
  if (opt_.spin == SpinMode::Smeared)
    return enumerate_channels(Lmax, equal_masses);
  // Uncoupled pieces: every (L, S, J) is its own channel.
  std::vector<Channel> out;
  for (int L = 0; L <= Lmax; ++L)
    for (int S = 0; S <= 1; ++S)
      for (int J = std::abs(L - S); J <= L + S; ++J) {
        Channel ch;
        ch.J = J;
        ch.P = (L % 2 == 0) ? -1 : 1;
        ch.C = equal_masses ? ((L + S) % 2 == 0 ? 1 : -1) : 0;
        ch.comps = {{L, S}};
        out.push_back(ch);
      }
  return out;
}

ChannelMatrices NRModel::build(const Channel& ch, const Vec& th,
                               MatrixCache& cache, bool derivatives) const {
  if (th.size() != kNRParams)
    throw std::invalid_argument("NRModel::build: theta has the wrong size");
  if (!(th(kM1) > 0 && th(kM2) > 0 && th(kSigma) > 0 && th(kAlphaS) > 0))
    throw std::invalid_argument(
        "NRModel::build: need m1, m2, sigma, alpha_s > 0");
  if (opt_.screened && !(th(kMuS) >= 1e-5))
    throw std::invalid_argument("NRModel::build: Model B needs mu_s >= 1e-5");
  if (opt_.relativistic_kinetic && opt_.running == Running::Constant)
    throw std::invalid_argument(
        "NRModel::build: the relativistic kinetic energy needs a smeared "
        "(running) Coulomb potential: sqrt(p^2+m^2) - kappa/r is unbounded "
        "below for kappa > 2/pi");
  // Resolution: a smeared contact inside H must be wider than the narrowest
  // Gaussian, else its energy is set by r1, not by the physics (the fit would
  // then tune the basis). In first order a point contact is legitimate.
  if (opt_.spin == SpinMode::Smeared)
    for (const auto& comp : ch.comps)
      if (th(kSigma) * cache.spec(comp.L).r1 > 1.0)
        throw std::invalid_argument(
            "NRModel::build: contact smearing 1/sigma below the basis "
            "resolution r1");

  const int nd = derivatives ? kNRParams : 0;
  const Dual m1 = Dual::var(th(kM1), kM1), m2 = Dual::var(th(kM2), kM2),
             als = Dual::var(th(kAlphaS), kAlphaS), b = Dual::var(th(kB), kB),
             c = Dual::var(th(kC), kC), sig = Dual::var(th(kSigma), kSigma),
             mu = Dual::var(th(kMuS), kMuS);
  const double sg = th(kSigma), mus = th(kMuS);

  const Dual inv_mr = Dual(1.0) / m1 + Dual(1.0) / m2;
  const Dual kappa = Dual(4.0 / 3.0) * als;
  const Dual cV =
      Dual(0.25) / (m1 * m1) + Dual(0.25) / (m2 * m2) + Dual(1.0) / (m1 * m2);
  const Dual cS = Dual(0.25) / (m1 * m1) + Dual(0.25) / (m2 * m2);
  const Dual Bco = Dual(0.25) * (Dual(1.0) / (m1 * m1) -
                                 Dual(1.0) / (m2 * m2));  // antisymmetric LS

  // coupling terms and the spin-term widths tau_k(sigma) with dtau/dsigma
  std::vector<Term> terms;
  if (opt_.running == Running::Gaussian && !opt_.shape_w.empty()) {
    if (opt_.shape_w.size() != opt_.shape_gamma.size())
      throw std::invalid_argument("NRModel: shape_w and shape_gamma differ");
    for (std::size_t k = 0; k < opt_.shape_w.size(); ++k)
      if (opt_.shape_w[k] > 0)
        terms.push_back({opt_.shape_w[k], opt_.shape_gamma[k], false});
  } else if (opt_.running == Running::Gaussian) {
    const auto a = gi_alpha(), g = gi_gamma();
    double tot = 0;
    for (double x : a) tot += x;
    for (std::size_t k = 0; k < a.size(); ++k)
      terms.push_back({a[k] / tot, g[k], false});
  } else {
    terms.push_back({1.0, 0.0, true});
  }
  std::vector<Dual> tau;
  std::vector<double> dtau, tauv;
  for (const Term& t : terms) {
    if (t.point) {
      tau.push_back(sig);
      dtau.push_back(1.0);
    } else {
      tau.push_back(sig * Dual(t.gamma) /
                    sqrt(Dual(t.gamma * t.gamma) + sig * sig));
      dtau.push_back(std::pow(t.gamma, 3) /
                     std::pow(t.gamma * t.gamma + sg * sg, 1.5));
    }
    tauv.push_back(tau.back().v);
  }

  // block layout
  ChannelMatrices out;
  std::vector<int> off;
  int N = 0;
  for (const Component& cp : ch.comps) {
    const int n = cache.basis(cp.L).size();
    out.L.push_back(cp.L);
    out.block_size.push_back(n);
    off.push_back(N);
    N += n;
  }
  out.H = Mat::Zero(N, N);
  out.S = Mat::Zero(N, N);
  const bool first_order = opt_.spin == SpinMode::FirstOrder;
  const bool in_H = opt_.spin == SpinMode::Smeared;
  const bool smeared_spin =
      in_H || (first_order && opt_.running == Running::Gaussian);
  if (first_order) out.W = Mat::Zero(N, N);
  out.n_dp = nd;
  for (int a = 0; a < nd; ++a) {
    out.dH[a] = Mat::Zero(N, N);
    if (first_order) out.dW[a] = Mat::Zero(N, N);
  }
  auto HB = [&](int i, int j) {
    return Block{&out.H, &out.dH, off[i], off[j], nd};
  };
  auto WB = [&](int i, int j) {
    return Block{&out.W, &out.dW, off[i], off[j], nd};
  };

  // scalar-potential spin-orbit matrix  <V_S'/r>:  b <1/r> (A);  b <e^{-mu
  // r}/r> (B, mu-dependent)
  auto vsp = [&](int Lb, int Lk, const Dual& coef, const Block& blk) {
    if (!opt_.screened) {
      blk.add(coef * b, cache.r_power(Lb, Lk, -1));
    } else {
      const Mat dM = -cache.exp_power(Lb, Lk, mus,
                                      0);  // d/dmu <e^{-mu r}/r> = -<e^{-mu r}>
      blk.add(coef * b, cache.exp_power(Lb, Lk, mus, -1), kMuS, &dM);
    }
  };

  // exact smeared-Coulomb matrices of one coupling term and their derivatives
  // with respect to tau
  auto gauss = [&](int Lb, int Lk, double t, int k) -> const Mat& {
    return cache.gaussian(Lb, Lk, t * t, k);
  };
  auto coulomb_spin = [&](const Block& blk, int Lb, int Lk, const Dual& co_G1,
                          const Dual& co_T, bool want_G1, bool want_T) {
    for (std::size_t k = 0; k < terms.size(); ++k) {
      const double t = tauv[k];
      if (want_G1) {
        const Mat G1 =
            gprime_over_r_matrix(cache.basis(Lb), cache.basis(Lk), t);
        const Mat dG1 = -(4.0 * t * t / kSqrtPi) * gauss(Lb, Lk, t, 0);
        blk.add(co_G1 * Dual(terms[k].w), G1, kSigma, &dG1, dtau[k]);
      }
      if (want_T) {
        const Mat T = tensor_matrix(cache.basis(Lb), cache.basis(Lk), t);
        const Mat dT = -(8.0 * std::pow(t, 4) / kSqrtPi) * gauss(Lb, Lk, t, 2);
        blk.add(co_T * Dual(terms[k].w), T, kSigma, &dT, dtau[k]);
      }
    }
  };

  for (std::size_t i = 0; i < ch.comps.size(); ++i) {
    const int L = ch.comps[i].L, S = ch.comps[i].S, J = ch.J;
    const Block hb = HB(static_cast<int>(i), static_cast<int>(i));
    out.S.block(off[i], off[i], out.block_size[i], out.block_size[i]) =
        cache.overlap(L);
    const Basis& B = cache.basis(L);

    // spin-independent part
    if (opt_.relativistic_kinetic) {
      // sqrt(p^2 + m^2) - m = p^2 / (sqrt(p^2 + m^2) + m) for each quark, and
      // its mass derivative m/E - 1 = -p^2 / (E (E + m)); the rest masses
      // stay in the mass offset
      for (int q : {static_cast<int>(kM1), static_cast<int>(kM2)}) {
        const double m = th(q);
        const Mat K = momentum_diagonal(
            B, [m](double p) { return p * p / (std::sqrt(p * p + m * m) + m); },
            64);
        const Mat dK = momentum_diagonal(
            B,
            [m](double p) {
              const double E = std::sqrt(p * p + m * m);
              return -p * p / (E * (E + m));
            },
            64);
        hb.add(Dual(1.0), K, q, &dK);
      }
    } else {
      hb.add(inv_mr * Dual(0.5), cache.p2(L));
    }
    hb.add(c, cache.overlap(L));
    if (opt_.running == Running::Constant) {
      hb.add(-kappa, cache.r_power(L, L, -1));
    } else {
      for (const Term& t : terms)
        hb.add(-kappa * Dual(t.w), erf_over_r(B, B, t.gamma));
    }
    if (!opt_.screened) {
      hb.add(b, cache.r_power(L, L, 1));
    } else {
      const Dual bm = b / mu;
      const Mat dE =
          -cache.exp_power(L, L, mus, 1);  // d/dmu <e^{-mu r}> = -<r e^{-mu r}>
      hb.add(bm, cache.overlap(L));
      hb.add(-bm, cache.exp_power(L, L, mus, 0), kMuS, &dE);
    }
    // smeared contact term (all spin modes except None): (32 alpha_s/(9 m1 m2))
    // sum_k w_k rho_{tau_k} s1.s2
    if (opt_.spin != SpinMode::None && (opt_.contact_in_all_L || L == 0)) {
      for (std::size_t k = 0; k < terms.size(); ++k) {
        const double t = tauv[k];
        const Dual coef = Dual(32.0 / (9.0 * kSqrtPi) * terms[k].w * s1s2(S)) *
                          als * tau[k] * tau[k] * tau[k] / (m1 * m2);
        const Mat dG = -2.0 * t * gauss(L, L, t, 2);
        hb.add(coef, gauss(L, L, t, 0), kSigma, &dG, dtau[k]);
      }
    }
    // L.S and tensor, diagonal in the component
    const double o_ls = (L > 0 && S == 1) ? ls(L, S, J) : 0.0;
    const double o_t = (L > 0 && S == 1) ? s12_diag(L, J) : 0.0;
    if (first_order && !smeared_spin && L >= 1 && S == 1) {
      // point-like 1/r^3 (finite for L >= 1): a = cV kappa <r^-3> - cS
      // <V_S'/r>;  t = kappa/(4 m1 m2) <r^-3>
      const Block wb = WB(static_cast<int>(i), static_cast<int>(i));
      wb.add(Dual(o_ls) * cV * kappa, cache.r_power(L, L, -3));
      vsp(L, L, -Dual(o_ls) * cS, wb);
      wb.add(Dual(o_t) * kappa / (Dual(4.0) * m1 * m2),
             cache.r_power(L, L, -3));
    }
    if (smeared_spin && L >= 1 && S == 1) {
      const Block target =
          in_H ? hb : WB(static_cast<int>(i), static_cast<int>(i));
      coulomb_spin(target, L, L, -Dual(o_ls) * cV * kappa,
                   -Dual(o_t) * kappa / (Dual(12.0) * m1 * m2), true, true);
      vsp(L, L, -Dual(o_ls) * cS, target);
    }
  }

  // off-diagonal blocks of coupled channels
  if (in_H && ch.comps.size() == 2) {
    const Component &c1 = ch.comps[0], &c2 = ch.comps[1];
    const int L1 = c1.L, L2 = c2.L;
    const Block hb = HB(0, 1);
    if (c1.S == 1 && c2.S == 1 && L1 != L2) {  // tensor: L = J-1 <-> J+1
      coulomb_spin(hb, L1, L2, Dual(0.0),
                   -Dual(s12_offdiag(ch.J)) * kappa / (Dual(12.0) * m1 * m2),
                   false, true);
    } else if (L1 == L2 && c1.S == 0 && c2.S == 1 &&
               L1 >= 1) {  // antisymmetric spin-orbit, m1 != m2
      coulomb_spin(hb, L1, L1, -Dual(ls_minus(L1)) * Bco * kappa, Dual(0.0),
                   true, false);
      vsp(L1, L1, -Dual(ls_minus(L1)) * Bco, hb);
    }
    out.H.block(off[1], off[0], out.block_size[1], out.block_size[0]) =
        out.H.block(off[0], off[1], out.block_size[0], out.block_size[1])
            .transpose();
    for (int a = 0; a < nd; ++a)
      out.dH[a].block(off[1], off[0], out.block_size[1], out.block_size[0]) =
          out.dH[a]
              .block(off[0], off[1], out.block_size[0], out.block_size[1])
              .transpose();
  }

  out.H = 0.5 * (out.H + out.H.transpose());
  for (int a = 0; a < nd; ++a)
    out.dH[a] = 0.5 * (out.dH[a] + out.dH[a].transpose());
  if (first_order) {
    out.W = 0.5 * (out.W + out.W.transpose());
    for (int a = 0; a < nd; ++a)
      out.dW[a] = 0.5 * (out.dW[a] + out.dW[a].transpose());
  }
  return out;
}

}  // namespace qgem
