#include "solve/spectrum.hpp"

#include <cmath>
#include <stdexcept>

namespace qgem {

SolvedChannel solve_channel(const Model& model, const Channel& ch,
                            const Vec& theta, MatrixCache& cache,
                            bool derivatives) {
  SolvedChannel out;
  out.channel = ch;
  out.mats = model.build(ch, theta, cache, derivatives);
  double eta = 1e-12;
  for (int L : out.mats.L) eta = std::max(eta, cache.spec(L).eta);
  if (out.mats.X.size() > 0) {
    // the model assembled H in a fixed orthonormal basis (non-local closure):
    // reuse exactly that basis
    const Mat& X = out.mats.X;
    Mat Hp = X.transpose() * out.mats.H * X;
    Hp = 0.5 * (Hp + Hp.transpose());
    Eigen::SelfAdjointEigenSolver<Mat> eh(Hp);
    out.gep = {eh.eigenvalues(),
               X * eh.eigenvectors(),
               X,
               static_cast<int>(X.cols()),
               0.0,
               1.0};
  } else {
    out.gep = solve_gep(out.mats.H, out.mats.S, eta);
  }
  const double off = model.mass_offset(theta);
  const int nl = static_cast<int>(out.gep.E.size());
  for (int k = 0; k < nl; ++k) {
    Level lv;
    lv.id = {ch.label(), k};
    lv.E = out.gep.E(k);
    lv.c = out.gep.C.col(k);
    if (out.mats.has_W()) lv.shift = lv.c.dot(out.mats.W * lv.c);
    lv.M = off + lv.E + lv.shift;
    int p0 = 0;
    for (int bsz : out.mats.block_size) {
      const Vec cb = lv.c.segment(p0, bsz);
      lv.prob.push_back(cb.dot(out.mats.S.block(p0, p0, bsz, bsz) * cb));
      p0 += bsz;
    }
    out.levels.push_back(std::move(lv));
  }
  return out;
}

Vec SolvedChannel::mass_gradient(int k, const Model& model,
                                 const Vec& theta) const {
  if (mats.n_dp == 0)
    throw std::logic_error(
        "mass_gradient: channel was solved without derivatives");
  const int np = mats.n_dp;
  Vec g = model.mass_offset_gradient(theta);
  const Vec& ck = levels.at(k).c;
  const bool W = mats.has_W();
  for (int a = 0; a < np; ++a) {
    g(a) += ck.dot(mats.dH[a] * ck);
    if (W) {
      g(a) += ck.dot(mats.dW[a] * ck);
      // eigenvector derivative (7.35), S' = 0:  c_k' = sum_{b != k} c_b (c_b^T
      // dH c_k) / (E_k - E_b)
      const Vec dHc = mats.dH[a] * ck;
      const Vec Wc = mats.W * ck;
      double acc = 0.0;
      for (int b = 0; b < gep.C.cols(); ++b) {
        if (b == k) continue;
        const double gap = gep.E(k) - gep.E(b);
        if (std::fabs(gap) < 1e-14) continue;
        const Vec cb = gep.C.col(b);
        acc += cb.dot(dHc) * cb.dot(Wc) / gap;
      }
      g(a) += 2.0 * acc;
    }
  }
  return g;
}

int track_level(const Mat& S, const Vec& c_before, const SolvedChannel& now) {
  int best = 0;
  double bestv = -1.0;
  const Vec Sc = S * c_before;
  for (int b = 0; b < now.gep.C.cols(); ++b) {
    const double o = std::fabs(now.gep.C.col(b).dot(Sc));
    if (o > bestv) {
      bestv = o;
      best = b;
    }
  }
  return best;
}

}  // namespace qgem
