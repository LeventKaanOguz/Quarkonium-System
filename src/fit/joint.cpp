#include "fit/joint.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace qgem {

JointProblem::JointProblem(Chi2Problem& a, Chi2Problem& b,
                           const std::vector<std::string>& shared)
    : a_(a), b_(b) {
  const auto na = a.params().free_names(), nb = b.params().free_names();
  for (const auto& s : shared) {
    const auto ia = std::find(na.begin(), na.end(), s);
    const auto ib = std::find(nb.begin(), nb.end(), s);
    if (ia == na.end() || ib == nb.end())
      throw std::invalid_argument("--shared: '" + s +
                                  "' is not free in both systems");
    const auto ta = a.params().specs()[a.params().index(s)].tr;
    const auto tb = b.params().specs()[b.params().index(s)].tr;
    if (ta != tb)
      throw std::invalid_argument("--shared: '" + s +
                                  "' has different transforms");
  }
  int next = a.n_free();
  for (const auto& n : nb) {
    const bool sh = std::find(shared.begin(), shared.end(), n) != shared.end();
    if (sh) {
      b_src_.push_back(
          static_cast<int>(std::find(na.begin(), na.end(), n) - na.begin()));
      ++k_;
    } else {
      b_src_.push_back(next++);
    }
  }
}

Vec JointProblem::u_b(const Vec& u) const {
  Vec ub(b_.n_free());
  for (int j = 0; j < ub.size(); ++j) ub(j) = u(b_src_[j]);
  return ub;
}

Vec JointProblem::join(const Vec& ua, const Vec& ub, bool shared_from_b) const {
  Vec u(n_free());
  u.head(a_.n_free()) = ua;
  for (int j = 0; j < ub.size(); ++j)
    if (b_src_[j] >= a_.n_free() || shared_from_b) u(b_src_[j]) = ub(j);
  return u;
}

std::vector<std::string> JointProblem::names() const {
  auto out = a_.params().free_names();
  const auto nb = b_.params().free_names();
  for (std::size_t j = 0; j < nb.size(); ++j)
    if (b_src_[j] >= a_.n_free()) out.push_back(nb[j] + "(b)");
  return out;
}

Vec JointProblem::dtheta_du(const Vec& u) const {
  Vec d(n_free());
  const Mat Ka = a_.params().dtheta_du(u_a(u));
  const Mat Kb = b_.params().dtheta_du(u_b(u));
  for (int j = 0; j < a_.n_free(); ++j)
    d(j) = Ka(a_.params().free_indices()[j], j);
  for (int j = 0; j < b_.n_free(); ++j)
    if (b_src_[j] >= a_.n_free())
      d(b_src_[j]) = Kb(b_.params().free_indices()[j], j);
  return d;
}

Result<Vec> JointProblem::residuals(const Vec& u) {
  auto za = a_.residuals(u_a(u));
  if (!za) return za;
  auto zb = b_.residuals(u_b(u));
  if (!zb) return zb;
  Vec z(za->size() + zb->size());
  z << *za, *zb;
  return z;
}

Result<Mat> JointProblem::jacobian(const Vec& u) {
  auto Ja = a_.jacobian(u_a(u));
  if (!Ja) return Ja;
  auto Jb = b_.jacobian(u_b(u));
  if (!Jb) return Jb;
  Mat J = Mat::Zero(Ja->rows() + Jb->rows(), n_free());
  J.topLeftCorner(Ja->rows(), Ja->cols()) = *Ja;
  for (int j = 0; j < Jb->cols(); ++j)
    J.block(Ja->rows(), b_src_[j], Jb->rows(), 1) += Jb->col(j);
  return J;
}

Vec leverages(const Mat& J) {
  Eigen::SelfAdjointEigenSolver<Mat> es(J.transpose() * J);
  const Vec lam = es.eigenvalues();
  const double tol = lam.maxCoeff() * 1e-14;
  Vec inv = Vec::Zero(lam.size());
  for (int i = 0; i < lam.size(); ++i)
    if (lam(i) > tol) inv(i) = 1 / lam(i);
  const Mat W = J * es.eigenvectors();  // N x p
  return (W.array().square().rowwise() * inv.transpose().array())
      .rowwise()
      .sum();
}

JointFit fit_joint(JointProblem& jp, Vec u, double sa, double sb,
                   const FitOptions& opt) {
  JointFit out;
  const int Na = jp.a().n_data();
  for (int it = 0; it < 30; ++it) {
    out.iterations = it + 1;
    jp.a().set_sigma_th(sa);
    jp.b().set_sigma_th(sb);
    out.fit = levenberg_marquardt(jp, u, opt);
    if (!out.fit.J.size()) return out;
    u = out.fit.u;
    const Vec h = leverages(out.fit.J);
    out.nu_a = Na - h.head(Na).sum();
    out.nu_b = (h.size() - Na) - h.tail(h.size() - Na).sum();
    const auto ra = jp.a().raw_residuals(jp.u_a(u));
    const auto rb = jp.b().raw_residuals(jp.u_b(u));
    if (!ra || !rb) break;
    const double na = solve_theory_floor(*ra, jp.a().exp_covariance(), out.nu_a,
                                         jp.a().floor_weights());
    const double nb = solve_theory_floor(*rb, jp.b().exp_covariance(), out.nu_b,
                                         jp.b().floor_weights());
    const bool done = std::fabs(na - sa) < 1e-9 && std::fabs(nb - sb) < 1e-9;
    sa = na, sb = nb;
    if (done) break;
  }
  jp.a().set_sigma_th(sa);
  jp.b().set_sigma_th(sb);
  out.fit = levenberg_marquardt(jp, u, opt);
  out.sigma_a = sa, out.sigma_b = sb;
  if (out.fit.J.size()) {
    const Vec h = leverages(out.fit.J);
    out.nu_a = Na - h.head(Na).sum();
    out.nu_b = (h.size() - Na) - h.tail(h.size() - Na).sum();
  }
  return out;
}

}  // namespace qgem
