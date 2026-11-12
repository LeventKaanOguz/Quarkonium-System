#include "fit/least_squares.hpp"

#include <algorithm>
#include <cmath>
#include <limits>

namespace qgem {

FitResult levenberg_marquardt(LeastSquaresProblem& prob, Vec u,
                              const FitOptions& opt) {
  FitResult out;
  auto z = prob.residuals(u);
  if (!z) {
    out.u = u;
    out.message = "initial residuals failed: " + z.error().message;
    return out;
  }
  double chi2 = z->squaredNorm(), lambda = opt.lambda0;
  bool stop = false;
  Mat J;
  for (int it = 0; it < opt.max_iter && !stop; ++it) {
    out.iterations = it + 1;
    auto Jr = prob.jacobian(u);
    if (!Jr) {
      out.message = "jacobian failed: " + Jr.error().message;
      break;
    }
    J = *Jr;
    const Mat F = J.transpose() * J;
    const Vec g = -J.transpose() *
                  *z;  // J = dz/du, so the Gauss-Newton gradient is J^T z
    const Vec D = F.diagonal().cwiseMax(1e-300);
    bool accepted = false;
    while (true) {
      Mat A = F;
      A.diagonal() += lambda * D;
      const Vec step = A.ldlt().solve(g);
      auto zn = prob.residuals(u + step);
      const double chi2n =
          zn ? zn->squaredNorm() : std::numeric_limits<double>::infinity();
      if (std::isfinite(chi2n) && chi2n < chi2) {
        const double rel = (chi2 - chi2n) / std::max(chi2, 1e-300);
        double rstep = 0;
        for (int a = 0; a < u.size(); ++a)
          rstep =
              std::max(rstep, std::fabs(step(a)) / (std::fabs(u(a)) + 1e-8));
        u += step;
        *z = *zn;
        chi2 = chi2n;
        lambda = std::max(lambda / 10, 1e-12);
        accepted = true;
        if (rel < opt.ftol && rstep < opt.xtol) {
          out.converged = true;
          stop = true;
        }
        break;
      }
      lambda *= 10;
      if (lambda > 1e12) {
        // no downhill step at any damping: at a minimum to working precision
        out.converged = true;
        out.message = "no further decrease (lambda > 1e12)";
        stop = true;
        break;
      }
    }
    if (!accepted) break;
  }
  if (!out.converged && out.message.empty()) out.message = "max_iter reached";
  auto Jf = prob.jacobian(u);
  out.u = u;
  out.z = *z;
  out.chi2 = chi2;
  out.J = Jf ? *Jf : J;
  return out;
}

double solve_theory_floor(const Vec& Delta, const Mat& V_exp, double target,
                          const Vec& w) {
  const int n = static_cast<int>(Delta.size());
  const Vec w2 = w.size() ? Vec(w.array().square()) : Vec::Ones(n);
  auto chi2 = [&](double x) {
    Mat V = V_exp;
    V.diagonal() += x * w2;
    return Delta.dot(V.llt().solve(Delta));
  };
  if (target <= 0) return std::numeric_limits<double>::infinity();
  if (chi2(0.0) <= target) return 0.0;
  // chi2(x) <= sum Delta_i^2 / (x w_i^2), so this hi brackets the root
  double hi = 2.0 * (Delta.array().square() / w2.array()).sum() / target;
  while (chi2(hi) > target) hi *= 2;
  double lo = 0.0;
  for (int it = 0; it < 200 && hi - lo > 1e-15 * hi; ++it) {
    const double mid = 0.5 * (lo + hi);
    (chi2(mid) > target ? lo : hi) = mid;
  }
  return std::sqrt(0.5 * (lo + hi));
}

Uncertainty analyse(const FitResult& fit, const ParamSet& params) {
  Uncertainty U;
  const Mat F = fit.J.transpose() * fit.J;
  const int p = static_cast<int>(F.rows());
  Eigen::SelfAdjointEigenSolver<Mat> es(F);
  const Vec lam = es.eigenvalues();
  const double tol = lam.maxCoeff() * 1e-14;
  Vec inv = Vec::Zero(p);
  for (int i = 0; i < p; ++i) inv(i) = (lam(i) > tol) ? 1.0 / lam(i) : 0.0;
  U.C = es.eigenvectors() * inv.asDiagonal() * es.eigenvectors().transpose();
  // K = d theta_free / du (square: only the free physical parameters)
  const Mat Kfull = params.dtheta_du(fit.u);
  Mat K(p, p);
  const auto& fi = params.free_indices();
  for (int a = 0; a < p; ++a) K.row(a) = Kfull.row(fi[a]);
  U.Ctheta = K * U.C * K.transpose();
  U.sigma = U.Ctheta.diagonal().cwiseMax(0).cwiseSqrt();
  U.rho = Mat::Zero(p, p);
  for (int a = 0; a < p; ++a)
    for (int b = 0; b < p; ++b)
      U.rho(a, b) = U.Ctheta(a, b) / (U.sigma(a) * U.sigma(b));
  U.sigma_conditional = Vec(p);
  // conditional errors in theta: Fisher in theta = K^-T F K^-1
  const Mat Kinv = K.inverse();
  const Mat Ftheta = Kinv.transpose() * F * Kinv;
  for (int a = 0; a < p; ++a)
    U.sigma_conditional(a) = 1.0 / std::sqrt(Ftheta(a, a));
  // dimensionless Fisher: relative variables d theta / theta
  Vec th(p);
  const Vec thetafull = params.theta(fit.u);
  for (int a = 0; a < p; ++a) th(a) = thetafull(fi[a]);
  const Mat Frel = th.asDiagonal() * Ftheta * th.asDiagonal();
  Eigen::SelfAdjointEigenSolver<Mat> er(Frel);
  U.fisher_eigenvalues = er.eigenvalues();
  U.fisher_eigenvectors = er.eigenvectors();
  U.condition =
      U.fisher_eigenvalues(p - 1) / std::max(U.fisher_eigenvalues(0), 1e-300);
  return U;
}

FitResult minimize(LeastSquaresProblem& problem, const Vec& u0,
                   const FitOptions& opt) {
  return opt.engine ? opt.engine(problem, u0, opt)
                    : levenberg_marquardt(problem, u0, opt);
}

}  // namespace qgem
