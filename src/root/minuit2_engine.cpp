#include "root/minuit2_engine.hpp"

#include <Minuit2/FCNBase.h>
#include <Minuit2/FunctionMinimum.h>
#include <Minuit2/MinosError.h>
#include <Minuit2/MnHesse.h>
#include <Minuit2/MnMigrad.h>
#include <Minuit2/MnMinos.h>
#include <Minuit2/MnPrint.h>
#include <Minuit2/MnUserParameterState.h>
#include <Minuit2/MnUserParameters.h>

#include <cmath>
#include <limits>

namespace qgem {
namespace {

using ROOT::Minuit2::FCNBase;

// chi^2(u) = |z(u)|^2 and its gradient 2 J^T z. A failed evaluation (a
// parameter outside its physical domain) returns a large value, which
// Migrad's line search treats as "step too far".
class Chi2Fcn final : public FCNBase {
 public:
  Chi2Fcn(LeastSquaresProblem& p, bool grad) : p_(&p), grad_(grad) {}
  double operator()(const std::vector<double>& v) const override {
    auto z = p_->residuals(to_vec(v));
    return z ? z->squaredNorm() : 1e30;
  }
  double Up() const override { return 1.0; }
  bool HasGradient() const override { return grad_; }
  std::vector<double> Gradient(const std::vector<double>& v) const override {
    const Vec u = to_vec(v);
    auto z = p_->residuals(u);
    auto J = p_->jacobian(u);
    std::vector<double> g(v.size(), 0.0);
    if (!z || !J) return g;
    const Vec gv = 2.0 * J->transpose() * (*z);
    for (std::size_t i = 0; i < g.size(); ++i) g[i] = gv(static_cast<int>(i));
    return g;
  }

 private:
  static Vec to_vec(const std::vector<double>& v) {
    return Eigen::Map<const Vec>(v.data(), static_cast<int>(v.size()));
  }
  LeastSquaresProblem* p_;
  bool grad_;
};

ROOT::Minuit2::MnUserParameters start(const Vec& u0) {
  ROOT::Minuit2::MnUserParameters up;
  for (int i = 0; i < u0.size(); ++i)
    up.Add("u" + std::to_string(i), u0(i),
           0.05 * std::max(1.0, std::fabs(u0(i))));
  return up;
}

}  // namespace

FitResult minuit2_minimize(LeastSquaresProblem& problem, const Vec& u0,
                           const FitOptions& /*opt*/, const MinuitOptions& mo) {
  ROOT::Minuit2::MnPrint::SetGlobalLevel(0);
  Chi2Fcn fcn(problem, mo.analytic_gradient);
  ROOT::Minuit2::MnMigrad migrad(fcn, start(u0),
                                 ROOT::Minuit2::MnStrategy(mo.strategy));
  const ROOT::Minuit2::FunctionMinimum fm = migrad(mo.max_fcn, mo.tolerance);
  FitResult r;
  r.u = Vec(u0.size());
  for (int i = 0; i < u0.size(); ++i) r.u(i) = fm.UserState().Value(i);
  r.iterations = fm.NFcn();
  r.converged = fm.IsValid();
  r.message = fm.IsValid() ? "Minuit2 Migrad converged"
                           : "Minuit2 Migrad: invalid minimum";
  auto z = problem.residuals(r.u);
  auto J = problem.jacobian(r.u);
  if (!z || !J) {
    r.converged = false;
    r.message = "Minuit2: residuals failed at the minimum";
    return r;
  }
  r.z = *z;
  r.J = *J;
  r.chi2 = r.z.squaredNorm();
  return r;
}

MinuitErrors minuit2_errors(LeastSquaresProblem& problem, const Vec& u_hat,
                            bool minos, const MinuitOptions& mo) {
  ROOT::Minuit2::MnPrint::SetGlobalLevel(0);
  Chi2Fcn fcn(problem, mo.analytic_gradient);
  const ROOT::Minuit2::MnStrategy strat(mo.strategy);
  ROOT::Minuit2::MnMigrad migrad(fcn, start(u_hat), strat);
  ROOT::Minuit2::FunctionMinimum fm = migrad(mo.max_fcn, mo.tolerance);
  ROOT::Minuit2::MnHesse hesse(strat);
  hesse(fcn, fm);
  const int p = static_cast<int>(u_hat.size());
  MinuitErrors e;
  e.valid = fm.IsValid() && fm.HasPosDefCovar();
  e.edm = fm.Edm();
  e.u = Vec(p);
  // a flat direction (e.g. a parameter the model does not use) leaves Hesse
  // without a covariance: report NaN rather than read an empty matrix
  e.cov_u = Mat::Constant(p, p, std::numeric_limits<double>::quiet_NaN());
  const bool have_cov =
      fm.UserState().HasCovariance() &&
      fm.UserState().Covariance().Nrow() == static_cast<unsigned>(p);
  for (int i = 0; i < p; ++i) {
    e.u(i) = fm.UserState().Value(i);
    if (have_cov)
      for (int j = 0; j < p; ++j)
        e.cov_u(i, j) = fm.UserState().Covariance()(i, j);
  }
  e.minos_lo = Vec::Zero(p);
  e.minos_hi = Vec::Zero(p);
  e.minos_valid.assign(p, false);
  if (minos && e.valid) {
    ROOT::Minuit2::MnMinos mn(fcn, fm, strat);
    for (int i = 0; i < p; ++i) {
      const ROOT::Minuit2::MinosError me = mn.Minos(i);
      e.minos_lo(i) = me.Lower();
      e.minos_hi(i) = me.Upper();
      e.minos_valid[i] = me.IsValid();
    }
  }
  e.nfcn = fm.NFcn();
  return e;
}

}  // namespace qgem
