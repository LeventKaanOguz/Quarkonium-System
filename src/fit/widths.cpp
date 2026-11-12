#include "fit/widths.hpp"

#include <cmath>
#include <limits>
#include <map>

#include "fit/joint.hpp"
#include "models/alpha_s.hpp"
#include "obs/observables.hpp"

namespace qgem {

std::vector<WidthChannel> width_channels(const Dataset& ds) {
  std::vector<WidthChannel> out;
  for (std::size_t i = 0; i < ds.d.size(); ++i) {
    const std::string& lab = ds.d[i].state.channel;
    const int k = static_cast<int>(i);
    if (lab.rfind("3S1", 0) == 0 || lab.rfind("3D1", 0) == 0)
      out.push_back({k, "e+ e-", "ee"});
    else if (lab.rfind("1S0", 0) == 0)
      out.push_back({k, "gamma gamma", "gg_S"});
    else if (lab == "3P0")
      out.push_back({k, "gamma gamma", "gg_P0"});
    else if (lab == "3P2")
      out.push_back({k, "gamma gamma", "gg_P2"});
  }
  return out;
}

std::vector<WidthValue> compute_widths(const Model& model, const Vec& theta,
                                       const Dataset& ds, MatrixCache& cache,
                                       double eQ,
                                       const std::vector<WidthChannel>& ch) {
  const double nan = std::numeric_limits<double>::quiet_NaN();
  std::vector<WidthValue> out(ch.size(), {nan, nan, nan, nan, nan, nan});
  std::map<std::string, SolvedChannel> solved;
  const double mQ = theta(0);
  for (std::size_t k = 0; k < ch.size(); ++k) {
    const Datum& d = ds.d[ch[k].datum];
    try {
      const auto [c, comp] = resolve_channel(d.state.channel, model, theta, 3);
      auto it = solved.find(d.state.channel);
      if (it == solved.end())
        it =
            solved
                .emplace(d.state.channel, solve_channel(model, c, theta, cache))
                .first;
      const SolvedChannel& sc = it->second;
      const int lv = pick_level(sc, comp, d.state.n);
      if (lv >= static_cast<int>(sc.levels.size())) continue;
      const double M = sc.level(lv).M;
      double lo = 0;
      if (ch[k].key == "ee")
        lo = width_ee_SD(origin_derivative(sc, lv, cache, 0),
                         origin_derivative(sc, lv, cache, 2), M, mQ, eQ);
      else if (ch[k].key == "gg_S")
        lo = width_gammagamma_S(R0_squared(sc, lv, cache), M, eQ);
      else {
        const double Rp = origin_derivative(sc, lv, cache, 1);
        lo = width_gammagamma_P(ch[k].key == "gg_P0" ? 0 : 2, Rp * Rp, mQ, eQ);
      }
      const double m1 = theta(0);
      out[k] = {lo,
                lo * qcd_factor(ch[k].key, alpha_s_msbar(M / 2)),
                lo * qcd_factor(ch[k].key, alpha_s_msbar(M / 4)),
                lo * qcd_factor(ch[k].key, alpha_s_msbar(M)),
                p2_expectation(sc, lv, cache) / (m1 * m1),
                M};
    } catch (const std::exception&) {
    }
  }
  return out;
}

double width_theory_error(const WidthValue& w) {
  return std::hypot(0.5 * std::fabs(w.nlo_hi - w.nlo_lo), w.v2 * w.nlo);
}

Mat width_theory_covariance(const std::vector<WidthValue>& w,
                            const std::vector<WidthChannel>& ch) {
  const int n = static_cast<int>(w.size());
  Mat V = Mat::Zero(n, n);
  for (int i = 0; i < n; ++i) {
    const double si = 0.5 * std::fabs(w[i].nlo_hi - w[i].nlo_lo);
    for (int j = 0; j < n; ++j)
      if (ch[i].key == ch[j].key)
        V(i, j) += si * 0.5 * std::fabs(w[j].nlo_hi - w[j].nlo_lo);
    V(i, i) += std::pow(w[i].v2 * w[i].nlo, 2);
  }
  return V;
}

WidthProblem::WidthProblem(const Model& model, const Dataset& ds,
                           ParamSet params, MatrixCache cache, double eQ,
                           WidthData data, const Vec& u_ref)
    : model_(model),
      ds_(ds),
      params_(std::move(params)),
      cache_(std::move(cache)),
      eQ_(eQ),
      data_(std::move(data)) {
  const auto w =
      compute_widths(model_, params_.theta(u_ref), ds_, cache_, eQ_, data_.ch);
  V_ = width_theory_covariance(w, data_.ch);
  V_.diagonal() += data_.s_exp.array().square().matrix();
  if (!V_.allFinite())
    throw std::runtime_error("WidthProblem: a width cannot be computed");
  Linv_ =
      Eigen::LLT<Mat>(V_).matrixL().solve(Mat::Identity(V_.rows(), V_.cols()));
}

Result<Vec> WidthProblem::predictions(const Vec& u) {
  const Vec th = params_.theta(u);
  if (!params_.within_bounds(th)) return fail("outside the parameter bounds");
  const auto w = compute_widths(model_, th, ds_, cache_, eQ_, data_.ch);
  Vec g(static_cast<int>(w.size()));
  for (std::size_t k = 0; k < w.size(); ++k) {
    if (!std::isfinite(w[k].nlo)) return fail("width not computable");
    g(static_cast<int>(k)) = w[k].nlo;
  }
  return g;
}

Result<Vec> WidthProblem::residuals(const Vec& u) {
  const auto g = predictions(u);
  if (!g) return g;
  return Vec(Linv_ * (data_.y - *g));
}

Result<Mat> WidthProblem::jacobian(const Vec& u) {
  // central differences; one-sided where a step leaves the parameter bounds
  // (e.g. mu_s at its lower bound)
  const auto z0 = residuals(u);
  if (!z0) return fail("width Jacobian: " + z0.error().message);
  Mat J(n_data(), n_free());
  for (int a = 0; a < n_free(); ++a) {
    const double h = 1e-5 * std::max(1.0, std::fabs(u(a)));
    Vec up = u, um = u;
    up(a) += h, um(a) -= h;
    const auto zp = residuals(up), zm = residuals(um);
    if (zp && zm)
      J.col(a) = (*zp - *zm) / (2 * h);
    else if (zp)
      J.col(a) = (*zp - *z0) / h;
    else if (zm)
      J.col(a) = (*z0 - *zm) / h;
    else
      return fail("width Jacobian: " + zp.error().message);
  }
  return J;
}

Result<Vec> MassWidthProblem::residuals(const Vec& u) {
  const auto zm = m_.residuals(u);
  if (!zm) return zm;
  const auto zw = w_.residuals(u);
  if (!zw) return zw;
  Vec z(zm->size() + zw->size());
  z << *zm, *zw;
  return z;
}

Result<Mat> MassWidthProblem::jacobian(const Vec& u) {
  const auto Jm = m_.jacobian(u);
  if (!Jm) return Jm;
  const auto Jw = w_.jacobian(u);
  if (!Jw) return Jw;
  Mat J(Jm->rows() + Jw->rows(), n_free());
  J << *Jm, *Jw;
  return J;
}

MassWidthFit fit_masses_and_widths(MassWidthProblem& p, Vec u, double s,
                                   const FitOptions& opt) {
  MassWidthFit out;
  const int Nm = p.masses().n_data();
  for (int it = 0; it < 30; ++it) {
    p.masses().set_sigma_th(s);
    out.fit = levenberg_marquardt(p, u, opt);
    if (!out.fit.J.size()) return out;
    u = out.fit.u;
    out.nu_m = Nm - leverages(out.fit.J).head(Nm).sum();
    const auto r = p.masses().raw_residuals(u);
    if (!r) break;
    const double next = solve_theory_floor(
        *r, p.masses().exp_covariance(), out.nu_m, p.masses().floor_weights());
    const bool done = std::fabs(next - s) < 1e-9;
    s = next;
    if (done) break;
  }
  p.masses().set_sigma_th(s);
  out.fit = levenberg_marquardt(p, u, opt);
  out.sigma_th = s;
  if (out.fit.J.size()) {
    out.nu_m = Nm - leverages(out.fit.J).head(Nm).sum();
    out.chi2_m = out.fit.z.head(Nm).squaredNorm();
    out.chi2_w = out.fit.z.tail(out.fit.z.size() - Nm).squaredNorm();
  }
  return out;
}

}  // namespace qgem
