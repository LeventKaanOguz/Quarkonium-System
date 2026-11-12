#include "fit/chi2_problem.hpp"

#include <cmath>
#include <set>
#include <stdexcept>

namespace qgem {
// A datum names a channel ("3P0") or one partial wave of a coupled channel
// ("3S1" of 3S1-3D1). Returns the channel to solve and, for the latter, the
// index of the component that labels the state (-1: the channel itself).
std::pair<Channel, int> resolve_channel(const std::string& label,
                                        const Model& model, const Vec& theta,
                                        int Lmax) {
  const bool equal = std::fabs(theta(0) - theta(1)) < 1e-12 * theta(0);
  const auto chans = model.channels(Lmax, equal);
  for (const Channel& ch : chans)
    if (ch.label() == label) return {ch, -1};
  for (const Channel& ch : chans)
    for (std::size_t i = 0; i < ch.comps.size(); ++i) {
      Channel single;
      single.J = ch.J;
      single.comps = {ch.comps[i]};
      if (ch.comps.size() > 1 && single.label() == label)
        return {ch, static_cast<int>(i)};
    }
  throw std::invalid_argument("unknown channel '" + label + "' for model " +
                              model.name());
}

// The n-th level of the channel that is dominated by component `comp` (all
// levels if comp < 0).
int pick_level(const SolvedChannel& sc, int comp, int n) {
  if (comp < 0) return n;
  int seen = 0;
  for (int k = 0; k < static_cast<int>(sc.levels.size()); ++k)
    if (sc.levels[k].prob[comp] > 0.5 && seen++ == n) return k;
  return static_cast<int>(sc.levels.size());  // not available
}

namespace {
// highest L among the channels named by the data
int needed_Lmax(const Dataset& d) {
  int L = 0;
  const char* letters = "SPDFGHIK";
  for (const Datum& x : d.d)
    for (char c : x.state.channel)
      for (int l = 0; l < 8; ++l)
        if (c == letters[l]) L = std::max(L, l);
  return L;
}
}  // namespace

Chi2Problem::Chi2Problem(const Model& model, Dataset data, ParamSet params,
                         MatrixCache cache, double sigma_th)
    : model_(model),
      data_(std::move(data)),
      params_(std::move(params)),
      cache_(std::move(cache)),
      sigma_th_(sigma_th) {
  n_fit_ = data_.n_fit();
  {
    // parameters that move with the free ones (free + tied): the only
    // derivatives the fit needs
    const Mat K = params_.dtheta_du(params_.u0());
    std::vector<int> active;
    for (int i = 0; i < K.rows(); ++i)
      if (K.row(i).cwiseAbs().maxCoeff() > 0) active.push_back(i);
    model_.restrict_derivatives(active);
  }
  if (n_fit_ == 0) throw std::invalid_argument("Chi2Problem: no data to fit");
  refresh_whitener();
}

void Chi2Problem::refresh_whitener() {
  const Mat V = covariance();
  Eigen::LLT<Mat> llt(V);
  if (llt.info() != Eigen::Success)
    throw std::runtime_error(
        "Chi2Problem: data covariance is not positive definite");
  Linv_ = llt.matrixL().solve(Mat::Identity(n_fit_, n_fit_));
}

Mat Chi2Problem::covariance() const {
  Mat V = data_.V_exp();
  V.diagonal().array() +=
      sigma_th_ * sigma_th_ * floor_weights().array().square();
  return V;
}

Vec Chi2Problem::floor_weights() const {
  Vec w(n_fit_);
  int i = 0;
  for (const Datum* d : data_.fitted()) w(i++) = d->floor_w;
  return w;
}

Mat Chi2Problem::exp_covariance() const { return data_.V_exp(); }

Vec Chi2Problem::exp_errors() const {
  return data_.V_exp().diagonal().cwiseSqrt();
}

void Chi2Problem::set_sigma_th(double s) {
  sigma_th_ = s;
  refresh_whitener();
}

Vec Chi2Problem::measurements() const {
  Vec y(n_fit_);
  int i = 0;
  for (const Datum* d : data_.fitted()) y(i++) = d->y;
  return y;
}

void Chi2Problem::set_measurements(const Vec& y) {
  int i = 0;
  for (Datum& d : data_.d)
    if (d.in_fit) d.y = y(i++);
}

// mu (and its gradient wrt the full physical theta) for the fitted data
Result<Vec> Chi2Problem::mu_fitted(const Vec& theta, Mat* dmu) {
  if (!params_.within_bounds(theta))
    return fail("parameters outside their bounds");
  const int Lmax = needed_Lmax(data_);
  Vec mu(n_fit_);
  if (dmu) *dmu = Mat::Zero(n_fit_, theta.size());
  std::map<std::string, SolvedChannel> solved;
  try {
    int i = 0;
    for (const Datum* d : data_.fitted()) {
      auto it = solved.find(d->state.channel);
      const auto [ch, comp] =
          resolve_channel(d->state.channel, model_, theta, Lmax);
      if (it == solved.end())
        it = solved
                 .emplace(
                     d->state.channel,
                     solve_channel(model_, ch, theta, cache_, dmu != nullptr))
                 .first;
      const SolvedChannel& sc = it->second;
      // Physical domain: a channel whose lowest level falls below zero mass has
      // collapsed (an over-strong point-like contact or spin-orbit term makes
      // H unbounded below; the basis alone then sets the energy). Such points
      // are outside the model, like a parameter bound.
      if (!sc.levels.empty() && sc.level(0).M < 0)
        return fail("collapsed level in channel " + d->state.channel);
      int k = pick_level(sc, comp, d->state.n);
      const int k_order = k;
      // level tracking by maximum overlap with the reference eigenvector (same
      // basis => same S)
      const std::string key =
          d->state.channel + "#" + std::to_string(d->state.n);
      if (refs_set_) {
        auto r = ref_vec_.find(key);
        if (r != ref_vec_.end() && r->second.size() == sc.mats.S.rows())
          k = track_level(sc.mats.S, r->second, sc);
      }
      if (k >= static_cast<int>(sc.levels.size()))
        return fail("level " + key + " not available");
      // In an uncoupled channel the n-th radial state is the n-th level: if
      // tracking disagrees, a spurious level has appeared below it.
      if (ch.comps.size() == 1 && k != k_order)
        return fail("level order changed in " + key + " (spurious level)");
      mu(i) = sc.level(k).M;
      if (dmu) dmu->row(i) = sc.mass_gradient(k, model_, theta).transpose();
      if (!refs_set_) ref_vec_[key] = sc.level(k).c;
      ++i;
    }
    refs_set_ = true;
  } catch (const std::exception& e) {
    return fail(e.what());
  }
  if (!mu.allFinite()) return fail("non-finite prediction");
  return mu;
}

Result<Vec> Chi2Problem::raw_residuals(const Vec& u) {
  auto mu = mu_fitted(params_.theta(u), nullptr);
  if (!mu) return std::unexpected(mu.error());
  return measurements() - *mu;
}

void Chi2Problem::add_prior(int theta_index, double mean, double sigma) {
  if (!(sigma > 0)) throw std::invalid_argument("add_prior: sigma > 0");
  priors_.push_back({theta_index, mean, sigma});
}

Result<Vec> Chi2Problem::residuals(const Vec& u) {
  auto r = raw_residuals(u);
  if (!r) return r;
  Vec z(n_data());
  z.head(n_fit_) = Linv_ * *r;
  const Vec th = params_.theta(u);
  for (std::size_t k = 0; k < priors_.size(); ++k)
    z(n_fit_ + k) = (priors_[k].mean - th(priors_[k].index)) / priors_[k].sigma;
  return z;
}

Result<Mat> Chi2Problem::jacobian(const Vec& u) {
  Mat dmu;
  auto mu = mu_fitted(params_.theta(u), &dmu);
  if (!mu) return std::unexpected(mu.error());
  const Mat K = params_.dtheta_du(u);
  Mat J(n_data(), K.cols());
  J.topRows(n_fit_) = -Linv_ * dmu * K;
  for (std::size_t k = 0; k < priors_.size(); ++k)
    J.row(n_fit_ + k) = -K.row(priors_[k].index) / priors_[k].sigma;
  return J;
}

Result<Vec> Chi2Problem::predictions(const Vec& theta) {
  const int Lmax = needed_Lmax(data_);
  Vec mu(data_.d.size());
  std::map<std::string, SolvedChannel> solved;
  try {
    for (std::size_t i = 0; i < data_.d.size(); ++i) {
      const Datum& d = data_.d[i];
      auto it = solved.find(d.state.channel);
      const auto [ch, comp] =
          resolve_channel(d.state.channel, model_, theta, Lmax);
      if (it == solved.end())
        it = solved
                 .emplace(d.state.channel,
                          solve_channel(model_, ch, theta, cache_, false))
                 .first;
      const int k = pick_level(it->second, comp, d.state.n);
      if (k >= static_cast<int>(it->second.levels.size()))
        return fail("level " + d.state.channel + " not available");
      mu(static_cast<int>(i)) = it->second.level(k).M;
    }
  } catch (const std::exception& e) {
    return fail(e.what());
  }
  return mu;
}

Result<std::pair<Vec, Mat>> Chi2Problem::predictions_with_gradients(
    const Vec& theta) {
  const int Lmax = needed_Lmax(data_);
  const int n = static_cast<int>(data_.d.size());
  Vec mu(n);
  Mat g(n, theta.size());
  std::map<std::string, SolvedChannel> solved;
  try {
    for (int i = 0; i < n; ++i) {
      const Datum& d = data_.d[i];
      auto it = solved.find(d.state.channel);
      const auto [ch, comp] =
          resolve_channel(d.state.channel, model_, theta, Lmax);
      if (it == solved.end())
        it = solved
                 .emplace(d.state.channel,
                          solve_channel(model_, ch, theta, cache_, true))
                 .first;
      const int k = pick_level(it->second, comp, d.state.n);
      if (k >= static_cast<int>(it->second.levels.size()))
        return fail("level " + d.state.channel + " not available");
      mu(i) = it->second.level(k).M;
      g.row(i) = it->second.mass_gradient(k, model_, theta).transpose();
    }
  } catch (const std::exception& e) {
    return fail(e.what());
  }
  return std::make_pair(mu, g);
}

}  // namespace qgem
