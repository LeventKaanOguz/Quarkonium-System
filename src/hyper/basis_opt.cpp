#include "hyper/basis_opt.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <map>
#include <string>

#include "fit/chi2_problem.hpp"
#include "math/optimize.hpp"
#include "stats/stats.hpp"

namespace qgem {

double dataset_trace(const Model& model, const Vec& theta, const Dataset& data,
                     const BasisSpec& spec) {
  const double inf = std::numeric_limits<double>::infinity();
  if (!(spec.r1 >= 0.002 && spec.r1 <= 1.0 && spec.rmax >= 3.0 &&
        spec.rmax <= 400.0 && spec.rmax > 1.5 * spec.r1))
    return inf;
  int Lmax = 0;
  for (const Datum& x : data.d)
    for (char c : x.state.channel)
      for (int l = 0; l < 8; ++l)
        if (c == "SPDFGHIK"[l]) Lmax = std::max(Lmax, l);
  try {
    MatrixCache cache(spec);
    std::map<std::string, std::pair<SolvedChannel, int>> solved;  // label -> K
    for (const Datum& x : data.d) {
      const auto [ch, comp] =
          resolve_channel(x.state.channel, model, theta, Lmax);
      auto it = solved.find(ch.label());
      if (it == solved.end())
        it = solved
                 .emplace(
                     ch.label(),
                     std::make_pair(solve_channel(model, ch, theta, cache), 0))
                 .first;
      const SolvedChannel& sc = it->second.first;
      if (sc.levels.empty() || sc.level(0).M < 0) return inf;
      const int k = pick_level(sc, comp, x.state.n);
      if (k >= static_cast<int>(sc.levels.size())) return inf;
      it->second.second = std::max(it->second.second, k + 1);
    }
    double f = 0;
    for (const auto& [label, v] : solved)
      for (int k = 0; k < v.second; ++k) f += v.first.level(k).M;
    return f;
  } catch (const std::exception&) {
    return inf;
  }
}

BasisOptResult optimize_dataset_basis(const Model& model, const Vec& theta,
                                      const Dataset& data, BasisSpec start) {
  int evals = 0;
  auto objective = [&](const Vec& x) {
    ++evals;
    BasisSpec s = start;
    s.r1 = std::exp(x(0));
    s.rmax = std::exp(x(1));
    return dataset_trace(model, theta, data, s);
  };
  Vec x0(2);
  x0 << std::log(start.r1), std::log(start.rmax);
  BasisOptResult r;
  r.start_objective = objective(x0);
  NelderMeadOptions o;
  o.step = 0.3;
  o.xtol = 1e-3;  // 0.1 % in the ranges
  o.ftol = 1e-9;  // GeV
  const NelderMeadResult nm = nelder_mead(objective, x0, o);
  r.spec = start;
  r.spec.r1 = std::exp(nm.x(0));
  r.spec.rmax = std::exp(nm.x(1));
  r.objective = nm.f;
  r.evaluations = evals;
  if (!(r.objective <= r.start_objective)) {  // never worse than the start
    r.spec = start;
    r.objective = r.start_objective;
  }
  return r;
}

Vec dataset_levels(const Model& model, const Vec& theta, const Dataset& data,
                   const BasisSpec& spec) {
  const double nan = std::numeric_limits<double>::quiet_NaN();
  Vec out = Vec::Constant(static_cast<int>(data.d.size()), nan);
  int Lmax = 0;
  for (const Datum& x : data.d)
    for (char c : x.state.channel)
      for (int l = 0; l < 8; ++l)
        if (c == "SPDFGHIK"[l]) Lmax = std::max(Lmax, l);
  MatrixCache cache(spec);
  std::map<std::string, SolvedChannel> solved;
  for (std::size_t i = 0; i < data.d.size(); ++i) {
    const Datum& x = data.d[i];
    try {
      const auto [ch, comp] =
          resolve_channel(x.state.channel, model, theta, Lmax);
      auto it = solved.find(ch.label());
      if (it == solved.end())
        it = solved.emplace(ch.label(), solve_channel(model, ch, theta, cache))
                 .first;
      const SolvedChannel& sc = it->second;
      if (sc.levels.empty() || sc.level(0).M < 0) continue;
      const int k = pick_level(sc, comp, x.state.n);
      if (k < static_cast<int>(sc.levels.size()))
        out(static_cast<int>(i)) = sc.level(k).M;
    } catch (const std::exception&) {
    }
  }
  return out;
}

BasisSizeChoice choose_basis_size(const Model& model, const Vec& theta,
                                  const Dataset& data, BasisSpec spec, int n0,
                                  int steps, double target) {
  BasisSizeChoice c;
  std::vector<Vec> M;
  spec.n = n0;
  for (int k = 0; k <= steps; ++k) {
    c.n.push_back(spec.n);
    M.push_back(dataset_levels(model, theta, data, spec));
    spec.n = 2 * spec.n - 1;
  }
  const int nd = static_cast<int>(data.d.size());
  for (std::size_t k = 0; k + 1 < M.size(); ++k) {
    double worst = 0;
    std::string who;
    for (int i = 0; i < nd; ++i) {
      const double d = M[k](i) - M[k + 1](i);
      if (!std::isfinite(d)) {
        worst = std::numeric_limits<double>::infinity();
        who = data.d[i].label;
        continue;
      }
      c.max_increase = std::max(c.max_increase, -d);
      if (std::fabs(d) > worst) worst = std::fabs(d), who = data.d[i].label;
    }
    c.max_change.push_back(worst);
    c.worst.push_back(who);
    if (!c.chosen && worst < target) c.chosen = c.n[k];
  }
  if (c.chosen) {
    const int k = static_cast<int>(std::find(c.n.begin(), c.n.end(), c.chosen) -
                                   c.n.begin());
    for (int i = 0; i < nd; ++i) {
      const double d = M[k](i) - M.back()(i);
      if (std::isfinite(d))
        c.gap_to_largest = std::max(c.gap_to_largest, std::fabs(d));
    }
  }
  return c;
}

NestedSequence nested_convergence(const Model& model, const Vec& theta,
                                  const Channel& ch, const BasisSpec& start,
                                  const Observable& obs, int steps) {
  NestedSequence s;
  BasisSpec spec = start;
  for (int k = 0; k <= steps; ++k) {
    MatrixCache cache(spec);
    const SolvedChannel sc = solve_channel(model, ch, theta, cache);
    s.n.push_back(spec.n);
    s.value.push_back(obs(sc, cache));
    spec.n = 2 * spec.n - 1;
  }
  const int m = static_cast<int>(s.value.size());
  s.sigma_num = std::fabs(s.value[m - 1] - s.value[m - 2]);
  s.aitken = aitken_limit(s.value[m - 3], s.value[m - 2], s.value[m - 1]);
  s.monotone = true;
  for (int k = 1; k < m; ++k)
    if (s.value[k] > s.value[k - 1] + 1e-12) s.monotone = false;
  return s;
}

}  // namespace qgem
