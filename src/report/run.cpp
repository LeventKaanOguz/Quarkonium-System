#include "report/run.hpp"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <limits>
#include <map>
#include <memory>
#include <numbers>
#include <optional>
#include <ostream>
#include <sstream>

#include "core/constants.hpp"
#include "data/dataset.hpp"
#include "fit/driver.hpp"
#include "fit/joint.hpp"
#include "fit/widths.hpp"
#include "hyper/basis_opt.hpp"
#include "models/alpha_s.hpp"
#include "models/gi_model.hpp"
#include "models/nr_model.hpp"
#include "obs/observables.hpp"
#include "report/json.hpp"
#include "stats/rng.hpp"
#include "stats/stats.hpp"

namespace qgem {
namespace {

struct System {
  std::string name, label, states;
  double m_start, alpha_start, b_start, sigma_start, sigma0_start, charge;
};

const System& system_by_name(const std::string& n) {
  static const std::vector<System> all = {
      {"ccbar", "charmonium", "config/ccbar_fit.states", 1.48, 0.55, 0.14, 1.1,
       1.8, 2.0 / 3.0},
      {"bbar", "bottomonium", "config/bbar_fit.states", 4.75, 0.38, 0.17, 1.6,
       1.8, -1.0 / 3.0},
      {"bc", "B_c", "config/bc.states", 0, 0, 0, 0, 0, 0}};
  for (const System& s : all)
    if (s.name == n) return s;
  throw std::runtime_error("unknown system '" + n + "' (ccbar, bbar, bc)");
}

// p-value for print: fixed below 1e-3 would read as zero
std::string fmt_p(double p) {
  if (!std::isfinite(p)) return "-";
  if (p < 1e-300) return "< 1e-300";
  std::ostringstream s;
  if (p >= 1e-3)
    s << std::fixed << std::setprecision(3) << p;
  else
    s << std::scientific << std::setprecision(1) << p;
  return s.str();
}

struct Built {
  std::unique_ptr<Model> model;
  ParamSet params;
  bool gi = false;
};

// Remove the --fix parameters from the free list and set their values; add
// the --free parameters.
ParamSet apply_fixed(const RunOptions& o, const ParamSet& ps) {
  if (o.fixed.empty() && o.freed.empty()) return ps;
  std::vector<std::string> free;
  for (const auto& n : ps.free_names()) {
    bool f = false;
    for (const auto& [name, v] : o.fixed) f = f || name == n;
    if (!f) free.push_back(n);
  }
  for (const auto& n : o.freed) {
    ps.index(n);  // throws if the model has no such parameter
    if (std::find(free.begin(), free.end(), n) == free.end()) free.push_back(n);
  }
  Vec t = ps.theta0();
  for (const auto& [name, v] : o.fixed) t(ps.index(name)) = v;
  return ps.with_values(t).with_free(free);
}

Built build_model(const RunOptions& o, const System& sys) {
  Built b;
  const bool scr = (o.model == "B" || o.model == "D");
  if (o.model == "C" || o.model == "D") {
    b.gi = true;
    auto m = std::make_unique<GIModel>(
        GIModelOptions{.screened = scr, .nq = o.gi_nq});
    const GiRunning r = o.running == "const"  ? GiRunning::Constant
                        : o.running == "free" ? GiRunning::Free
                                              : GiRunning::Fixed;
    const double mq = sys.name == "ccbar" ? 1.628 : 4.977;
    b.params = m->godfrey_isgur_params(mq, mq, r, o.with_c);
    std::vector<std::string> free = {"m1", "b", "sigma0"};
    if (o.with_c) free.push_back("c");
    if (r == GiRunning::Free) {
      free.push_back("alpha1");
      free.push_back("alpha2");
      free.push_back("alpha3");
    }
    if (r == GiRunning::Constant) free.push_back("alpha1");
    if (scr && o.fit_mu_s) free.push_back("mu_s");
    b.params = b.params.with_free(free);
    {
      Vec t = b.params.theta0();
      t(b.params.index("mu_s")) = o.mu_s;
      b.params = b.params.with_values(t);
    }
    b.model = std::move(m);
  } else if (o.model == "A" || o.model == "B") {
    if (o.running != "const" && o.running != "gi" && o.running != "frozen" &&
        o.running != "apt")
      throw std::runtime_error(
          "models A/B: --running const|gi|frozen|apt (free is C/D only)");
    NRModelOptions no;
    no.screened = scr;
    no.running = o.running == "const" ? Running::Constant : Running::Gaussian;
    if (o.running == "frozen" || o.running == "apt") {
      const auto shape =
          o.running == "frozen"
              ? make_frozen_coupling(o.lambda_qcd, o.gluon_mass, o.nf)
              : make_apt_coupling(o.lambda_qcd, o.nf);
      const GaussianShape g = fit_gaussian_shape(*shape, default_shape_gamma());
      no.shape_w = g.w;
      no.shape_gamma = g.gamma;
      no.shape_label = o.running == "frozen" ? "frozen" : "APT";
    }
    if (o.kinetic != "nr" && o.kinetic != "rel")
      throw std::runtime_error("--kinetic nr|rel");
    no.relativistic_kinetic = o.kinetic == "rel";
    no.spin = o.spin == "none"      ? SpinMode::None
              : o.spin == "smeared" ? SpinMode::Smeared
                                    : SpinMode::FirstOrder;
    auto m = std::make_unique<NRModel>(no);
    b.params = m->default_params(sys.m_start, sys.m_start, sys.alpha_start,
                                 sys.b_start, 0.0, sys.sigma_start);
    std::vector<ParamSpec> sp = b.params.specs();
    sp[1].tied_to = "m1";
    b.params = ParamSet(sp);
    std::vector<std::string> free = {"m1", "alpha_s", "b", "sigma"};
    if (o.with_c) free.push_back("c");
    if (scr && o.fit_mu_s) free.push_back("mu_s");
    b.params = b.params.with_free(free);
    {
      Vec t = b.params.theta0();
      t(b.params.index("mu_s")) = o.mu_s;
      b.params = b.params.with_values(t);
    }
    b.model = std::move(m);
  } else {
    throw std::runtime_error("--model must be A, B, C or D");
  }
  if (b.gi && o.running != "const" && o.running != "gi" && o.running != "free")
    throw std::runtime_error("models C/D: --running const|gi|free");
  b.params = apply_fixed(o, b.params);
  return b;
}

// Solve the channel of a datum at theta and return it with the level index.
struct Solved {
  SolvedChannel sc;
  int level;
  int comp;
};
Solved solve_state(const Model& model, const Vec& theta, MatrixCache& cache,
                   const Datum& d) {
  auto [ch, comp] = resolve_channel(d.state.channel, model, theta, 3);
  Solved s{solve_channel(model, ch, theta, cache), 0, comp};
  s.level = pick_level(s.sc, comp, d.state.n);
  if (s.level >= static_cast<int>(s.sc.levels.size()))
    throw std::runtime_error("level " + d.state.channel + " not available");
  return s;
}

std::string clean_label(std::string s) {  // PDG writes "B_c()+"
  for (auto p = s.find("()"); p != std::string::npos; p = s.find("()"))
    s.erase(p, 2);
  return s;
}

double percentile(std::vector<double> v, double q) {
  if (v.empty()) return 0;
  std::sort(v.begin(), v.end());
  const double x = q * (static_cast<double>(v.size()) - 1);
  const auto i = static_cast<std::size_t>(x);
  const double f = x - static_cast<double>(i);
  return i + 1 < v.size() ? (1 - f) * v[i] + f * v[i + 1] : v[i];
}

// Summary of samples of the free parameters given in u variables: mean, std,
// 16% and 84% quantiles of the physical parameters.
std::string summarize_samples(const Mat& usamp, const ParamSet& ps,
                              const std::string& title, std::ostream& out) {
  const auto names = ps.free_names();
  const int p = static_cast<int>(names.size());
  const int n = static_cast<int>(usamp.rows());
  json::Object j;
  Vec mean(p), sd(p), q16(p), q84(p);
  for (int k = 0; k < p; ++k) {
    std::vector<double> v(n);
    for (int r = 0; r < n; ++r)
      v[r] = ps.theta(Vec(usamp.row(r).transpose()))(ps.free_indices()[k]);
    double m = 0;
    for (double x : v) m += x;
    m /= std::max(1, n);
    double var = 0;
    for (double x : v) var += (x - m) * (x - m);
    mean(k) = m;
    sd(k) = std::sqrt(var / std::max(1, n - 1));
    q16(k) = percentile(v, 0.1587);
    q84(k) = percentile(v, 0.8413);
  }
  out << "\n" << title << " (" << n << " samples):\n";
  for (int k = 0; k < p; ++k)
    out << "  " << std::left << std::setw(9) << names[k] << std::right
        << std::fixed << std::setprecision(5) << std::setw(11) << mean(k)
        << " +- " << std::setw(9) << sd(k) << "   68% [" << q16(k) << ", "
        << q84(k) << "]\n"
        << std::defaultfloat;
  j.integer("n", n)
      .put("mean", json::array(mean))
      .put("std", json::array(sd))
      .put("q16", json::array(q16))
      .put("q84", json::array(q84));
  return j.render(6);
}

std::string state_name(const Datum& d) {
  return d.state.channel +
         (d.state.n ? "(" + std::to_string(d.state.n + 1) + ")" : "");
}

std::string fmt(double v, int w, int p) {
  std::ostringstream s;
  s << std::fixed << std::setw(w) << std::setprecision(p) << v;
  return s.str();
}

// B_c parameters from the two fits: log-interpolation in the reduced mass (as
// for the running coupling), linear if a parameter is not positive; masses are
// the fitted quark masses.
Vec interpolate_bc(const Vec& cc, const Vec& bb) {
  Vec th = cc;
  const double mc = cc(0), mb = bb(0);
  th(0) = mc;
  th(1) = mb;
  const double rc = mc / 2, rb = mb / 2, rbc = mc * mb / (mc + mb);
  const double t =
      (std::log(rbc) - std::log(rc)) / (std::log(rb) - std::log(rc));
  for (int i = 2; i < th.size(); ++i) {
    if (cc(i) == bb(i))
      th(i) = cc(i);
    else if (cc(i) > 0 && bb(i) > 0)
      th(i) = std::exp((1 - t) * std::log(cc(i)) + t * std::log(bb(i)));
    else
      th(i) = (1 - t) * cc(i) + t * bb(i);
  }
  return th;
}

}  // namespace

int run(const RunOptions& opt, std::ostream& out) {
  const std::string pdg_path = opt.pdg_file.empty()
                                   ? opt.root + "/data/pdg-2026.0.sqlite"
                                   : opt.pdg_file;
  PdgDatabase pdg(pdg_path);
  std::map<std::string, std::pair<Vec, Built*>>
      fitted;  // theta of each fitted system
  // the separate fits, kept for the joint fit (--shared)
  struct Kept {
    std::unique_ptr<Chi2Problem> prob;
    Vec u;
    double sigma_th = 0, chi2 = 0;
  };
  std::map<std::string, Kept> kept;
  bool joint_done = false;
  std::vector<std::unique_ptr<Built>> keep;
  json::Object json_all;
  std::string json_systems = "[";
  int nsys = 0;

  const std::map<std::string, std::string> coupling = {
      {"const", "constant"},
      {"gi", "running (Godfrey-Isgur three-Gaussian form)"},
      {"free", "running (Godfrey-Isgur form, free amplitudes)"},
      {"frozen", "running (frozen, gluon mass)"},
      {"apt", "running (analytic perturbation theory)"}};
  const auto cit = coupling.find(opt.running);
  out << "# qgem: model " << opt.model << ", alpha_s "
      << (cit != coupling.end() ? cit->second : opt.running)
      << ", additive constant c " << (opt.with_c ? "FREE" : "fixed to 0")
      << "\n# basis: n = " << opt.basis.n << " Gaussians, ranges "
      << opt.basis.r1 << " ... " << opt.basis.rmax
      << " GeV^-1 (geometric, ratio "
      << fmt(std::pow(opt.basis.rmax / opt.basis.r1,
                      1.0 / std::max(1, opt.basis.n - 1)),
             0, 3)
      << ")\n# PDG file: " << pdg_path << "\n\n";

  for (const std::string& sname : opt.systems) {
    BasisSpec basis = opt.basis;  // per system: --optimize-basis may move it
    const System& sys = system_by_name(sname);
    const bool is_bc = sname == "bc";
    Dataset ds = load_dataset(opt.root + "/" + sys.states, pdg);
    if (const auto fw = opt.floor_weights.find(sname);
        fw != opt.floor_weights.end()) {
      for (const auto& [label, w] : fw->second) {
        bool found = false;
        for (Datum& x : ds.d)
          if (clean_label(x.label) == label) x.floor_w = w, found = true;
        if (!found)
          throw std::runtime_error("--floor-weights: no state '" + label +
                                   "' in " + sname);
      }
    }
    out << "================ " << sys.label << " (" << ds.d.size()
        << " PDG states) ================\n";

    Vec theta;
    Built* built = nullptr;
    double sigma_th = 0, chi2 = 0;
    Uncertainty U;
    LooResult loo;
    Mat Cfull;
    int nu = 0;
    ParamSet ps;
    std::unique_ptr<Chi2Problem> prob;
    std::string ptable;
    json::Object jan;  // analyses of this system (JSON)
    FitResult best_fit;
    Mat post;  // posterior draws (u, ln sigma_th) of --mcmc

    if (!is_bc) {
      keep.push_back(std::make_unique<Built>(build_model(opt, sys)));
      built = keep.back().get();
      ps = built->params;
      // A fit from the default start, then `multistart` random restarts
      // (each with its own problem, so level tracking starts afresh). Score:
      // the theory floor (iterated) or chi^2 (fixed floor); lower is better.
      const bool iter_floor = opt.fixed_floor < 0;
      auto score = [&](const FloorFit& f) {
        if (f.fit.J.size() == 0) return std::numeric_limits<double>::infinity();
        return iter_floor ? f.sigma_th : f.fit.chi2;
      };
      // Nested start: when the model extends a smaller one (free c, free
      // screening mass, free GI amplitudes), fit the smaller model first and
      // release the extra parameters from its minimum. LM only goes downhill,
      // so chi^2(extended) <= chi^2(nested) and the Wilks tests are valid.
      {
        std::vector<std::string> ext = {"c", "mu_s"};
        if (built->gi && opt.running == "free")
          for (const char* a : {"alpha1", "alpha2", "alpha3"}) ext.push_back(a);
        std::vector<std::string> base;
        for (const auto& n : ps.free_names())
          if (std::find(ext.begin(), ext.end(), n) == ext.end())
            base.push_back(n);
        if (base.size() < ps.free_names().size()) {
          const ParamSet ps0 = ps.with_free(base);
          Chi2Problem p0(*built->model, ds, ps0, MatrixCache(basis), 0.010);
          const FloorFit f0 =
              fit_with_floor(p0, ps0.u0(), 0.010, opt.fit, opt.fixed_floor);
          if (f0.fit.J.size()) {
            ps = ps.with_values(ps0.theta(f0.fit.u));
            out << "nested start: the model without the extra parameters "
                   "gives "
                << (iter_floor
                        ? "sigma_th = " + fmt(1e3 * f0.sigma_th, 0, 3) + " MeV"
                        : "chi2 = " + fmt(f0.fit.chi2, 0, 4))
                << "; released from there\n";
            jan.number("nested_start_score",
                       iter_floor ? f0.sigma_th : f0.fit.chi2);
          }
        }
      }
      prob = std::make_unique<Chi2Problem>(*built->model, ds, ps,
                                           MatrixCache(basis), 0.010);
      FloorFit ff =
          fit_with_floor(*prob, ps.u0(), 0.010, opt.fit, opt.fixed_floor);
      std::vector<double> start_scores = {score(ff)};
      int best_start = 0;
      // One more start from theta; keeps it if it beats the best so far.
      auto try_start = [&](const Vec& th) {
        if (!ps.within_bounds(th)) {
          start_scores.push_back(std::numeric_limits<double>::infinity());
          return;
        }
        auto trial = std::make_unique<Chi2Problem>(*built->model, ds, ps,
                                                   MatrixCache(basis), 0.010);
        FloorFit f;
        try {
          f = fit_with_floor(*trial, ps.u(th), 0.010, opt.fit, opt.fixed_floor);
        } catch (const std::exception&) {
          start_scores.push_back(std::numeric_limits<double>::infinity());
          return;
        }
        start_scores.push_back(score(f));
        if (score(f) < score(ff) - 1e-9) {
          ff = std::move(f);
          prob = std::move(trial);
          best_start = static_cast<int>(start_scores.size()) - 1;
        }
      };
      // Seeds: known minima (another minimiser, the nested model's best fit).
      int n_seeds = 0;
      if (const auto it = opt.starts.find(sname); it != opt.starts.end())
        for (const auto& seed : it->second) {
          Vec th = ps.theta0();
          for (int i : ps.free_indices())
            for (const auto& [name, v] : seed)
              if (ps.specs()[i].name == name) th(i) = v;
          try_start(th);
          ++n_seeds;
        }
      // Random restarts around the best point so far.
      const Vec centre = ff.fit.J.size() ? ps.theta(ff.fit.u) : ps.theta0();
      Rng rng(opt.seed, std::hash<std::string>{}(sname));
      for (int k = 1; k <= opt.multistart; ++k) {
        Vec th = centre;
        for (int i : ps.free_indices())
          th(i) = th(i) != 0 ? th(i) * std::exp(0.3 * rng.normal())
                             : 0.05 * rng.normal();
        try_start(th);
      }
      // Variational basis: minimise the level sum of the data's channels over
      // (r1, rmax) at the fitted theta, refit, repeat until the ranges settle.
      if (opt.optimize_basis && ff.fit.J.size()) {
        for (int it = 0; it < 4; ++it) {
          const BasisOptResult bo = optimize_dataset_basis(
              *built->model, ps.theta(ff.fit.u), ds, basis);
          out << "basis optimisation " << it + 1
              << ": r1 = " << fmt(bo.spec.r1, 0, 4)
              << ", rmax = " << fmt(bo.spec.rmax, 0, 2) << " GeV^-1; level sum "
              << fmt(bo.start_objective, 0, 7) << " -> "
              << fmt(bo.objective, 0, 7) << " GeV (" << bo.evaluations
              << " evaluations)\n";
          const bool moved =
              std::fabs(std::log(bo.spec.r1 / basis.r1)) > 0.01 ||
              std::fabs(std::log(bo.spec.rmax / basis.rmax)) > 0.01;
          if (!moved) break;
          auto refit = std::make_unique<Chi2Problem>(
              *built->model, ds, ps, MatrixCache(bo.spec), 0.010);
          FloorFit f = fit_with_floor(*refit, ff.fit.u, ff.sigma_th, opt.fit,
                                      opt.fixed_floor);
          if (!f.fit.J.size()) break;  // keep the previous basis and fit
          basis = bo.spec;
          ff = std::move(f);
          prob = std::move(refit);
        }
      }
      prob->set_sigma_th(ff.sigma_th);
      best_fit = ff.fit;
      if (opt.multistart > 0 || n_seeds > 0) {
        out << "multistart: " << start_scores.size() << " starts (default, "
            << n_seeds << (n_seeds == 1 ? " seed, " : " seeds, ")
            << opt.multistart << " random), best = start " << best_start
            << "; scores (" << (iter_floor ? "sigma_th MeV" : "chi2") << "):";
        for (double x : start_scores)
          out << " " << fmt(iter_floor ? 1e3 * x : x, 0, 2);
        out << "\n";
        Vec sc(static_cast<int>(start_scores.size()));
        for (int k = 0; k < sc.size(); ++k) sc(k) = start_scores[k];
        jan.put("multistart_scores", json::array(sc))
            .integer("multistart_best", best_start)
            .integer("multistart_seeds", n_seeds);
      }
      if (ff.fit.J.size() == 0) {
        out << "fit failed: " << ff.fit.message << "\n\n";
        continue;
      }
      if (!ff.fit.converged)
        out << "warning: fit did not converge (" << ff.fit.message << ")\n";
      theta = ps.theta(ff.fit.u);
      sigma_th = ff.sigma_th;
      chi2 = ff.fit.chi2;
      nu = prob->n_data() - prob->n_free();
      U = analyse(ff.fit, ps);
      loo = leave_one_out(ff.fit);
      const Mat K = ps.dtheta_du(ff.fit.u);
      Cfull = K * U.C * K.transpose();
      fitted[sname] = {theta, built};
      // Goodness of fit. With the floor iterated, chi2/nu = 1 by
      // construction and sigma_th itself measures the misfit; with a fixed
      // floor chi2 is a test statistic (conditional on that floor). The fit
      // with experimental errors alone tests the model at the precision of
      // the data and motivates sigma_th.
      double chi2_exp = std::numeric_limits<double>::quiet_NaN();
      prob->set_sigma_th(0.0);
      if (const auto z = prob->residuals(ff.fit.u))
        chi2_exp = z->head(prob->n_measurements()).squaredNorm();
      prob->set_sigma_th(ff.sigma_th);  // the fit's weights, for the refits
      // The chi^2 test inverted: "the model is accurate to sigma" is accepted
      // for every sigma at which chi2(sigma) lies in the central chi^2_nu
      // interval. Parameters are held at the fit (exact to first order).
      const double nan = std::numeric_limits<double>::quiet_NaN();
      double ci68[2] = {nan, nan}, ci95[2] = {nan, nan};
      if (iter_floor && nu > 0)
        if (const auto raw = prob->raw_residuals(ff.fit.u)) {
          const Mat Ve = prob->exp_covariance();
          const Vec w = prob->floor_weights();
          auto at = [&](double q) {
            return solve_theory_floor(*raw, Ve, chi2_quantile(q, nu), w);
          };
          ci68[0] = at(0.5 + 0.6827 / 2), ci68[1] = at(0.5 - 0.6827 / 2);
          ci95[0] = at(0.5 + 0.9545 / 2), ci95[1] = at(0.5 - 0.9545 / 2);
        }
      out << "model " << built->model->name() << "\n  " << prob->n_data()
          << " fitted states, " << prob->n_free()
          << " free parameters, nu = " << nu << "\n";
      if (iter_floor)
        out << "  theory floor sigma_th = " << fmt(1e3 * sigma_th, 0, 3)
            << " MeV (chi2/nu = 1 by construction)\n"
            << "  chi2 test of 'the model is accurate to sigma_th' accepts "
            << "sigma_th in [" << fmt(1e3 * ci68[0], 0, 2) << ", "
            << fmt(1e3 * ci68[1], 0, 2) << "] MeV at 68% CL, ["
            << fmt(1e3 * ci95[0], 0, 2) << ", " << fmt(1e3 * ci95[1], 0, 2)
            << "] MeV at 95% CL\n";
      else
        out << "  fixed theory floor sigma_th = " << fmt(1e3 * sigma_th, 0, 3)
            << " MeV: chi2 = " << fmt(chi2, 0, 3) << " for nu = " << nu
            << ", P(chi2_nu >= chi2) = " << fmt_p(chi2_sf(chi2, nu))
            << " (conditional on the floor)\n";
      out << "  experimental errors only (sigma_th = 0): chi2 = "
          << fmt(chi2_exp, 0, 1) << " for nu = " << nu << ", P "
          << (chi2_sf(chi2_exp, nu) < 1e-300 ? "" : "= ")
          << fmt_p(chi2_sf(chi2_exp, nu))
          << " (the model does not reach the precision of the data: "
             "hence the floor)\n\n";
      jan.put("sigma_th_ci68_GeV",
              json::array(Vec((Vec(2) << ci68[0], ci68[1]).finished())))
          .put("sigma_th_ci95_GeV",
               json::array(Vec((Vec(2) << ci95[0], ci95[1]).finished())));
      jan.number("chi2_exp_only", chi2_exp)
          .number("gof_p", iter_floor ? std::numeric_limits<double>::quiet_NaN()
                                      : chi2_sf(chi2, nu));
      const auto names = ps.free_names();
      out << "free parameters:\n";
      for (int k = 0; k < static_cast<int>(names.size()); ++k)
        out << "  " << std::left << std::setw(9) << names[k] << std::right
            << fmt(theta(ps.free_indices()[k]), 11, 5) << " +- "
            << fmt(U.sigma(k), 9, 5) << "\n";
      out << "fixed parameters:";
      for (int i = 0; i < ps.n_total(); ++i)
        if (!ps.specs()[i].free && ps.specs()[i].tied_to.empty())
          out << " " << ps.specs()[i].name << "=" << theta(i);
      out << "\n";
      const int p = static_cast<int>(names.size());
      if (p > 1) {
        out << "\ncorrelations:\n          ";
        for (const auto& n : names) out << std::setw(9) << n;
        out << "\n";
        for (int i = 0; i < p; ++i) {
          out << "  " << std::left << std::setw(8) << names[i] << std::right;
          for (int j = 0; j < p; ++j) out << fmt(U.rho(i, j), 9, 3);
          out << "\n";
        }
        out << "condition number of the relative Fisher matrix: "
            << std::scientific << std::setprecision(2) << U.condition
            << std::defaultfloat
            << (U.condition > 1e10 ? "  (DEGENERATE: fix a parameter, add a "
                                     "prior or more data)"
                                   : "")
            << "\n";
      }
      if (opt.nested_target > 0) {
        const BasisSizeChoice bc = choose_basis_size(
            *built->model, theta, ds, basis, 10, 4, opt.nested_target);
        out << "basis convergence on nested sizes (r1 = " << basis.r1
            << ", rmax = " << basis.rmax
            << " GeV^-1; largest change of any state to the next size):\n";
        for (std::size_t k = 0; k < bc.max_change.size(); ++k)
          out << "  n = " << std::setw(3) << bc.n[k] << " -> " << std::setw(3)
              << bc.n[k + 1] << ": " << fmt(1e3 * bc.max_change[k], 10, 4)
              << " MeV  (" << clean_label(bc.worst[k]) << ")\n";
        out << "  smallest n below " << fmt(1e3 * opt.nested_target, 0, 3)
            << " MeV: "
            << (bc.chosen ? std::to_string(bc.chosen) : std::string("none"))
            << "; its largest distance to n = " << bc.n.back() << ": "
            << fmt(1e3 * bc.gap_to_largest, 0, 4)
            << " MeV; largest rise of a level along the sequence: "
            << fmt(1e3 * bc.max_increase, 0, 4) << " MeV\n";
        std::vector<std::string> worst;
        for (const auto& w : bc.worst) worst.push_back(clean_label(w));
        Vec ns(static_cast<int>(bc.n.size())),
            mc(static_cast<int>(bc.max_change.size()));
        for (std::size_t k = 0; k < bc.n.size(); ++k) ns(k) = bc.n[k];
        for (std::size_t k = 0; k < bc.max_change.size(); ++k)
          mc(k) = bc.max_change[k];
        json::Object nj;
        nj.put("n", json::array(ns))
            .put("max_change_GeV", json::array(mc))
            .put("worst", json::array(worst))
            .integer("chosen", bc.chosen)
            .number("target_GeV", opt.nested_target)
            .number("gap_to_largest_GeV", bc.gap_to_largest)
            .number("max_increase_GeV", bc.max_increase);
        jan.put("nested", nj.render(6));
      }
      jan.number("basis_r1", basis.r1)
          .number("basis_rmax", basis.rmax)
          .number("basis_trace_GeV",
                  dataset_trace(*built->model, theta, ds, basis));
      jan.put("names", json::array(names))
          .put("value", json::array(Vec(theta(ps.free_indices()))))
          .put("sigma", json::array(U.sigma))
          .put("corr", json::array(U.rho))
          .number("condition", U.condition)
          .integer("converged", ff.fit.converged)
          .integer("iterations", ff.fit.iterations)
          .put("leverage", json::array(loo.leverage))
          .put("loo_pull", json::array(loo.pull));

      // Standardised out-of-sample misses t = miss / sqrt(sigma_th^2 w^2 +
      // s_exp^2 + sigma_par^2), every term from the refit without the
      // held-out data. With the scale estimated on nu_left degrees of
      // freedom, t ~ t_{nu_left} if the error model is right: E[t^2] =
      // nu/(nu - 2) (infinite for nu <= 2), coverage from t_coverage.
      const Vec s_exp = prob->exp_errors(), w_fit = prob->floor_weights();
      auto calibration = [&](const Vec& miss, const Vec& sth, const Vec& spar,
                             const std::vector<int>& fit_pos,
                             const std::vector<int>& nu_left, Vec& t) {
        t = Vec::Constant(miss.size(), nan);
        int n = 0, in1 = 0, in2 = 0, nu_lo = 1 << 30, nu_hi = 0;
        double t2 = 0, et2 = 0, e1 = 0, e2 = 0;
        for (int i = 0; i < miss.size(); ++i) {
          const int f = fit_pos[i], nl = nu_left[i];
          const double v = sth(i) * sth(i) * w_fit(f) * w_fit(f) +
                           s_exp(f) * s_exp(f) + spar(i) * spar(i);
          if (!std::isfinite(miss(i)) || !std::isfinite(v) || v <= 0 || nl < 1)
            continue;
          t(i) = miss(i) / std::sqrt(v);
          ++n, t2 += t(i) * t(i), in1 += std::fabs(t(i)) < 1,
              in2 += std::fabs(t(i)) < 2;
          et2 += nl > 2 ? nl / (nl - 2.0)
                        : std::numeric_limits<double>::infinity();
          e1 += t_coverage(1, nl), e2 += t_coverage(2, nl);
          nu_lo = std::min(nu_lo, nl), nu_hi = std::max(nu_hi, nl);
        }
        std::ostringstream o;
        if (n)
          o << "RMS t = " << fmt(std::sqrt(t2 / n), 0, 2) << " (expected "
            << (std::isfinite(et2) ? fmt(std::sqrt(et2 / n), 0, 2)
                                   : std::string("none: nu <= 2"))
            << "); |t| < 1: " << in1 << "/" << n << " (expected "
            << fmt(100 * e1 / n, 0, 0) << "%), |t| < 2: " << in2 << "/" << n
            << " (expected " << fmt(100 * e2 / n, 0, 0)
            << "%); scale from nu = " << nu_lo
            << (nu_hi > nu_lo ? "-" + std::to_string(nu_hi) : "");
        return o.str();
      };
      // ---- exact leave-one-out (refit without each state, floor re-solved)
      if (opt.exact_loo) {
        const LooRefits lr =
            exact_leave_one_out_refits(*prob, ff.fit.u, {}, iter_floor);
        const Vec& miss = lr.miss;
        std::vector<int> pos(static_cast<std::size_t>(miss.size()));
        for (std::size_t i = 0; i < pos.size(); ++i)
          pos[i] = static_cast<int>(i);
        Vec t;
        const std::string cal =
            calibration(miss, lr.sigma_th, lr.sigma_par, pos,
                        std::vector<int>(pos.size(), nu - 1), t);
        out << "\nexact leave-one-out (refit without the state"
            << (iter_floor ? ", floor re-solved" : "")
            << ", predict it): RMS = "
            << fmt(1e3 * std::sqrt(miss.squaredNorm() / miss.size()), 0, 2)
            << " MeV, max |miss| = "
            << fmt(1e3 * miss.cwiseAbs().maxCoeff(), 0, 2) << " MeV\n"
            << "  standardised misses: " << cal << "\n";
        jan.put("loo_exact_miss_GeV", json::array(miss))
            .put("loo_sigma_th_GeV", json::array(lr.sigma_th))
            .put("loo_sigma_par_GeV", json::array(lr.sigma_par))
            .put("loo_t", json::array(t));
        // delete-one jackknife of the free parameters (physical units)
        Mat freecols(lr.theta.rows(), p);
        for (int k = 0; k < p; ++k)
          freecols.col(k) = lr.theta.col(ps.free_indices()[k]);
        const Vec jk = jackknife_sigma(freecols);
        out << "jackknife errors (from the leave-one-out refits):";
        for (int k = 0; k < p; ++k)
          out << " " << names[k] << " " << fmt(jk(k), 0, 5);
        out << "\n";
        jan.put("jackknife_sigma", json::array(jk));
      }
      // ---- leave-multiplet-out: all states of one (n, L) multiplet removed
      // together and predicted (needs at least p data left)
      if (opt.lmo) {
        std::map<std::string, std::vector<int>> groups;
        const auto& dd = prob->data().d;
        for (std::size_t k = 0; k < dd.size(); ++k)
          if (dd[k].in_fit && dd[k].state.channel.size() > 1)
            groups[std::to_string(dd[k].state.n + 1) + dd[k].state.channel[1]]
                .push_back(static_cast<int>(k));
        // position of each datum among the fitted data
        std::map<int, int> fit_pos;
        for (std::size_t k = 0, i = 0; k < dd.size(); ++k)
          if (dd[k].in_fit)
            fit_pos[static_cast<int>(k)] = static_cast<int>(i++);
        std::string arr = "[";
        double ss = 0;
        int nm = 0;
        std::vector<double> all_miss, all_sth, all_spar;
        std::vector<int> all_pos, all_nu;
        out << "\nleave-multiplet-out (refit without the whole multiplet"
            << (iter_floor ? ", floor re-solved" : "") << "):\n";
        for (const auto& [gname, idx] : groups) {
          if (prob->n_data() - static_cast<int>(idx.size()) < p) continue;
          const GroupRefit gr =
              leave_group_out(*prob, ff.fit.u, idx, {}, iter_floor);
          std::vector<std::string> labels;
          for (int k : idx) labels.push_back(clean_label(dd[k].label));
          out << "  " << std::left << std::setw(4) << gname << std::right;
          for (int g = 0; g < gr.miss.size(); ++g) {
            out << "  " << labels[g] << " " << fmt(1e3 * gr.miss(g), 0, 1);
            if (std::isfinite(gr.miss(g))) {
              ss += gr.miss(g) * gr.miss(g);
              ++nm;
            }
            all_miss.push_back(gr.miss(g));
            all_sth.push_back(gr.sigma_th);
            all_spar.push_back(gr.sigma_par(g));
            all_pos.push_back(fit_pos.at(idx[g]));
            all_nu.push_back(nu - static_cast<int>(idx.size()));
          }
          out << " MeV; floor of the refit "
              << (std::isfinite(gr.sigma_th)
                      ? fmt(1e3 * gr.sigma_th, 0, 2) + " MeV"
                      : std::string("undetermined (nu = 0)"))
              << "\n";
          json::Object o;
          o.str("group", gname)
              .put("labels", json::array(labels))
              .put("miss_GeV", json::array(gr.miss))
              .put("sigma_par_GeV", json::array(gr.sigma_par))
              .number("sigma_th_GeV", gr.sigma_th);
          arr += (arr.size() > 1 ? ", " : "") + o.render(6);
        }
        arr += "]";
        if (nm) {
          auto vec = [](const std::vector<double>& v) {
            return Vec(
                Eigen::Map<const Vec>(v.data(), static_cast<int>(v.size())));
          };
          Vec t;
          out << "  RMS = " << fmt(1e3 * std::sqrt(ss / nm), 0, 2)
              << " MeV over " << nm << " states; standardised misses: "
              << calibration(vec(all_miss), vec(all_sth), vec(all_spar),
                             all_pos, all_nu, t)
              << "\n";
        }
        jan.put("lmo", arr);
      }
      // ---- random-split validation: repeated K-fold and the learning curve.
      // Each split refits without the held-out states (needs >= p left) and
      // scores the held-out (test) and the kept (train) states.
      if (opt.kfold > 1 || opt.learning_curve) {
        const auto& dd = prob->data().d;
        std::vector<int> fit_idx;
        for (std::size_t k = 0; k < dd.size(); ++k)
          if (dd[k].in_fit) fit_idx.push_back(static_cast<int>(k));
        const int N = static_cast<int>(fit_idx.size());
        // one random stream per method: results do not depend on which
        // methods are switched on
        Rng rng_kfold(opt.seed + 7, std::hash<std::string>{}(sname + "/kfold"));
        Rng rng_curve(opt.seed + 7, std::hash<std::string>{}(sname + "/curve"));
        auto shuffled = [&](Rng& rng) {
          std::vector<int> v = fit_idx;
          for (int i = N - 1; i > 0; --i)
            std::swap(v[i],
                      v[static_cast<int>(rng.uniform() * (i + 1)) % (i + 1)]);
          return v;
        };
        // One split: refit without `test`; squared errors of the kept (train)
        // and held-out (test) states, and the split's own test RMS (for a
        // median that is robust to splits that lose a whole multiplet).
        std::vector<double> per_split;
        auto split = [&](const std::vector<int>& test, double& tr, int& ntr,
                         double& te, int& nte) {
          const auto pred = refit_without(*prob, ff.fit.u, test);
          if (!pred) return;
          double se = 0;
          int ns = 0;
          for (int k : fit_idx) {
            const double e = dd[k].y - (*pred)(k);
            if (!std::isfinite(e)) continue;
            if (std::find(test.begin(), test.end(), k) != test.end()) {
              te += e * e, ++nte, se += e * e, ++ns;
            } else {
              tr += e * e, ++ntr;
            }
          }
          if (ns) per_split.push_back(std::sqrt(se / ns));
        };
        auto median = [](std::vector<double> v) {
          if (v.empty()) return std::numeric_limits<double>::quiet_NaN();
          std::sort(v.begin(), v.end());
          const std::size_t m = v.size() / 2;
          return v.size() % 2 ? v[m] : 0.5 * (v[m - 1] + v[m]);
        };
        auto root = [](double s, int n) {
          return n ? std::sqrt(s / n)
                   : std::numeric_limits<double>::quiet_NaN();
        };
        if (opt.kfold > 1) {
          const int K = std::min(opt.kfold, N), R = 10;
          double tr = 0, te = 0;
          int ntr = 0, nte = 0;
          per_split.clear();
          for (int r = 0; r < R; ++r) {
            const std::vector<int> v = shuffled(rng_kfold);
            for (int f = 0; f < K; ++f) {
              std::vector<int> test;
              for (int i = f; i < N; i += K) test.push_back(v[i]);
              if (N - static_cast<int>(test.size()) < p) continue;
              split(test, tr, ntr, te, nte);
            }
          }
          const double trr = root(tr, ntr), ter = root(te, nte),
                       tmed = median(per_split);
          out << "\n"
              << K << "-fold cross-validation (" << R
              << " random splits): train RMS = " << fmt(1e3 * trr, 0, 2)
              << " MeV, test RMS = " << fmt(1e3 * ter, 0, 2)
              << " MeV (median over folds " << fmt(1e3 * tmed, 0, 2) << ")\n";
          json::Object cv;
          cv.integer("K", K)
              .integer("repeats", R)
              .number("train_rms_GeV", trr)
              .number("test_rms_GeV", ter)
              .number("test_median_GeV", tmed);
          jan.put("kfold", cv.render(6));
        }
        if (opt.learning_curve) {
          const int nm = std::max(0, N - 1 - p);
          Vec sizes(nm), trs(nm), tes(nm), tmd(nm);
          out << "learning curve (RMS of model - PDG, MeV; 20 random subsets "
                 "per size):\n  fitted states    train     test   test "
                 "median\n";
          for (int m = p + 1, j = 0; m <= N - 1; ++m, ++j) {
            double tr = 0, te = 0;
            int ntr = 0, nte = 0;
            per_split.clear();
            for (int r = 0; r < 20; ++r) {
              const std::vector<int> v = shuffled(rng_curve);
              split(std::vector<int>(v.begin() + m, v.end()), tr, ntr, te, nte);
            }
            sizes(j) = m;
            trs(j) = root(tr, ntr);
            tes(j) = root(te, nte);
            tmd(j) = median(per_split);
            out << std::setw(15) << m << fmt(1e3 * trs(j), 9, 1)
                << fmt(1e3 * tes(j), 9, 1) << fmt(1e3 * tmd(j), 13, 1) << "\n";
          }
          json::Object lc;
          lc.put("train_size", json::array(sizes))
              .put("train_rms_GeV", json::array(trs))
              .put("test_rms_GeV", json::array(tes))
              .put("test_median_GeV", json::array(tmd));
          jan.put("learning_curve", lc.render(6));
        }
      }
      // ---- parametric bootstrap
      if (opt.bootstrap > 0) {
        const BootstrapResult br = bootstrap_fit(*prob, ff.fit.u, opt.bootstrap,
                                                 opt.seed, opt.threads);
        Mat ok(br.u.rows(), br.u.cols());
        int n = 0;
        for (int r = 0; r < br.u.rows(); ++r)
          if (br.u.row(r).allFinite()) ok.row(n++) = br.u.row(r);
        ok.conservativeResize(n, Eigen::NoChange);
        json::Object jb;
        jb.integer("failed", br.failed)
            .put("summary",
                 summarize_samples(ok, ps, "parametric bootstrap", out));
        jan.put("bootstrap", jb.render(4));
      }
      // ---- Bayesian posterior: the parameters and the theory floor
      // (marginalised, scale-invariant prior), with posterior predictive checks
      if (opt.mcmc > 0) {
        const McmcResult mr =
            sample_posterior_floor(*prob, ff.fit, ff.sigma_th, ps, 4, opt.mcmc,
                                   opt.mcmc / 4, opt.seed);
        int n = 0;
        for (const Chain& c : mr.chains)
          n += static_cast<int>(c.samples.rows());
        post.resize(n, ps.n_free() + 1);
        int r0 = 0;
        for (const Chain& c : mr.chains) {
          post.middleRows(r0, c.samples.rows()) = c.samples;
          r0 += static_cast<int>(c.samples.rows());
        }
        json::Object jm;
        jm.put("summary",
               summarize_samples(post.leftCols(ps.n_free()), ps,
                                 "Bayesian posterior (4 chains; flat in u, "
                                 "theory floor marginalised)",
                                 out))
            .put("rhat", json::array(mr.rhat))
            .put("n_eff", json::array(mr.n_eff));
        std::vector<double> sth(static_cast<std::size_t>(n));
        for (int r = 0; r < n; ++r) sth[r] = std::exp(post(r, ps.n_free()));
        const double s16 = percentile(sth, 0.1587), s50 = percentile(sth, 0.5),
                     s84 = percentile(sth, 0.8413);
        const PosteriorCheck pc = posterior_predictive(*prob, post);
        out << "  sigma_th (posterior median, 68%): " << fmt(1e3 * s50, 0, 2)
            << " MeV [" << fmt(1e3 * s16, 0, 2) << ", " << fmt(1e3 * s84, 0, 2)
            << "]\n  R-hat (parameters, ln sigma_th):";
        for (int k = 0; k < mr.rhat.size(); ++k)
          out << " " << fmt(mr.rhat(k), 0, 3);
        out << "   n_eff:";
        for (int k = 0; k < mr.n_eff.size(); ++k)
          out << " " << fmt(mr.n_eff(k), 0, 0);
        out << "\n  posterior predictive p (masses, " << pc.draws
            << " draws): T = chi2 " << fmt(pc.ppp_chi2, 0, 3)
            << ", T = max |z| " << fmt(pc.ppp_max, 0, 3)
            << " (small: the data are more discrepant than the model's own "
               "replicas)\n";
        jm.number("sigma_th_median_GeV", s50)
            .number("sigma_th_q16_GeV", s16)
            .number("sigma_th_q84_GeV", s84)
            .number("ppp_chi2", pc.ppp_chi2)
            .number("ppp_max", pc.ppp_max)
            .integer("ppp_draws", pc.draws);
        jan.put("mcmc", jm.render(4));
      }
      // ---- independent minimiser (Minuit2 Hesse/Minos) at the same floor
      if (opt.cross_check)
        jan.put("cross_check", opt.cross_check(*prob, ff.fit, ps, out));
    } else {
      if (!fitted.count("ccbar") || !fitted.count("bbar")) {
        out << "B_c needs the charmonium and bottomonium fits (use --systems "
               "ccbar,bbar,bc)\n\n";
        continue;
      }
      built = fitted["ccbar"].second;
      theta = interpolate_bc(fitted["ccbar"].first, fitted["bbar"].first);
      out << "predicted from the " << (joint_done ? "joint fit" : "two fits")
          << ": m_c = " << fmt(theta(0), 0, 4)
          << ", m_b = " << fmt(theta(1), 0, 4)
          << "; the other parameters interpolate in ln(reduced mass) between "
             "the charmonium and bottomonium values\n";
    }

    // ---- every PDG state
    // ----------------------------------------------------------------------------------------
    MatrixCache cache(basis);
    Vec snum = Vec::Zero(ds.d.size());
    Vec mu(ds.d.size()), spar = Vec::Zero(ds.d.size());
    Mat gr, pred_cov;
    if (!is_bc) {
      Chi2Problem full(*built->model, ds, ps, MatrixCache(basis), sigma_th);
      auto pg = full.predictions_with_gradients(theta);
      if (!pg)
        throw std::runtime_error("prediction failed: " + pg.error().message);
      mu = pg->first;
      for (int i = 0; i < mu.size(); ++i)
        spar(i) = std::sqrt(std::max(
            0.0,
            propagate_variance(Vec(pg->second.row(i).transpose()), Cfull)));
      // parametric covariance of every prediction (correlated splittings,
      // ratios and differences need it, not only spar)
      pred_cov = pg->second * Cfull * pg->second.transpose();
      snum = numerical_errors(*prob, theta, basis);
      // numerical_errors covers all data of the problem
    } else {
      for (std::size_t i = 0; i < ds.d.size(); ++i) {
        const Solved s = solve_state(*built->model, theta, cache, ds.d[i]);
        // the interpolated parameters can leave the physical domain (see
        // Chi2Problem): report a collapsed channel instead of its level
        const bool collapsed = s.sc.level(0).M < 0;
        if (collapsed)
          out << "warning: the B_c channel collapses at the interpolated "
                 "parameters (no prediction)\n";
        mu(i) = collapsed ? std::numeric_limits<double>::quiet_NaN()
                          : s.sc.level(s.level).M;
      }
    }

    // per-state diagnostics at the fitted parameters (not for the B_c):
    // <v^2> = <p^2>/m1^2, r_rms (fm), <p^4>(1/8m1^3 + 1/8m2^3) (GeV)
    Mat diag = Mat::Constant(static_cast<int>(ds.d.size()), 3,
                             std::numeric_limits<double>::quiet_NaN());
    std::vector<WaveFunction> wfs(ds.d.size());
    if (!is_bc)
      for (std::size_t i = 0; i < ds.d.size(); ++i) {
        try {
          const Solved st = solve_state(*built->model, theta, cache, ds.d[i]);
          const double m1 = theta(0), m2 = theta(1);
          diag(i, 0) = p2_expectation(st.sc, st.level, cache) / (m1 * m1);
          diag(i, 1) =
              std::sqrt(moment(st.sc, st.level, cache, 2)) * hbarc_GeV_fm;
          diag(i, 2) = p4_expectation(st.sc, st.level, cache) *
                       (1 / (8 * m1 * m1 * m1) + 1 / (8 * m2 * m2 * m2));
          wfs[i] = wave_function(st.sc, st.level, cache);
        } catch (const std::exception&) {
        }
      }
    // Parametric band of every wave function: u(r) = r R(r) of the dominant
    // component for 40 draws u ~ N(u_hat, C_u), each draw's sign aligned with
    // the central curve; 16/50/84 % quantiles on 100 points up to 4 r_rms.
    std::vector<std::string> wband(ds.d.size());
    if (!is_bc && best_fit.J.size() && U.C.size() && U.condition < 1e10) {
      Eigen::SelfAdjointEigenSolver<Mat> es(U.C);
      const Mat A = es.eigenvectors() *
                    es.eigenvalues().cwiseMax(0.0).cwiseSqrt().asDiagonal();
      Rng rng(opt.seed + 11, std::hash<std::string>{}(sname + "/wfband"));
      const int K = 40, G = 100;
      std::vector<Vec> draws;
      for (int k = 0; k < K; ++k) {
        Vec xi(A.cols());
        for (int j = 0; j < xi.size(); ++j) xi(j) = rng.normal();
        draws.push_back(ps.theta(Vec(best_fit.u + A * xi)));
      }
      for (std::size_t i = 0; i < ds.d.size(); ++i) {
        const WaveFunction& w0 = wfs[i];
        if (w0.c.empty() || !std::isfinite(diag(static_cast<int>(i), 1)))
          continue;
        const int L = w0.L[w0.dominant];
        const double rmax = 4 * diag(static_cast<int>(i), 1) / hbarc_GeV_fm;
        Vec r(G), u0(G);
        for (int g = 0; g < G; ++g) {
          r(g) = rmax * (g + 1) / G;
          u0(g) = r(g) * radial_value(cache.basis(L), w0.c[w0.dominant], r(g));
        }
        std::vector<Vec> us;
        for (const Vec& th : draws) {
          try {
            const Solved sd = solve_state(*built->model, th, cache, ds.d[i]);
            const WaveFunction wk = wave_function(sd.sc, sd.level, cache);
            int dk = 0;  // the component with the dominant L
            for (std::size_t q = 0; q < wk.L.size(); ++q)
              if (wk.L[q] == L) dk = static_cast<int>(q);
            Vec u(G);
            for (int g = 0; g < G; ++g)
              u(g) = r(g) * radial_value(cache.basis(L), wk.c[dk], r(g));
            if (u.dot(u0) < 0) u = -u;
            us.push_back(u);
          } catch (const std::exception&) {
          }
        }
        if (us.size() < 10) continue;
        Vec lo(G), mid(G), hi(G);
        for (int g = 0; g < G; ++g) {
          std::vector<double> v;
          for (const Vec& u : us) v.push_back(u(g));
          lo(g) = percentile(v, 0.1587), mid(g) = percentile(v, 0.5),
          hi(g) = percentile(v, 0.8413);
        }
        json::Object jb;
        jb.put("r_fm", json::array(Vec(r * hbarc_GeV_fm)))
            .put("u_central", json::array(u0))
            .put("u_q16", json::array(lo))
            .put("u_q50", json::array(mid))
            .put("u_q84", json::array(hi))
            .integer("draws", static_cast<int>(us.size()));
        wband[i] = jb.render(6);
      }
    }
    // Truncation test (A/B with the non-relativistic kinetic energy): the
    // first omitted term -<p^4>(1/8m1^3 + 1/8m2^3) is added to the model's own
    // values and refitted; what the parameters cannot absorb, R, is the
    // a-priori misfit of the truncation, and sigma_a its floor.
    const bool omits_p4 =
        !is_bc && (opt.model == "A" || opt.model == "B") && opt.kinetic == "nr";
    Truncation tr;
    std::string tr_skip;
    if (omits_p4) {
      Vec shift = diag.col(2);
      bool fit_ok = best_fit.J.size() > 0;
      for (std::size_t i = 0; i < ds.d.size(); ++i)
        if (!std::isfinite(shift(static_cast<int>(i)))) {
          if (ds.d[i].in_fit) fit_ok = false;
          shift(static_cast<int>(i)) = 0;  // a prediction without <p^4>
        }
      if (!fit_ok)
        tr_skip = "<p^4> unavailable for a fitted state";
      else if (tr = absorb_shift(*prob, best_fit.u, shift); !tr.R.size())
        tr_skip = "the refit of the shifted values failed";
      for (std::size_t i = 0; i < ds.d.size() && tr.R.size(); ++i)
        if (!std::isfinite(diag(static_cast<int>(i), 2)))
          tr.R(static_cast<int>(i)) = std::numeric_limits<double>::quiet_NaN();
    }
    out << "\nstate          channel     PDG (GeV) model (GeV) M"
           "-PDG(MeV)  par(MeV)  num(MeV)     pull lever.    role\n";
    json::Object jsys;
    std::string jrows = "[";
    int fit_row = 0;
    for (std::size_t i = 0; i < ds.d.size(); ++i) {
      const Datum& x = ds.d[i];
      const double diff = 1e3 * (mu(i) - x.y);
      double pull = 0, h = 0;
      const double floor_i = sigma_th * x.floor_w;
      const double sig =
          std::sqrt(spar(i) * spar(i) + x.err * x.err + floor_i * floor_i);
      pull = sig > 0 ? (x.y - mu(i)) / sig : 0;
      const bool show_pull = !is_bc;
      if (x.in_fit && !is_bc) {
        h = loo.leverage(fit_row);
        pull = loo.pull(fit_row);
        ++fit_row;
      }
      out << std::left << std::setw(15) << clean_label(x.label) << std::setw(10)
          << state_name(x) << std::right << fmt(x.y, 11, 5) << fmt(mu(i), 12, 5)
          << fmt(diff, 11, 1)
          << (is_bc ? "         -" : fmt(1e3 * spar(i), 10, 1))
          << fmt(1e3 * snum(i), 10, 3)
          << (show_pull ? fmt(pull, 9, 2) : "        -")
          << (x.in_fit ? fmt(h, 7, 2) : "      -") << "    "
          << (x.in_fit ? "fit" : "prediction") << "\n";
      json::Object r;
      r.str("label", clean_label(x.label))
          .str("channel", x.state.channel)
          .integer("n", x.state.n)
          .number("pdg_GeV", x.y)
          .number("model_GeV", mu(i))
          .number("sigma_par_GeV", spar(i))
          .number("sigma_num_GeV", snum(i))
          .number("sigma_exp_GeV", std::hypot(x.err, x.sys_common))
          .number("floor_w", x.floor_w)
          .number("pull", pull)
          .integer("in_fit", x.in_fit);
      if (!is_bc) {
        r.number("v2", diag(i, 0))
            .number("r_rms_fm", diag(i, 1))
            .number("p4_term_GeV", diag(i, 2));
        if (tr.R.size()) r.number("trunc_GeV", tr.R(static_cast<int>(i)));
        if (const WaveFunction& w = wfs[i]; !w.c.empty()) {
          std::string comps = "[";
          for (std::size_t k = 0; k < w.c.size(); ++k) {
            json::Object jc;
            jc.integer("L", w.L[k])
                .integer("S", w.S[k])
                .number("weight", w.weight[k])
                .put("c", json::array(w.c[k]));
            comps += (k ? ", " : "") + jc.render(8);
          }
          json::Object jw;
          jw.number("norm", w.norm)
              .integer("nodes", w.nodes)
              .put("components", comps + "]");
          r.put("wavefunction", jw.render(6));
        }
        if (!wband[i].empty()) r.put("wavefunction_band", wband[i]);
      }
      jrows += (i ? ", " : "") + r.render(4);
    }
    if (!is_bc) {
      // internal consistency of the potential description, per state: the
      // velocity of each quark (the expansion parameter of the
      // non-relativistic and Breit-Fermi terms), the size (linear confinement
      // and its string breaking act at large r), and the first neglected
      // kinetic term <p^4>(1/8m1^3 + 1/8m2^3), an a-priori truncation error
      out << "  size and velocity: <v^2> = <p^2>/m_Q^2; r_rms; neglected "
             "kinetic term <p^4>/8m^3 per quark (MeV)\n";
      for (std::size_t i = 0; i < ds.d.size(); ++i)
        out << "    " << std::left << std::setw(15)
            << clean_label(ds.d[i].label) << std::right
            << " <v^2> = " << fmt(diag(i, 0), 5, 3)
            << "   r_rms = " << fmt(diag(i, 1), 5, 2)
            << " fm   p4 term = " << fmt(1e3 * diag(i, 2), 6, 1) << "\n";
    }
    jrows += "]";
    if (!is_bc) {
      // wave functions: components and weights, norm, radial nodes (must
      // equal the radial quantum number), and the expansion coefficients
      static const char* kL = "SPDFGH";
      out << "  wave functions (R_L(r) = sum_n c_n N_nL r^L exp(-nu_n r^2), "
          << basis.n << " Gaussians, ranges " << basis.r1 << " ... "
          << basis.rmax << " GeV^-1): components, norm, radial nodes\n";
      int mismatch = 0;
      for (std::size_t i = 0; i < ds.d.size(); ++i) {
        const WaveFunction& w = wfs[i];
        if (w.c.empty()) continue;
        out << "    " << std::left << std::setw(15)
            << clean_label(ds.d[i].label) << std::right;
        for (std::size_t k = 0; k < w.c.size(); ++k)
          out << (k ? " + " : " ") << 2 * w.S[k] + 1 << kL[w.L[k]] << " "
              << fmt(100 * w.weight[k], 0, 2) << "%";
        const bool ok = w.nodes == ds.d[i].state.n;
        mismatch += !ok;
        out << "; norm " << fmt(w.norm, 0, 8) << "; nodes " << w.nodes
            << (ok ? ""
                   : "  (EXPECTED " + std::to_string(ds.d[i].state.n) + ")")
            << "\n";
        if (opt.wavefunctions)
          for (std::size_t k = 0; k < w.c.size(); ++k) {
            out << "      c(" << 2 * w.S[k] + 1 << kL[w.L[k]] << ") =";
            for (int j = 0; j < w.c[k].size(); ++j)
              out << " " << std::setprecision(6) << w.c[k](j);
            out << std::defaultfloat << "\n";
          }
      }
      out << "  radial nodes: "
          << (mismatch ? std::to_string(mismatch) +
                             " state(s) differ from the "
                             "radial quantum number"
                       : std::string("all equal the radial "
                                     "quantum number"))
          << "\n";
      // the basis, for reconstructing the wave functions from the JSON
      json::Object jb;
      std::string norms = "{";
      for (int L = 0; L <= 3; ++L) {
        const Basis& b = cache.basis(L);
        Vec nl(b.size());
        for (int j = 0; j < b.size(); ++j) nl(j) = b.norm(j);
        norms += (L ? ", " : "") + std::string("\"") + std::to_string(L) +
                 "\": " + json::array(nl);
      }
      jb.integer("n", basis.n)
          .number("r1_GeVinv", basis.r1)
          .number("rmax_GeVinv", basis.rmax)
          .put("nu_GeV2", json::array(cache.basis(0).nu))
          .put("norm_by_L", norms + "}");
      jsys.put("basis", jb.render(4));
    }
    if (omits_p4 && !tr.R.size())
      out << "  truncation test skipped: " << tr_skip << "\n";
    if (tr.R.size()) {
      // data - model against R: the omitted term explains the misfit if the
      // residual pattern is R (lambda = 1) and chi2/nu ~ 1 at sigma_a.
      std::vector<int> fit_idx, pred_idx;
      for (std::size_t i = 0; i < ds.d.size(); ++i)
        (ds.d[i].in_fit ? fit_idx : pred_idx).push_back(static_cast<int>(i));
      const int nf = static_cast<int>(fit_idx.size());
      Vec Rf(nf), D(nf), Pf(nf);
      for (int i = 0; i < nf; ++i) {
        Rf(i) = tr.R(fit_idx[i]);
        D(i) = ds.d[fit_idx[i]].y - mu(fit_idx[i]);
        Pf(i) = diag(fit_idx[i], 2);
      }
      const Eigen::LLT<Mat> L(prob->covariance());
      const Vec Rw = L.matrixL().solve(Rf), Dw = L.matrixL().solve(D);
      const double lam = Rw.dot(Dw) / Rw.squaredNorm();
      const double lam_err = 1 / Rw.norm();
      Chi2Problem at(*prob);  // the fit at the a-priori floor
      at.set_sigma_th(tr.sigma_a);
      const FitResult fa = levenberg_marquardt(at, best_fit.u);
      const double chi2_a =
          fa.J.size() ? fa.chi2 : std::numeric_limits<double>::quiet_NaN();
      auto rms_of = [&](const std::vector<int>& idx, auto f) {
        double s2 = 0;
        int n = 0;
        for (int i : idx)
          if (std::isfinite(f(i))) s2 += f(i) * f(i), ++n;
        return n ? std::sqrt(s2 / n) : std::numeric_limits<double>::quiet_NaN();
      };
      out << "  truncation test (first omitted term -<p^4>(1/8m1^3 + 1/8m2^3), "
             "absorbed by a refit):\n"
          << "    omitted term over the fitted states: "
          << fmt(1e3 * Pf.minCoeff(), 0, 1) << " to "
          << fmt(1e3 * Pf.maxCoeff(), 0, 1) << " MeV; unabsorbed part R: RMS "
          << fmt(1e3 * Rf.norm() / std::sqrt(nf), 0, 2) << " MeV\n"
          << "    a-priori floor sigma_a = " << fmt(1e3 * tr.sigma_a, 0, 2)
          << " MeV (sum R^2/(s^2 + sigma_a^2) = nu); fitted sigma_th = "
          << fmt(1e3 * sigma_th, 0, 2) << " MeV\n"
          << "    chi2 at sigma_a (refit) = " << fmt(chi2_a, 0, 2)
          << " for nu = " << nu << ": chi2/nu = " << fmt(chi2_a / nu, 0, 2)
          << ", P(chi2_nu >= chi2) = " << fmt_p(chi2_sf(chi2_a, nu)) << "\n"
          << "    pattern: (PDG - model) = lambda R, lambda = "
          << fmt(lam, 0, 2) << " +- " << fmt(lam_err, 0, 2)
          << " (1: the omitted term is the misfit; 0: unrelated)\n"
          << "    predictions: RMS of R "
          << fmt(1e3 * rms_of(pred_idx, [&](int i) { return tr.R(i); }), 0, 1)
          << " MeV against their RMS miss "
          << fmt(1e3 *
                     rms_of(pred_idx, [&](int i) { return mu(i) - ds.d[i].y; }),
                 0, 1)
          << " MeV\n";
      json::Object jt;
      jt.number("sigma_a_GeV", tr.sigma_a)
          .number("chi2", chi2_a)
          .integer("nu", nu)
          .number("p", chi2_sf(chi2_a, nu))
          .number("lambda", lam)
          .number("lambda_err", lam_err);
      jan.put("truncation", jt.render(6));
    }
    if (!is_bc)
      out << "  par: parametric error (error propagation); num: basis error "
             "(n -> 2n-1); lever.: leverage h = dmu_i/dy_i of a fitted state\n"
             "  pull: (PDG - model) / combined error; for fitted states the "
             "leave-one-out pull, i.e. the state predicted without itself\n";

    // ---- widths from the wave functions at the origin
    // ----------------------------------------------------------
    if (is_bc) {
      out << "\n(B_c has no annihilation channel into e+e- or gamma gamma; its "
             "prediction carries the systematic of the parameter "
             "interpolation)\n\n";
      jsys.str("system", sys.name)
          .str("model", built->model->name())
          .put("theta_names", json::array(built->model->param_names()))
          .put("theta", json::array(theta))
          .put("states", jrows);
      json_systems += (nsys++ ? ", " : "") + jsys.render(2);
      continue;
    }
    // Widths: LO and NLO, with the parametric error of NLO (central
    // differences in u, propagated with the fit's covariance) and the standard
    // theory error (half the scale band, in quadrature with <v^2> Gamma_NLO).
    const double eQ = sys.charge;
    const auto wch = width_channels(ds);
    const auto wv = compute_widths(*built->model, theta, ds, cache, eQ, wch);
    const int nwc = static_cast<int>(wch.size()), pf = ps.n_free();
    Mat Gw = Mat::Zero(nwc, pf);  // d Gamma_NLO / d u
    if (best_fit.J.size())
      for (int a = 0; a < pf; ++a) {
        const double h = 1e-5 * std::max(1.0, std::fabs(best_fit.u(a)));
        Vec up = best_fit.u, um = best_fit.u;
        up(a) += h, um(a) -= h;
        const auto wp =
            compute_widths(*built->model, ps.theta(up), ds, cache, eQ, wch);
        const auto wm =
            compute_widths(*built->model, ps.theta(um), ds, cache, eQ, wch);
        for (int k = 0; k < nwc; ++k)
          Gw(k, a) = (wp[k].nlo - wm[k].nlo) / (2 * h);
      }
    const Mat Vpar_w =
        U.C.size() ? Mat(Gw * U.C * Gw.transpose()) : Mat::Zero(nwc, nwc);
    std::vector<std::optional<PdgWidth>> pw(nwc);
    for (int k = 0; k < nwc; ++k)
      pw[k] = pdg.partial_width(ds.d[wch[k].datum].mcid, wch[k].process);
    out << "\nannihilation widths (keV; M = model mass): LO (Van Royen-"
           "Weisskopf); NLO = LO x one-loop QCD factor, MS-bar alpha_s at mu = "
           "M/2; errors: parametric, theory (scale band mu = M/4 .. M, half "
           "width, (+) <v^2> NLO); PDG (* = summary-table BR x total width); "
           "pull = (PDG - NLO)/total error\n";
    out << "state          process            LO        NLO     +-par      "
           "+-th   band [M/4, M]             PDG               pull\n";
    int nw = 0;
    std::string jw = "[";
    for (int k = 0; k < nwc; ++k) {
      const Datum& x = ds.d[wch[k].datum];
      const WidthValue& w = wv[k];
      if (!std::isfinite(w.nlo)) continue;
      const double par = std::sqrt(std::max(0.0, Vpar_w(k, k)));
      const double th = width_theory_error(w);
      out << std::left << std::setw(15) << clean_label(x.label) << std::setw(13)
          << wch[k].process << std::right << fmt(1e6 * w.lo, 9, 4)
          << fmt(1e6 * w.nlo, 11, 4) << fmt(1e6 * par, 10, 4)
          << fmt(1e6 * th, 10, 4) << "   ["
          << fmt(1e6 * std::min(w.nlo_lo, w.nlo_hi), 7, 4) << ", "
          << fmt(1e6 * std::max(w.nlo_lo, w.nlo_hi), 7, 4) << "]";
      json::Object jo;
      jo.str("label", clean_label(x.label))
          .str("process", wch[k].process)
          .number("lo_keV", 1e6 * w.lo)
          .number("nlo_keV", 1e6 * w.nlo)
          .number("nlo_err_par_keV", 1e6 * par)
          .number("nlo_err_th_keV", 1e6 * th)
          .number("nlo_band_lo_keV", 1e6 * std::min(w.nlo_lo, w.nlo_hi))
          .number("nlo_band_hi_keV", 1e6 * std::max(w.nlo_lo, w.nlo_hi))
          .number("alpha_s_mu", alpha_s_msbar(w.M / 2))
          .number("v2", w.v2)
          .number("M_GeV", w.M)
          .integer("in_fit", x.in_fit);
      if (pw[k]) {
        const double tot = std::sqrt(std::pow(pw[k]->err_keV, 2) +
                                     1e12 * (par * par + th * th));
        out << fmt(pw[k]->value_keV, 10, 4) << " +- "
            << fmt(pw[k]->err_keV, 6, 4) << (pw[k]->derived ? "*" : " ")
            << fmt((pw[k]->value_keV - 1e6 * w.nlo) / tot, 8, 2);
        jo.number("pdg_keV", pw[k]->value_keV)
            .number("pdg_err_keV", pw[k]->err_keV)
            .integer("pdg_derived", pw[k]->derived)
            .number("pull", (pw[k]->value_keV - 1e6 * w.nlo) / tot);
      } else {
        out << "          -";
      }
      out << "\n";
      jw += (nw ? ", " : "") + jo.render(6);
      ++nw;
    }
    if (nw == 0) out << "  (no annihilation channel among these states)\n";
    jw += "]";

    // Decay benchmark: the measured widths of the fitted states, which the
    // mass fit never saw; Pearson chi^2 with V = V_exp + V_theory (scale part
    // correlated within a process) + V_par (correlated through the fit).
    std::vector<int> bk;
    for (int k = 0; k < nwc; ++k)
      if (pw[k] && ds.d[wch[k].datum].in_fit && std::isfinite(wv[k].nlo))
        bk.push_back(k);
    if (!bk.empty()) {
      const int nb = static_cast<int>(bk.size());
      std::vector<WidthChannel> bch;
      std::vector<WidthValue> bwv;
      for (int k : bk) bch.push_back(wch[k]), bwv.push_back(wv[k]);
      Mat V = width_theory_covariance(bwv, bch);
      Vec D(nb);
      for (int i = 0; i < nb; ++i) {
        D(i) = 1e-6 * pw[bk[i]]->value_keV - wv[bk[i]].nlo;
        V(i, i) += std::pow(1e-6 * pw[bk[i]]->err_keV, 2);
        for (int j = 0; j < nb; ++j) V(i, j) += Vpar_w(bk[i], bk[j]);
      }
      const double chi2_dec = D.dot(V.llt().solve(D));
      std::vector<std::string> blab;
      Vec bpull(nb);
      for (int i = 0; i < nb; ++i) {
        blab.push_back(clean_label(ds.d[wch[bk[i]].datum].label));
        bpull(i) = D(i) / std::sqrt(V(i, i));
      }
      out << "decay benchmark (measured widths of the fitted states; the mass "
             "fit never saw them): chi2 = "
          << fmt(chi2_dec, 0, 2) << " for N = " << nb
          << ", chi2/N = " << fmt(chi2_dec / nb, 0, 2)
          << ", P(chi2_N >= chi2) = " << fmt_p(chi2_sf(chi2_dec, nb)) << "\n";
      json::Object jd;
      jd.number("chi2", chi2_dec)
          .integer("n", nb)
          .number("p", chi2_sf(chi2_dec, nb))
          .put("labels", json::array(blab))
          .put("pulls", json::array(bpull));
      // posterior predictive check of the decays: the widths of posterior
      // draws (parametric uncertainty from the draws), V = V_exp + V_theory
      if (post.rows() > 0) {
        double ppp = 0;
        int nd = 0;
        const int step = std::max<int>(1, static_cast<int>(post.rows()) / 100);
        for (int r = 0; r < post.rows(); r += step) {
          const Vec th = ps.theta(Vec(post.row(r).head(pf).transpose()));
          const auto w = compute_widths(*built->model, th, ds, cache, eQ, bch);
          Mat Vd = width_theory_covariance(w, bch);
          Vec Dd(nb);
          bool ok = Vd.allFinite();
          for (int i = 0; i < nb && ok; ++i) {
            Dd(i) = 1e-6 * pw[bk[i]]->value_keV - w[i].nlo;
            Vd(i, i) += std::pow(1e-6 * pw[bk[i]]->err_keV, 2);
            ok = std::isfinite(Dd(i));
          }
          if (!ok) continue;
          ppp += chi2_sf(Dd.dot(Vd.llt().solve(Dd)), nb);
          ++nd;
        }
        if (nd) {
          out << "  posterior predictive p (decays, " << nd
              << " draws): " << fmt(ppp / nd, 0, 3) << "\n";
          jd.number("ppp", ppp / nd).integer("ppp_draws", nd);
        }
      }
      jan.put("decay_benchmark", jd.render(6));

      // Combined fit of masses and these widths (--fit-widths)
      if (opt.fit_widths && best_fit.J.size()) {
        WidthData wd;
        wd.y.resize(nb), wd.s_exp.resize(nb);
        for (int i = 0; i < nb; ++i) {
          wd.ch.push_back(wch[bk[i]]);
          wd.y(i) = 1e-6 * pw[bk[i]]->value_keV;
          wd.s_exp(i) = 1e-6 * pw[bk[i]]->err_keV;
        }
        WidthProblem wpr(*built->model, ds, ps, MatrixCache(basis), eQ, wd,
                         best_fit.u);
        Chi2Problem pm(*prob);
        MassWidthProblem mw(pm, wpr);
        const MassWidthFit f = fit_masses_and_widths(mw, best_fit.u, sigma_th);
        if (!f.fit.J.size()) {
          out << "masses + widths fit failed: " << f.fit.message << "\n";
        } else {
          const Uncertainty Uc = analyse(f.fit, ps);
          const Vec tc = ps.theta(f.fit.u);
          out << "masses + widths fit (--fit-widths): mass floor sigma_th = "
              << fmt(1e3 * f.sigma_th, 0, 3) << " MeV (mass-only "
              << fmt(1e3 * sigma_th, 0, 3)
              << "); widths chi2 = " << fmt(f.chi2_w, 0, 2) << " for N = " << nb
              << ", P = " << fmt_p(chi2_sf(f.chi2_w, nb)) << "\n  parameters:";
          const auto nm = ps.free_names();
          Vec val(pf), shift(pf);
          for (int a = 0; a < pf; ++a) {
            val(a) = tc(ps.free_indices()[a]);
            shift(a) = (val(a) - theta(ps.free_indices()[a])) / U.sigma(a);
            out << " " << nm[a] << " = " << fmt(val(a), 0, 5) << " +- "
                << fmt(Uc.sigma(a), 0, 5) << " (" << fmt(shift(a), 0, 2)
                << " sigma from mass-only)";
          }
          out << "\n";
          json::Object jm;
          jm.number("sigma_th_GeV", f.sigma_th)
              .number("nu_mass", f.nu_m)
              .number("chi2_mass", f.chi2_m)
              .number("chi2_widths", f.chi2_w)
              .integer("n_widths", nb)
              .number("p_widths", chi2_sf(f.chi2_w, nb))
              .put("names", json::array(nm))
              .put("value", json::array(val))
              .put("sigma", json::array(Uc.sigma))
              .put("shift_sigma", json::array(shift));
          jan.put("mass_width_fit", jm.render(6));
        }
      }
    }
    out << "\n";
    jsys.str("system", sys.name)
        .str("model", built->model->name())
        .number("sigma_th_GeV", sigma_th)
        .number("chi2", chi2)
        .integer("n_data", prob->n_data())
        .integer("n_free", prob->n_free())
        .integer("nu", nu)
        .put("theta_names", json::array(built->model->param_names()))
        .put("analysis", jan.render(4))
        .put("widths", jw)
        .put("theta", json::array(theta))
        .put("theta_cov", json::array(Cfull))
        .put("pred_cov_GeV2", json::array(pred_cov))
        .put("V_exp_GeV2", json::array(ds.V_exp()))
        .put("states", jrows);
    json_systems += (nsys++ ? ", " : "") + jsys.render(2);
    kept[sname] = Kept{std::move(prob), best_fit.u, sigma_th, chi2};

    // ---- joint fit of c c-bar and b b-bar with shared parameters
    if (!opt.shared.empty() && !joint_done && kept.count("ccbar") &&
        kept.count("bbar")) {
      joint_done = true;
      Kept& kc = kept.at("ccbar");
      Kept& kb = kept.at("bbar");
      out << "================ joint fit of charmonium and bottomonium "
             "(shared:";
      for (const auto& n : opt.shared) out << " " << n;
      out << ") ================\n";
      Chi2Problem pc(*kc.prob), pb(*kb.prob);  // copies
      JointProblem jp(pc, pb, opt.shared);
      const bool iter = opt.fixed_floor < 0;
      // Universality (likelihood ratio): the joint fit at the separate fits'
      // floors against the separate fits, dchi2 ~ chi2_k for k shared. Two
      // starts: the shared values of either system; the lower chi2 is kept.
      pc.set_sigma_th(kc.sigma_th), pb.set_sigma_th(kb.sigma_th);
      FitResult fs;
      for (bool from_b : {false, true}) {
        const FitResult f =
            levenberg_marquardt(jp, jp.join(kc.u, kb.u, from_b));
        if (f.J.size() && (!fs.J.size() || f.chi2 < fs.chi2)) fs = f;
      }
      JointFit jf;
      if (fs.J.size()) {
        if (iter) {
          jf = fit_joint(jp, fs.u, kc.sigma_th, kb.sigma_th);
        } else {
          jf.fit = fs;
          jf.sigma_a = kc.sigma_th, jf.sigma_b = kb.sigma_th;
          const Vec h = leverages(fs.J);
          jf.nu_a = pc.n_data() - h.head(pc.n_data()).sum();
          jf.nu_b = pb.n_data() - h.tail(pb.n_data()).sum();
        }
      }
      if (!fs.J.size() || !jf.fit.J.size()) {
        out << "joint fit failed from both starts\n\n";
      } else {
        const double dchi2 = fs.chi2 - (kc.chi2 + kb.chi2);
        const double pval = chi2_sf(std::max(dchi2, 0.0), jp.k());
        pc.set_sigma_th(jf.sigma_a), pb.set_sigma_th(jf.sigma_b);
        const Vec u = jf.fit.u;
        const Vec tc = pc.params().theta(jp.u_a(u)),
                  tb = pb.params().theta(jp.u_b(u));
        // linear errors of the joint parameters
        Eigen::SelfAdjointEigenSolver<Mat> es(jf.fit.J.transpose() * jf.fit.J);
        const Vec lam = es.eigenvalues();
        Vec inv = Vec::Zero(lam.size());
        for (int i = 0; i < lam.size(); ++i)
          if (lam(i) > lam.maxCoeff() * 1e-14) inv(i) = 1 / lam(i);
        const Mat Cu = es.eigenvectors() * inv.asDiagonal() *
                       es.eigenvectors().transpose();
        const Vec d = jp.dtheta_du(u);
        const auto names = jp.names();
        Vec val(static_cast<int>(names.size())), err(val.size());
        const auto nc = pc.params().free_names();
        for (int j = 0; j < val.size(); ++j) {
          const bool is_b =
              names[j].size() > 3 &&
              names[j].compare(names[j].size() - 3, 3, "(b)") == 0;
          const std::string base =
              is_b ? names[j].substr(0, names[j].size() - 3) : names[j];
          val(j) =
              is_b ? tb(pb.params().index(base)) : tc(pc.params().index(base));
          err(j) = std::sqrt(std::max(0.0, Cu(j, j))) * std::fabs(d(j));
        }
        // RMS of model - PDG per system: fitted states and predictions
        auto rms_pair = [](Chi2Problem& pr, const Vec& th) {
          double sf = 0, sp = 0;
          int nf = 0, np = 0;
          if (const auto mu = pr.predictions(th))
            for (std::size_t k = 0; k < pr.data().d.size(); ++k) {
              const double e = (*mu)(static_cast<int>(k)) - pr.data().d[k].y;
              if (!std::isfinite(e)) continue;
              (pr.data().d[k].in_fit ? sf : sp) += e * e;
              ++(pr.data().d[k].in_fit ? nf : np);
            }
          const double nan = std::numeric_limits<double>::quiet_NaN();
          return std::make_pair(nf ? std::sqrt(sf / nf) : nan,
                                np ? std::sqrt(sp / np) : nan);
        };
        const auto rc = rms_pair(pc, tc), rb = rms_pair(pb, tb);
        out << "  " << pc.n_data() << " + " << pb.n_data() << " fitted states, "
            << jp.n_free() << " free parameters (" << jp.k() << " shared)\n"
            << "  floors: charmonium sigma_th = " << fmt(1e3 * jf.sigma_a, 0, 3)
            << " MeV (nu = " << fmt(jf.nu_a, 0, 2) << "), bottomonium "
            << fmt(1e3 * jf.sigma_b, 0, 3)
            << " MeV (nu = " << fmt(jf.nu_b, 0, 2) << "); separate fits "
            << fmt(1e3 * kc.sigma_th, 0, 3) << " and "
            << fmt(1e3 * kb.sigma_th, 0, 3) << " MeV\n"
            << "  (each floor solves chi2_s = N_s - sum of the leverages of "
               "its "
               "states)\n"
            << "  universality test: at the separate floors, chi2(joint) - "
               "chi2(separate) = "
            << fmt(dchi2, 0, 2) << " for " << jp.k() << " constraint"
            << (jp.k() > 1 ? "s" : "") << ", P(chi2_" << jp.k()
            << " >= dchi2) = " << fmt_p(pval) << "\n"
            << "  RMS model - PDG (MeV): charmonium fitted "
            << fmt(1e3 * rc.first, 0, 1) << ", predicted "
            << fmt(1e3 * rc.second, 0, 1) << "; bottomonium fitted "
            << fmt(1e3 * rb.first, 0, 1) << ", predicted "
            << fmt(1e3 * rb.second, 0, 1) << "\n"
            << "joint parameters ((b): bottomonium's own):\n";
        for (int j = 0; j < val.size(); ++j)
          out << "  " << std::left << std::setw(12) << names[j] << std::right
              << fmt(val(j), 11, 5) << " +- " << fmt(err(j), 9, 5) << "\n";
        out << "\n";
        fitted["ccbar"].first = tc;  // the B_c below uses the joint parameters
        fitted["bbar"].first = tb;
        json::Object jj;
        jj.put("shared", json::array(opt.shared))
            .put("names", json::array(names))
            .put("value", json::array(val))
            .put("sigma", json::array(err))
            .number("sigma_th_ccbar_GeV", jf.sigma_a)
            .number("sigma_th_bbar_GeV", jf.sigma_b)
            .number("nu_ccbar", jf.nu_a)
            .number("nu_bbar", jf.nu_b)
            .number("chi2", jf.fit.chi2)
            .number("dchi2", dchi2)
            .integer("k", jp.k())
            .number("p", pval)
            .number("ccbar_rms_fit_GeV", rc.first)
            .number("ccbar_rms_pred_GeV", rc.second)
            .number("bbar_rms_fit_GeV", rb.first)
            .number("bbar_rms_pred_GeV", rb.second);
        json_all.put("joint", jj.render(2));
      }
    }
  }
  json_systems += "]";
  if (!opt.json_out.empty()) {
    std::ofstream f(opt.json_out);
    json::Object cfg;
    cfg.str("model", opt.model)
        .str("running", opt.running)
        .integer("c_free", opt.with_c)
        .str("spin", opt.spin)
        .number("mu_s", opt.mu_s)
        .integer("fit_mu_s", opt.fit_mu_s)
        .number("fixed_floor_GeV", opt.fixed_floor)
        .str("kinetic", opt.kinetic)
        .str("fitter", opt.fitter)
        .integer("multistart", opt.multistart)
        .integer("basis_n", opt.basis.n)
        .number("basis_r1", opt.basis.r1)
        .number("basis_rmax", opt.basis.rmax)
        .number("lambda_qcd", opt.lambda_qcd)
        .number("gluon_mass", opt.gluon_mass)
        .integer("seed", static_cast<long>(opt.seed));
    json_all.put("config", cfg.render(2));
    json_all.put("systems", json_systems);
    f << json_all.render() << "\n";
    out << "wrote " << opt.json_out << "\n";
  }
  return 0;
}

}  // namespace qgem
