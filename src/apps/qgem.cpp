// qgem: command-line front end. With no arguments it runs the default
// configuration (model A, constant alpha_s, c = 0) and prints every PDG state.
// See `qgem --help` and README.md.
#include <cstdio>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include "report/json.hpp"
#include "report/run.hpp"
#ifdef QGEM_HAVE_ROOT
#include "root/minuit2_engine.hpp"
#endif

namespace {

constexpr const char* kUsage = R"(usage: qgem [options]

model
  --model A|B|C|D        A non-relativistic Cornell + Breit-Fermi   B A, screened
                         C Godfrey-Isgur relativized               D C, screened
  --running FORM         alpha_s: A/B const|gi|frozen|apt, C/D const|gi|free
  --c on|off             additive constant c free (default off: c = 0)
  --spin none|first|smeared   spin terms of A/B (default first)
  --mu-s GeV             screening mass of B/D (default 0.1)
  --fit-mu-s             fit mu_s instead of fixing it
  --kinetic nr|rel       A/B kinetic energy: p^2/2m_r, or sqrt(p^2+m^2) - m
                         (needs a running, i.e. smeared, coupling)
  --fix NAME=VALUE       remove a parameter from the fit (repeatable)
  --free NAME            add a parameter the model keeps fixed (repeatable)
  --lambda GeV --gluon-mass GeV   shape of the frozen/APT couplings (0.30, 0.50)
fit
  --systems ccbar,bbar,bc   (default all three; B_c is predicted)
  --floor MeV            fixed theory floor (default: iterate to chi2/nu = 1)
  --floor-weights SYS:LABEL=W,...   state-dependent floor W_i sigma_th
  --wf-band              68 % parametric band of every wave function (JSON)
  --fit-widths           also fit the measured annihilation widths of the
                         fitted states (default: masses only, the widths a
                         held-out benchmark)
  --wavefunctions        print the expansion coefficients of every wave
                         function (the JSON always has them)
  --shared NAMES         joint c c-bar + b b-bar fit with these parameters
                         common (e.g. b, or b,alpha_s): the universality test,
                         and the B_c from the joint parameters
  --fitter lm|minuit2    minimiser (default: minuit2 in the ROOT build, else lm;
                         lm is faster for models C/D but more often stops in
                         a local minimum)
  --minuit-check         Minuit2 Migrad+Hesse+Minos cross-check at the best fit
  --multistart K         K extra random starts, keep the best
  --start SYS:N=V,...    one more start for system SYS (repeatable)
  --loo                  exact leave-one-out refits (+ jackknife errors)
  --lmo                  leave-multiplet-out refits (all states of one n, L)
  --kfold K              repeated K-fold cross-validation (10 random splits)
  --learning-curve       train/test RMS against the number of fitted states
  --nested-check MeV     convergence of every state on nested bases
                         n = 10, 19, 37, 73, 145; smallest n below MeV
  --eta X                overlap-eigenvalue cut of the orthogonalisation
  --gi-nq N              quadrature nodes of the GI momentum integrals (64)
  --optimize-basis       variational (r1, rmax) per system, alternated with
                         the fit (minimises the sum of the levels used)
  --bootstrap B          parametric bootstrap replicas
  --mcmc STEPS           Metropolis steps per chain (4 chains)
numerics / io
  --basis n,r1,rmax      Gaussian basis (default 36,0.02,60)
  --threads T  --seed S  --out result.json  --pdg FILE  --root DIR
)";

std::vector<std::string> split(const std::string& s, char sep) {
  std::vector<std::string> out;
  std::stringstream ss(s);
  std::string t;
  while (std::getline(ss, t, sep)) out.push_back(t);
  return out;
}

#ifdef QGEM_HAVE_ROOT
// Minuit2 at the best fit: Migrad restarted from the LM minimum, Hesse
// covariance and Minos intervals, compared with the linear (J^T J)^-1 errors.
std::string minuit_cross_check(qgem::Chi2Problem& prob,
                               const qgem::FitResult& fit,
                               const qgem::ParamSet& ps, std::ostream& out) {
  using namespace qgem;
  const MinuitErrors e = minuit2_errors(prob, fit.u, true);
  const Uncertainty U = analyse(fit, ps);
  const auto names = ps.free_names();
  const int p = static_cast<int>(names.size());
  const Vec th = ps.theta(e.u);
  Vec dchi2(1);
  dchi2(0) = (prob.residuals(e.u)->squaredNorm() - fit.chi2);
  out << "\nMinuit2 cross-check (Migrad " << (e.valid ? "valid" : "INVALID")
      << ", EDM " << std::scientific << std::setprecision(1) << e.edm
      << std::defaultfloat << ", " << e.nfcn
      << " calls; chi2(Minuit) - chi2(LM) = " << std::setprecision(3)
      << dchi2(0) << "):\n"
      << "  param      Migrad value      linear sigma   Minos -/+ (physical)\n";
  Vec lo(p), hi(p), val(p);
  for (int k = 0; k < p; ++k) {
    const int i = ps.free_indices()[k];
    Vec ul = e.u, uh = e.u;
    ul(k) += e.minos_lo(k);
    uh(k) += e.minos_hi(k);
    val(k) = th(i);
    lo(k) = ps.theta(ul)(i) - th(i);
    hi(k) = ps.theta(uh)(i) - th(i);
    if (lo(k) > hi(k)) std::swap(lo(k), hi(k));  // decreasing transform
    out << "  " << std::left << std::setw(9) << names[k] << std::right
        << std::fixed << std::setprecision(6) << std::setw(14) << th(i)
        << std::setw(14) << U.sigma(k) << "   " << std::setprecision(5) << lo(k)
        << " / +" << hi(k) << (e.minos_valid[k] ? "" : "  (Minos invalid)")
        << "\n"
        << std::defaultfloat;
  }
  json::Object j;
  j.integer("valid", e.valid)
      .number("edm", e.edm)
      .integer("nfcn", e.nfcn)
      .number("dchi2_vs_lm", dchi2(0))
      .put("value", json::array(val))
      .put("minos_lo", json::array(lo))
      .put("minos_hi", json::array(hi))
      .put("hesse_cov_u", json::array(e.cov_u))
      .put("linear_cov_u", json::array(U.C));
  return j.render(6);
}
#endif

// The minimiser: Levenberg-Marquardt (built in) or Minuit2 Migrad (ROOT).
void set_fitter(qgem::RunOptions& o, const std::string& name) {
  if (name == "lm") {
    o.fit.engine = nullptr;
  } else if (name == "minuit2") {
#ifdef QGEM_HAVE_ROOT
    o.fit.engine = [](qgem::LeastSquaresProblem& p, const qgem::Vec& u,
                      const qgem::FitOptions& f) {
      return qgem::minuit2_minimize(p, u, f);
    };
#else
    throw std::runtime_error(
        "--fitter minuit2 needs the ROOT build (cmake -DQGEM_WITH_ROOT=ON)");
#endif
  } else {
    throw std::runtime_error("--fitter lm|minuit2");
  }
  o.fitter = name;
}

}  // namespace

int main(int argc, char** argv) {
  try {
    qgem::RunOptions o;
#ifdef QGEM_SOURCE_DIR
    o.root = QGEM_SOURCE_DIR;
#endif
#ifdef QGEM_HAVE_ROOT
    // Minuit2 by default: in the model sweep, Levenberg-Marquardt with nine
    // starts stopped in a local minimum in 33 of 134 fits, one Minuit2 start
    // in 1 (docs/report/REPORT.md, section 2.5). LM is faster: --fitter lm.
    set_fitter(o, "minuit2");
#endif
    for (int i = 1; i < argc; ++i) {
      const std::string k = argv[i];
      auto val = [&]() -> std::string {
        if (i + 1 >= argc) throw std::runtime_error("missing value for " + k);
        return argv[++i];
      };
      if (k == "--model")
        o.model = val();
      else if (k == "--running")
        o.running = val();
      else if (k == "--c") {
        const std::string v = val();
        if (v != "on" && v != "off") throw std::runtime_error("--c on|off");
        o.with_c = v == "on";
      } else if (k == "--spin")
        o.spin = val();
      else if (k == "--mu-s")
        o.mu_s = std::stod(val());
      else if (k == "--fit-mu-s")
        o.fit_mu_s = true;
      else if (k == "--fix") {
        const std::string v = val();
        const auto eq = v.find('=');
        if (eq == std::string::npos)
          throw std::runtime_error("--fix NAME=VALUE");
        o.fixed.emplace_back(v.substr(0, eq), std::stod(v.substr(eq + 1)));
      } else if (k == "--start") {
        // SYSTEM:NAME=V,NAME=V (repeatable: each is one more starting point)
        const std::string v = val();
        const auto colon = v.find(':');
        if (colon == std::string::npos)
          throw std::runtime_error("--start SYSTEM:NAME=VALUE,...");
        std::vector<std::pair<std::string, double>> seed;
        for (const auto& kv : split(v.substr(colon + 1), ',')) {
          const auto eq = kv.find('=');
          if (eq == std::string::npos)
            throw std::runtime_error("--start SYSTEM:NAME=VALUE,...");
          seed.emplace_back(kv.substr(0, eq), std::stod(kv.substr(eq + 1)));
        }
        o.starts[v.substr(0, colon)].push_back(std::move(seed));
      } else if (k == "--lambda")
        o.lambda_qcd = std::stod(val());
      else if (k == "--gluon-mass")
        o.gluon_mass = std::stod(val());
      else if (k == "--systems")
        o.systems = split(val(), ',');
      else if (k == "--shared")
        o.shared = split(val(), ',');
      else if (k == "--wavefunctions")
        o.wavefunctions = true;
      else if (k == "--fit-widths")
        o.fit_widths = true;
      else if (k == "--wf-band")
        o.wf_band = true;
      else if (k == "--basis") {
        const auto v = split(val(), ',');
        if (v.size() != 3) throw std::runtime_error("--basis n,r1,rmax");
        o.basis.n = std::stoi(v[0]);
        o.basis.r1 = std::stod(v[1]);
        o.basis.rmax = std::stod(v[2]);
      } else if (k == "--floor")
        o.fixed_floor = 1e-3 * std::stod(val());
      else if (k == "--fitter")
        set_fitter(o, val());
      else if (k == "--minuit-check") {
#ifdef QGEM_HAVE_ROOT
        o.cross_check = minuit_cross_check;
#else
        throw std::runtime_error("--minuit-check needs the ROOT build");
#endif
      } else if (k == "--multistart")
        o.multistart = std::stoi(val());
      else if (k == "--floor-weights") {
        // SYSTEM:LABEL=W,LABEL=W (relative floor of each listed state)
        const std::string v = val();
        const auto colon = v.find(':');
        if (colon == std::string::npos)
          throw std::runtime_error("--floor-weights SYSTEM:LABEL=W,...");
        auto& m = o.floor_weights[v.substr(0, colon)];
        for (const auto& kv : split(v.substr(colon + 1), ',')) {
          const auto eq = kv.rfind('=');
          const double w = std::stod(kv.substr(eq + 1));
          if (eq == std::string::npos || !(w > 0))
            throw std::runtime_error("--floor-weights: LABEL=W with W > 0");
          m[kv.substr(0, eq)] = w;
        }
      } else if (k == "--kinetic")
        o.kinetic = val();
      else if (k == "--free")
        o.freed.push_back(val());
      else if (k == "--nested-check")
        o.nested_target = 1e-3 * std::stod(val());
      else if (k == "--eta")
        o.basis.eta = std::stod(val());
      else if (k == "--gi-nq")
        o.gi_nq = std::stoi(val());
      else if (k == "--optimize-basis")
        o.optimize_basis = true;
      else if (k == "--kfold")
        o.kfold = std::stoi(val());
      else if (k == "--learning-curve")
        o.learning_curve = true;
      else if (k == "--lmo")
        o.lmo = true;
      else if (k == "--loo")
        o.exact_loo = true;
      else if (k == "--bootstrap")
        o.bootstrap = std::stoi(val());
      else if (k == "--mcmc")
        o.mcmc = std::stoi(val());
      else if (k == "--threads")
        o.threads = std::stoi(val());
      else if (k == "--seed")
        o.seed = std::stoul(val());
      else if (k == "--out")
        o.json_out = val();
      else if (k == "--pdg")
        o.pdg_file = val();
      else if (k == "--root")
        o.root = val();
      else if (k == "--help" || k == "-h") {
        std::fputs(kUsage, stdout);
        return 0;
      } else {
        throw std::runtime_error("unknown option '" + k + "' (try --help)");
      }
    }
    return qgem::run(o, std::cout);
  } catch (const std::exception& e) {
    std::fprintf(stderr, "error: %s\n", e.what());
    return 1;
  }
}
