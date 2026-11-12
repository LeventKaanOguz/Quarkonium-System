#pragma once
// The end-to-end run: choose a model and a running form of alpha_s (with or
// without the additive constant c), fit each system to every PDG state below
// threshold, and print every state found in the PDG file with its model value,
// uncertainties, pulls, numerical error and leptonic / two-photon widths. B_c
// is predicted from the two fits.
#include <functional>
#include <iosfwd>
#include <map>
#include <string>
#include <vector>

#include "basis/basis.hpp"
#include "fit/chi2_problem.hpp"

namespace qgem {

struct RunOptions {
  std::string model = "A";  ///< A, B (non-relativistic Cornell / screened), C,
                            ///< D (Godfrey-Isgur / screened)
  std::string running =
      "const";        ///< const | gi | free   (alpha_s: constant, Godfrey-Isgur
                      ///< running, running with free amplitudes)
  double mu_s = 0.1;  ///< screening mass of models B and D, GeV
  bool fit_mu_s = false;  ///< fit mu_s instead of keeping it fixed (the charm
                          ///< data alone drive it to zero)
  bool with_c = false;    ///< free additive constant c (nonphysical, degenerate
                          ///< with the quark mass)
  std::string spin = "first";  ///< models A/B: none | first | smeared
  std::vector<std::string> systems = {"ccbar", "bbar", "bc"};
  std::string root = ".";  ///< repository root (config/, data/)
  std::string pdg_file;    ///< default: <root>/data/pdg-2026.0.sqlite
  // (r1, rmax): in the flat variational optimum of the level sums of every
  // model; n = 37: the smallest member of the nested family 10, 19, 37, ...
  // whose states move by < 0.1 MeV to the next (report, numerical section)
  BasisSpec basis{.n = 37, .r1 = 0.05, .rmax = 15.0, .eta = 1e-12};
  int bootstrap = 0, mcmc = 0, threads = 1;
  unsigned long seed = 20261011;
  std::string json_out;  ///< optional machine-readable result

  /// Shape parameters of the frozen / APT couplings (models A/B with
  /// --running frozen|apt); the saturation value alpha_s(0) is fitted.
  double lambda_qcd = 0.30, gluon_mass = 0.50;
  int nf = 4;
  /// Theory floor: < 0 iterates it to chi^2/nu = 1 (the default); >= 0 keeps
  /// it fixed (GeV), which makes chi^2, AIC and BIC comparable across models.
  double fixed_floor = -1.0;
  int multistart = 0;  ///< extra random starts; the best fit is kept
  /// Parameters removed from the fit and fixed at a value (--fix NAME=VALUE).
  std::vector<std::pair<std::string, double>> fixed;
  /// Extra starting points per system (--start SYSTEM:NAME=V,...), e.g. the
  /// minimum found by another minimiser or by the nested model. Only free
  /// parameters are set; names the model lacks are ignored.
  std::map<std::string,
           std::vector<std::vector<std::pair<std::string, double>>>>
      starts;
  /// Parameters released into the fit although the model fixes them
  /// (--free NAME), e.g. the GI smearing slope s.
  std::vector<std::string> freed;
  /// state-dependent theory floor: system -> (state label -> relative floor
  /// w), floor of state i = w_i sigma_th (--floor-weights); unlisted: 1
  std::map<std::string, std::map<std::string, double>> floor_weights;
  /// A/B kinetic energy: "nr" (p^2/2 m_r) or "rel" (sqrt(p^2+m^2) - m)
  std::string kinetic = "nr";
  /// parameters shared by a joint c c-bar + b b-bar fit after the separate
  /// ones (--shared b): flavour independence of the potential; empty: none
  std::vector<std::string> shared;
  /// print the expansion coefficients of every wave function (always in JSON)
  bool wavefunctions = false;
  bool exact_loo = false;  ///< exact leave-one-out refits (cross-validation)
  /// leave-multiplet-out: refit without every state of one (n, L) multiplet
  /// and predict them (cross-validation against overfitting the splittings)
  bool lmo = false;
  /// repeated K-fold cross-validation (0 = off): random folds, 10 repeats
  int kfold = 0;
  /// learning curve: train and test RMS against the training-set size
  bool learning_curve = false;
  /// optimise (r1, rmax) per system variationally (hyper/basis_opt.hpp)
  bool optimize_basis = false;
  /// > 0: nested-convergence table of every datum at the fit, sizes
  /// 10, 19, 37, 73, 145, with this target (GeV) for the choice of n
  double nested_target = 0;
  int gi_nq = 64;  ///< Gauss-Laguerre nodes of the GI momentum integrals
  std::string fitter = "lm";  ///< label of the minimiser in `fit.engine`
  FitOptions fit;             ///< fit.engine empty: Levenberg-Marquardt
  /// Optional cross-check at the best fit (the ROOT/Minuit2 Hesse + Minos
  /// analysis); writes text to the stream and returns a rendered JSON object.
  std::function<std::string(Chi2Problem&, const FitResult&, const ParamSet&,
                            std::ostream&)>
      cross_check;
};

/// Returns 0 on success. Everything is written to `out`.
int run(const RunOptions& opt, std::ostream& out);

}  // namespace qgem
