#pragma once
// Annihilation widths as data: their predictions at LO and NLO, the standard
// theory error, the decay benchmark (a Pearson chi^2 of the measured widths
// against a fit that never saw them) and a least-squares problem for fits
// that include them.
#include <string>
#include <vector>

#include "data/dataset.hpp"
#include "fit/chi2_problem.hpp"
#include "models/cache.hpp"
#include "models/model.hpp"

namespace qgem {

/// An annihilation width of a datum: e+e- for 3S1/3D1 states, gamma gamma
/// for 1S0, 3P0 and 3P2. `key` names the one-loop QCD factor (qcd_factor).
struct WidthChannel {
  int datum = 0;  ///< index into Dataset::d
  std::string process, key;
};
std::vector<WidthChannel> width_channels(const Dataset& ds);

/// Widths in GeV: LO (Van Royen-Weisskopf), NLO (x the one-loop factor with
/// the MS-bar alpha_s at mu = M/2) and NLO at mu = M/4 and M; <v^2> of the
/// state. NaN if the state cannot be solved.
struct WidthValue {
  double lo, nlo, nlo_lo, nlo_hi, v2, M;
};
std::vector<WidthValue> compute_widths(const Model& model, const Vec& theta,
                                       const Dataset& ds, MatrixCache& cache,
                                       double eQ,
                                       const std::vector<WidthChannel>& ch);

/// Standard theory error of an NLO width: the truncation of the alpha_s
/// series, half the scale band (mu = M/4 ... M), in quadrature with the
/// O(v^2) relativistic correction of NRQCD power counting, <v^2> Gamma_NLO.
double width_theory_error(const WidthValue& w);
/// Theory covariance: the scale part fully correlated within a process (one
/// mu moves every width of the process together), the v^2 part diagonal.
Mat width_theory_covariance(const std::vector<WidthValue>& w,
                            const std::vector<WidthChannel>& ch);

/// Measured widths: PDG value and error (GeV) of the channels that have one.
struct WidthData {
  std::vector<WidthChannel> ch;
  Vec y, s_exp;
};

/// Least-squares problem of the measured widths: z = L^-1 (y - Gamma_NLO(u)),
/// V = diag(s_exp^2) + the theory covariance at u_ref, fixed during the fit
/// (the errors do not move with the parameters). Jacobian by central
/// differences in u.
class WidthProblem final : public LeastSquaresProblem {
 public:
  WidthProblem(const Model& model, const Dataset& ds, ParamSet params,
               MatrixCache cache, double eQ, WidthData data, const Vec& u_ref);
  int n_data() const override { return static_cast<int>(data_.y.size()); }
  int n_free() const override { return params_.n_free(); }
  Result<Vec> residuals(const Vec& u) override;
  Result<Mat> jacobian(const Vec& u) override;
  /// NLO widths of the measured channels at u (GeV).
  Result<Vec> predictions(const Vec& u);
  const Mat& covariance() const { return V_; }
  const WidthData& data() const { return data_; }

 private:
  const Model& model_;
  const Dataset& ds_;
  ParamSet params_;
  MatrixCache cache_;
  double eQ_;
  WidthData data_;
  Mat V_, Linv_;
};

/// Masses and widths in one least-squares problem (same parameters): the
/// mass residuals of `m` stacked on the width residuals of `w`.
class MassWidthProblem final : public LeastSquaresProblem {
 public:
  MassWidthProblem(Chi2Problem& m, WidthProblem& w) : m_(m), w_(w) {}
  int n_data() const override { return m_.n_data() + w_.n_data(); }
  int n_free() const override { return m_.n_free(); }
  Result<Vec> residuals(const Vec& u) override;
  Result<Mat> jacobian(const Vec& u) override;
  Chi2Problem& masses() { return m_; }
  WidthProblem& widths() { return w_; }

 private:
  Chi2Problem& m_;
  WidthProblem& w_;
};

/// The combined fit with the mass floor iterated on the mass block: chi2_m =
/// N_m - sum_{i in masses} h_ii (leverage share of the combined fit); the
/// width errors stay fixed, so the width block of chi^2 remains a test.
struct MassWidthFit {
  FitResult fit;
  double sigma_th = 0, nu_m = 0, chi2_m = 0, chi2_w = 0;
};
MassWidthFit fit_masses_and_widths(MassWidthProblem& p, Vec u, double sigma0,
                                   const FitOptions& opt = {});

}  // namespace qgem
