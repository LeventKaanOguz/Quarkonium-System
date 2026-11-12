#pragma once
// chi^2 of a model against a dataset: whitened residuals, Feynman-Hellmann
// Jacobian, level tracking, theory floor (notes 10.2, Alg. 11.7).
#include <map>
#include <memory>

#include "data/dataset.hpp"
#include "fit/least_squares.hpp"
#include "models/model.hpp"

namespace qgem {

/// A datum names a channel ("3P0") or one partial wave of a coupled channel
/// ("3S1" of 3S1-3D1). Returns the channel to solve and, for the latter, the
/// component that labels the state (-1: the channel itself).
std::pair<Channel, int> resolve_channel(const std::string& label,
                                        const Model& model, const Vec& theta,
                                        int Lmax);
/// Index of the n-th level dominated by component `comp` (any level if comp <
/// 0); levels.size() if unavailable.
int pick_level(const SolvedChannel& sc, int comp, int n);

class Chi2Problem final : public LeastSquaresProblem {
 public:
  /// `Lmax` is the highest partial wave solved (channels needed by the data are
  /// found automatically).
  Chi2Problem(const Model& model, Dataset data, ParamSet params,
              MatrixCache cache, double sigma_th = 0.0);

  int n_data() const override {
    return n_fit_ + static_cast<int>(priors_.size());
  }
  int n_free() const override { return params_.n_free(); }
  Result<Vec> residuals(const Vec& u) override;
  Result<Mat> jacobian(const Vec& u) override;

  /// Gaussian prior on a physical parameter (notes 10.9.4: penalties = priors):
  /// adds ((theta_a - mean)/sigma)^2 to chi^2 as an extra whitened residual. A
  /// prior is the Bayesian form of a penalty term.
  void add_prior(int theta_index, double mean, double sigma);
  /// relative floor of each fitted state (Datum::floor_w), fit order
  Vec floor_weights() const;
  int n_measurements() const {
    return n_fit_;
  }  ///< fitted data only (without priors)

  void set_sigma_th(double s);  ///< GeV, common theory floor added in
                                ///< quadrature to every fitted datum
  double sigma_th() const { return sigma_th_; }
  const ParamSet& params() const { return params_; }
  const Dataset& data() const { return data_; }
  Model const& model() const { return model_; }
  MatrixCache& cache() { return cache_; }

  /// Model predictions of all data (fitted and predicted), GeV, at the physical
  /// parameters theta.
  Result<Vec> predictions(const Vec& theta);
  /// Predictions of all data and their gradients dM/dtheta (rows: data,
  /// columns: the model's physical parameters).
  Result<std::pair<Vec, Mat>> predictions_with_gradients(const Vec& theta);
  /// y - mu for the fitted data at u.
  Result<Vec> raw_residuals(const Vec& u);
  /// The data covariance used for whitening: V_exp + sigma_th^2 I.
  Mat covariance() const;
  Vec exp_errors() const;      ///< sqrt(diag V_exp) of the fitted data
  Mat exp_covariance() const;  ///< V_exp of the fitted data (PDG only)
  /// Replace the measured values (parametric bootstrap); the covariance is
  /// unchanged.
  void set_measurements(const Vec& y_fitted);
  Vec measurements() const;

 private:
  Result<Vec> mu_fitted(const Vec& theta, Mat* dmu_dtheta);
  void refresh_whitener();

  const Model& model_;
  Dataset data_;
  ParamSet params_;
  MatrixCache cache_;
  double sigma_th_;
  int n_fit_ = 0;
  Mat Linv_;  // whitener: z = Linv (y - mu)
  std::map<std::string, Vec>
      ref_vec_;  // level-tracking reference per (channel, n)
  struct Prior {
    int index;
    double mean, sigma;
  };
  std::vector<Prior> priors_;
  bool refs_set_ = false;
};

}  // namespace qgem
