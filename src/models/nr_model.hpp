#pragma once
// Models A and B (notes 4.1, 4.2): non-relativistic kinetic energy, Coulomb +
// linear (A) or screened (B) potential, Breit-Fermi spin terms. Parameters
// theta = (m1, m2, alpha_s, b, c, sigma, mu_s).
#include "models/model.hpp"

namespace qgem {

enum class SpinMode {
  None,        ///< spin-independent only (no contact term): reference tests
  FirstOrder,  ///< smeared contact term in H; L.S and tensor added at first
               ///< order with the eigenvector (point-like 1/r^3, finite for L
               ///< >= 1). The worked example of notes 10.10.
  Smeared      ///< everything in H, built from the Gaussian-smeared Coulomb
           ///< potential (regular, bounded below): coupled channels (3S1-3D1,
           ///< ...) and, for m1 != m2, 1L_L-3L_L mixing.
};

/// Strong-coupling form. Constant: V = -4 alpha_s/(3r), spin terms built from
/// the Gaussian-smeared Coulomb potential (width sigma). Gaussian: the
/// Godfrey-Isgur running coupling alpha_s(Q) = alpha_s sum_k w_k exp(-Q^2/4
/// gamma_k^2) (weights w_k = alpha_k/0.6, so that alpha_s is the saturation
/// value alpha_s(0)); every potential is then an exact sum of erf(gamma_k r)/r
/// terms, smeared by sigma in the spin terms (widths compose, notes 6.7).
enum class Running { Constant, Gaussian };

struct NRModelOptions {
  Running running = Running::Constant;
  bool screened = false;  ///< Model B (adds mu_s)
  SpinMode spin = SpinMode::FirstOrder;
  bool contact_in_all_L =
      true;  ///< smeared contact term in every L (else only L = 0)
  /// Running::Gaussian with a shape other than Godfrey-Isgur's: weights
  /// (sum = 1) and ranges of the Gaussian sum (see fit_gaussian_shape), and a
  /// label for name(). Empty: the Godfrey-Isgur shape.
  std::vector<double> shape_w = {}, shape_gamma = {};
  std::string shape_label = {};
  /// semi-relativistic kinetic energy sum_i (sqrt(p^2 + m_i^2) - m_i) instead
  /// of p^2/2 m_r (notes 4.1, variants): all orders of the mass-velocity
  /// term. Requires a smeared (running) Coulomb potential: with a point
  /// Coulomb -kappa/r the operator is bounded below only for kappa <= 2/pi
  /// (Herbst).
  bool relativistic_kinetic = false;
};

enum NRParam { kM1 = 0, kM2, kAlphaS, kB, kC, kSigma, kMuS, kNRParams };

class NRModel final : public Model {
 public:
  explicit NRModel(NRModelOptions opt = {}) : opt_(opt) {}
  const NRModelOptions& options() const { return opt_; }

  std::string name() const override;
  std::vector<std::string> param_names() const override;
  std::vector<Channel> channels(int Lmax, bool equal_masses) const override;
  double mass_offset(const Vec& theta) const override {
    return theta(kM1) + theta(kM2);
  }
  Vec mass_offset_gradient(const Vec& theta) const override;
  ChannelMatrices build(const Channel& ch, const Vec& theta, MatrixCache& cache,
                        bool derivatives) const override;

  /// Convenience: default parameter specs (log transform on the positive ones).
  ParamSet default_params(double m1, double m2, double alpha_s, double b,
                          double c, double sigma, double mu_s = 0.1) const;

 private:
  NRModelOptions opt_;
};

}  // namespace qgem
