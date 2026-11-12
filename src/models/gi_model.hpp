#pragma once
// Models C and D (notes 4.3, 4.4): the Godfrey-Isgur relativized Hamiltonian
// and its screened ("modified") variant.
//   H = E1(p) + E2(p) + (1 + p^2/E1E2)^{1/2} G~(r) (1 + p^2/E1E2)^{1/2} + S~(r)
//   + H~c + H~t + H~so(v) + H~so(s)
// The eigenvalue is the meson mass itself (no offset). Non-local sandwiches
// f(p) W(r) f(p) use the closure F' W' F' in the orthonormalized basis
// (notes 7.7(ii)); the two-step method is available for validation.
#include <map>
#include <tuple>

#include "models/model.hpp"

namespace qgem {

enum GIParam {
  kGm1 = 0,
  kGm2,
  kGb,
  kGc,
  kGsigma0,
  kGs,
  kGepsC,
  kGepsT,
  kGepsSov,
  kGepsSos,
  kGalpha1,
  kGalpha2,
  kGalpha3,
  kGgamma1,
  kGgamma2,
  kGgamma3,
  kGmuS,
  kGIParams
};

struct GIModelOptions {
  bool screened =
      false;    ///< Model D: S~ -> S~_scr with the smeared exponential
  int nq = 64;  ///< Gauss-Laguerre nodes of the momentum integrals
};

/// How alpha_s runs inside the GI Hamiltonian: the three fixed Godfrey-Isgur
/// Gaussians, the same shape with free amplitudes alpha_k, or a constant
/// coupling (one amplitude, gamma_1 -> infinity: erf(gamma r)/r -> 1/r
/// exactly).
enum class GiRunning { Fixed, Free, Constant };

class GIModel final : public Model {
 public:
  explicit GIModel(GIModelOptions opt = {}) : opt_(opt) {}
  std::string name() const override {
    return opt_.screened ? "D (modified GI)" : "C (Godfrey-Isgur)";
  }
  std::vector<std::string> param_names() const override;
  std::vector<Channel> channels(int Lmax, bool equal_masses) const override;
  double mass_offset(const Vec&) const override { return 0.0; }
  Vec mass_offset_gradient(const Vec&) const override {
    return Vec::Zero(kGIParams);
  }
  void restrict_derivatives(const std::vector<int>& active) const override {
    active_ = active;
  }
  ChannelMatrices build(const Channel& ch, const Vec& theta, MatrixCache& cache,
                        bool derivatives) const override;

  /// The parameters of the original paper (notes 12.3, "recalled from memory --
  /// verify before use").
  /// `with_c`: keep the additive constant c free (start -0.253, as in the
  /// paper) or fix c = 0.
  ParamSet godfrey_isgur_params(double m1, double m2,
                                GiRunning running = GiRunning::Fixed,
                                bool with_c = true) const;
  /// Smearing width sigma_12 (6.8).
  static double sigma12(const Vec& theta);

 private:
  /// Momentum-space matrices of one basis, keyed by (kind, L, m_a, m_b, eps).
  /// They depend on the masses and the eps only, so one memo serves every
  /// finite-difference evaluation of build() that leaves those unchanged.
  using MomentumMemo =
      std::map<std::tuple<int, int, double, double, double>, Mat>;
  ChannelMatrices assemble(const Channel& ch, const Vec& theta,
                           MatrixCache& cache, MomentumMemo& memo) const;
  GIModelOptions opt_;
  mutable std::vector<int> active_;  // empty: differentiate everything
};

}  // namespace qgem
