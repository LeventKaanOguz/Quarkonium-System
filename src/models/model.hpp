#pragma once
// Abstract bound-state model: it assembles H, S (and an optional first-order
// matrix W) for one channel. Everything downstream (solver, fit, derivatives,
// uncertainties) is shared between models.
#include <array>
#include <string>
#include <vector>

#include "core/types.hpp"
#include "models/cache.hpp"
#include "models/dual.hpp"
#include "models/params.hpp"
#include "spin/spin.hpp"

namespace qgem {

struct ChannelMatrices {
  Mat H, S;
  Mat W;  ///< first-order spin matrix: M = offset + E + c^T W c (empty if
  Mat X;  ///< optional orthonormalizer (X^T S X = 1) the model built H in
          ///< (closure); empty: the solver computes it unused)
  std::vector<int> L;           ///< partial wave of each block
  std::vector<int> block_size;  ///< basis functions per block
  int n_dp =
      0;  ///< number of parameter derivatives stored below (0: not requested)
  std::array<Mat, kMaxP> dH, dW;  ///< dH/dtheta_a, dW/dtheta_a (S does not
                                  ///< depend on physics parameters)
  bool has_W() const { return W.size() > 0; }
};

class Model {
 public:
  virtual ~Model() = default;
  virtual std::string name() const = 0;
  virtual std::vector<std::string> param_names() const = 0;
  int n_params() const { return static_cast<int>(param_names().size()); }
  /// Channels the model solves. Models that treat a coupling to first order
  /// list the pieces separately.
  virtual std::vector<Channel> channels(int Lmax, bool equal_masses) const = 0;
  virtual double mass_offset(
      const Vec& theta) const = 0;  ///< m1 + m2 (A, B) or 0 (C, D)
  virtual Vec mass_offset_gradient(const Vec& theta) const = 0;
  /// Hint from the fitter: only these parameters will be differentiated (models
  /// that differentiate numerically use it to skip the rest). Default: ignored.
  virtual void restrict_derivatives(const std::vector<int>& /*active*/) const {}
  virtual ChannelMatrices build(const Channel& ch, const Vec& theta,
                                MatrixCache& cache, bool derivatives) const = 0;
};

}  // namespace qgem
