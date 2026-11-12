#pragma once
// Solving a channel: eigenlevels, masses, level identification and
// Feynman-Hellmann gradients.
#include <span>
#include <string>
#include <vector>

#include "core/types.hpp"
#include "linalg/gep.hpp"
#include "models/model.hpp"

namespace qgem {

struct StateId {
  std::string channel;  ///< Channel::label(), e.g. "3P0"
  int n = 0;            ///< radial index within the channel, 0 = lowest
  bool operator==(const StateId&) const = default;
};

struct Level {
  StateId id;
  double E =
      0.0;  ///< eigenvalue of H (without offset, without first-order shift)
  double shift = 0.0;        ///< first-order c^T W c
  double M = 0.0;            ///< offset + E + shift
  Vec c;                     ///< S-normalized eigenvector (all blocks stacked)
  std::vector<double> prob;  ///< probability of each component
};

struct SolvedChannel {
  Channel channel;
  ChannelMatrices mats;
  GepResult gep;
  std::vector<Level> levels;
  const Level& level(int n) const { return levels.at(n); }
  /// dM_k/dtheta (length n_params) by Feynman-Hellmann, incl. the first-order
  /// shift (eq. 7.34, 7.35). Requires `mats` built with derivatives.
  Vec mass_gradient(int k, const Model& model, const Vec& theta) const;
};

/// Build and diagonalize one channel. `derivatives` stores dH/dW for the
/// gradients.
SolvedChannel solve_channel(const Model& model, const Channel& ch,
                            const Vec& theta, MatrixCache& cache,
                            bool derivatives = false);

/// Index of the level in `now` that continues level a of `before` (maximum
/// |c_b^T S c_a|); equal bases.
int track_level(const Mat& S, const Vec& c_before, const SolvedChannel& now);

}  // namespace qgem
