#pragma once
// Named model parameters with fixing, tying and transforms (notes 10.3.5,
// Alg. 11.9: bounds by transform, never by clipping). theta = full physical
// vector; u = vector of the free, transformed parameters that the minimizer
// sees.
#include <limits>
#include <string>
#include <vector>

#include "core/types.hpp"

namespace qgem {

enum class Transform { Identity, Log };

struct ParamSpec {
  std::string name;
  double value = 0.0;
  bool free = true;
  Transform tr = Transform::Identity;
  double lo =
      -std::numeric_limits<double>::infinity();  ///< physical bounds (checked,
                                                 ///< not clipped)
  double hi = std::numeric_limits<double>::infinity();
  std::string tied_to = {};  ///< non-empty: theta[this] = theta[tied_to]; the
                             ///< parameter is not free itself
};

class ParamSet {
 public:
  ParamSet() = default;
  explicit ParamSet(std::vector<ParamSpec> specs);

  int n_total() const { return static_cast<int>(specs_.size()); }
  int n_free() const { return static_cast<int>(free_idx_.size()); }
  const std::vector<ParamSpec>& specs() const { return specs_; }
  int index(const std::string& name) const;  ///< throws if unknown
  std::vector<std::string> free_names() const;
  const std::vector<int>& free_indices() const { return free_idx_; }

  Vec theta0() const;  ///< starting physical vector
  Vec u0() const;      ///< starting free vector
  Vec theta(const Vec& u)
      const;  ///< u -> full theta (transforms, ties, fixed values)
  Vec u(const Vec& theta) const;      ///< inverse on the free parameters
  Mat dtheta_du(const Vec& u) const;  ///< n_total x n_free
  bool within_bounds(const Vec& theta) const;

  /// copy with a different set of free parameters (by name) and a new starting
  /// theta
  ParamSet with_values(const Vec& theta) const;
  ParamSet with_free(const std::vector<std::string>& free_names) const;

 private:
  std::vector<ParamSpec> specs_;
  std::vector<int> free_idx_;
  std::vector<int> tie_src_;  // -1 if not tied
};

}  // namespace qgem
