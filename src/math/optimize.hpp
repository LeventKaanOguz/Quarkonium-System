#pragma once
// Derivative-free minimization (Nelder-Mead) for the few continuous basis
// hyperparameters (notes Alg. 11.8).
#include <functional>

#include "core/types.hpp"

namespace qgem {

struct NelderMeadOptions {
  int max_iter = 400;
  double xtol = 1e-6;   ///< simplex size
  double ftol = 1e-12;  ///< spread of the function values
  double step = 0.3;    ///< initial simplex edge
};
struct NelderMeadResult {
  Vec x;
  double f = 0;
  int iterations = 0;
  bool converged = false;
};

NelderMeadResult nelder_mead(const std::function<double(const Vec&)>& f,
                             const Vec& x0, const NelderMeadOptions& opt = {});

}  // namespace qgem
