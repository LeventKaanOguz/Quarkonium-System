#pragma once
// Gauss quadrature by Golub-Welsch (notes Alg. 11.2). Rules are cached by (n,
// alpha): rebuilding one per matrix element costs three orders of magnitude.
#include "core/types.hpp"

namespace qgem {

struct QuadRule {
  Vec x, w;
};

/// n-point Gauss-Legendre rule on [0, 1].
const QuadRule& gauss_legendre01(int n);

/// n-point generalized Gauss-Laguerre rule for the weight t^alpha e^{-t} on [0,
/// inf), alpha > -1.
const QuadRule& gauss_laguerre(int n, double alpha);

}  // namespace qgem
