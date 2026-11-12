#pragma once
// Momentum-space matrix elements (notes 7.6, 7.7). A Gaussian of range nu in r
// is a Gaussian of range 1/(4 nu) in p, so any f(p) is diagonal there and its
// matrix elements are one-dimensional integrals done by generalized
// Gauss-Laguerre; sandwiches f(p) W(r) f(p) are formed by closure or by the
// two-step method.
#include <functional>

#include "basis/basis.hpp"
#include "core/types.hpp"

namespace qgem {

/// <R_i| f(p) |R_j> = S_ij / Gamma(L+3/2) sum_k w_k f( sqrt(x_k / beta_ij) ),
/// beta_ij = (nu_i+nu_j)/(4 nu_i nu_j), nodes/weights of the generalized
/// Gauss-Laguerre rule with alpha = L + 1/2 (notes 11.10). Same L for bra and
/// ket.
Mat momentum_diagonal(const Basis& b, const std::function<double(double)>& f,
                      int nq = 64);

/// Closure in the orthonormalized basis X (X^T S X = 1): <f W g> ~ (X^T F
/// X)(X^T W X)(X^T G X) with F, G the momentum-diagonal matrices and W the
/// local matrix, all in the original basis. Returns the primed result (an Nk x
/// Nk matrix in the orthonormal basis).
Mat sandwich_closure(const Mat& F, const Mat& W, const Mat& G, const Mat& Xbra,
                     const Mat& Xket);

/// Two-step ("exact") method: tabulate chi_n(r) = [f(p) R_n](r) = sqrt(2/pi)
/// Int p^2 j_l(p r) f(p) R~_n(p) dp on radial quadrature nodes, then <f W
/// f>_{n'n} = Int r^2 chi_n' W chi_n dr. Bra and ket may have different L.
/// Costs a Hankel transform per basis function; use it for validation and for
/// the final numbers.
Mat sandwich_two_step(const Basis& bra, const Basis& ket,
                      const std::function<double(double)>& f,
                      const std::function<double(double)>& W,
                      double rmax = 40.0);

}  // namespace qgem
