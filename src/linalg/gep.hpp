#pragma once
// Generalized symmetric eigenproblem H c = E S c with a nearly singular overlap
// S.
#include "core/types.hpp"

namespace qgem {

struct GepResult {
  Vec E;           ///< eigenvalues, ascending
  Mat C;           ///< eigenvectors as columns, c_k^T S c_k = 1
  Mat X;           ///< orthonormalizer, X^T S X = 1 (columns = kept directions)
  int n_kept = 0;  ///< number of directions kept
  double s_min = 0, s_max = 0;  ///< smallest kept / largest eigenvalue of S
};

/// X with X^T S X = 1, keeping eigen-directions of S with s_i > eta * s_max.
Mat orthonormalizer(const Mat& S, double eta);

/// Canonical orthogonalization (default). Throws if no direction survives the
/// threshold.
GepResult solve_gep(const Mat& H, const Mat& S, double eta = 1e-12);

/// Cholesky reduction S = L L^T; only for cross-checks (fails for
/// ill-conditioned S).
GepResult solve_gep_cholesky(const Mat& H, const Mat& S);

}  // namespace qgem
