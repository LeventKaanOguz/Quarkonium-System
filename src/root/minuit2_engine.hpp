#pragma once
// Optional CERN ROOT bridge (built only when ROOT is found): Minuit2 as a
// second, independent minimiser of the same objective sum_i z_i^2 that
// Levenberg-Marquardt minimises, plus Hesse (second-derivative covariance)
// and Minos (profile-likelihood, asymmetric) intervals.
//
// Nothing outside src/root/ and src/apps/ depends on ROOT: the core library,
// the unit tests and CI build without it.
#include <vector>

#include "fit/least_squares.hpp"

namespace qgem {

struct MinuitOptions {
  int strategy = 1;         ///< 0 fast, 1 default, 2 careful
  double tolerance = 1e-4;  ///< Migrad stops at EDM < 2e-3 * tolerance * Up
  bool analytic_gradient =
      true;              ///< use 2 J^T z (Feynman-Hellmann); else numerical
  unsigned max_fcn = 0;  ///< 0: Minuit's default
};

/// Migrad on sum z^2 (Up = 1, i.e. Delta chi^2 = 1 errors) from u0. Usable
/// as FitOptions::engine. The returned z and J are evaluated at the minimum
/// (so every downstream analysis is the same as for LM).
FitResult minuit2_minimize(LeastSquaresProblem& problem, const Vec& u0,
                           const FitOptions& opt, const MinuitOptions& mo = {});

struct MinuitErrors {
  bool valid = false;  ///< Migrad + Hesse converged with a positive-definite
                       ///< covariance
  double edm = 0;
  int nfcn = 0;
  Vec u;                   ///< Migrad minimum (u variables)
  Mat cov_u;               ///< Hesse covariance, u variables
  Vec minos_lo, minos_hi;  ///< Minos offsets (lo < 0 < hi), u variables
  std::vector<bool> minos_valid;
};

/// Migrad from u_hat, then Hesse and (optionally) Minos for every parameter.
MinuitErrors minuit2_errors(LeastSquaresProblem& problem, const Vec& u_hat,
                            bool minos, const MinuitOptions& mo = {});

}  // namespace qgem
