#pragma once
// Exact closed-form radial matrix elements of every smeared potential used by
// the models (notes 6.5, 7.5.5, 7.5.7). No quadrature: each element is a finite
// combination of Gamma functions, erfcx and the elementary integral
//   J_s(a, tau) = Int_0^tau (a + t^2)^{-s} dt,   s = k/2,   exact upward
//   recursion from asinh / atan.
// Real-range bases only (complex-range bases fall back to numerical quadrature
// in me/radial).
#include "basis/basis.hpp"
#include "core/types.hpp"

namespace qgem {

/// J_s(a, tau) for s = 1/2, 1, 3/2, ... (2s integer >= 1), a > 0, tau >= 0.
double inverse_quadratic_integral(int two_s, double a, double tau);

/// <g'/r>, g = erf(tau r)/r (the vector spin-orbit radial function of the
/// smeared Coulomb potential).
Mat gprime_over_r_matrix(const Basis& bra, const Basis& ket, double tau);
/// <T>, T = g'/r - g'' = 3 g'/r + (4 tau^3/sqrt(pi)) e^{-tau^2 r^2}  (tensor
/// radial function).
Mat tensor_matrix(const Basis& bra, const Basis& ket, double tau);
/// <r~> and <r~'/r> for the Gaussian-smeared linear potential (6.22); r~'/r is
/// regular at the origin.
Mat smeared_r_matrix(const Basis& bra, const Basis& ket, double sigma);
Mat smeared_r_prime_over_r_matrix(const Basis& bra, const Basis& ket,
                                  double sigma);
/// <(e^{-mu r})~> and <(e^{-mu r})~'/r> by the transfer lemmas (7.22), (7.23)
/// with exact Taylor jets in a. Equal L.
Mat smeared_exp_matrix(const Basis& b, double sigma, double mu);
Mat smeared_exp_prime_over_r_matrix(const Basis& b, double sigma, double mu);

}  // namespace qgem
