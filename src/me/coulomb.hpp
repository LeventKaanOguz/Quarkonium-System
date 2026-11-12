#pragma once
// Radial functions of the Gaussian-smeared Coulomb potential g(r) = erf(tau
// r)/r (notes 6.5.1, 6.5.3). Both are regular at the origin; the naive formulas
// cancel catastrophically for tau r < 1, so a power series is used there.
namespace qgem {

/// g'(r)/r  with g = erf(tau r)/r:  tau^3 h(z), z = tau r,  h = (2/sqrt(pi))
/// e^{-z^2}/z^2 - erf(z)/z^3. Limit at the origin: -4 tau^3 / (3 sqrt(pi)).
double coulomb_gp_over_r(double tau, double r);

/// T_tau(r) = g'/r - g''  (notes 6.20). Point limit: -3/r^3; regular, ~ -8
/// tau^5 r^2/(5 sqrt(pi)) at 0.
double coulomb_tensor(double tau, double r);

}  // namespace qgem
