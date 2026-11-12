#pragma once
// Gaussian-smeared confining potentials of Godfrey-Isgur (notes 6.5.4, 6.5.5)
// with the radial derivatives that enter the scalar (Thomas) spin-orbit term.
// Smearing: f~ = rho_sigma * f, rho_sigma = (sigma/sqrt(pi))^3 e^{-sigma^2
// r^2}.
namespace qgem {

/// r~(r) = e^{-s^2 r^2}/(sqrt(pi) s) + (r + 1/(2 s^2 r)) erf(s r) (6.22)
double smeared_r(double sigma, double r);
/// r~'(r)/r, regular at the origin (4 s/(3 sqrt(pi)) there): power series for s
/// r < 1.
double smeared_r_prime_over_r(double sigma, double r);

/// (e^{-mu r})~ (6.24) in the numerically stable erfcx form (6.27), and its
/// r-derivative (6.25).
double smeared_exp(double sigma, double mu, double r);
double smeared_exp_prime(double sigma, double mu, double r);

}  // namespace qgem
