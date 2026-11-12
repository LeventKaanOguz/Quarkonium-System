#pragma once
// Special functions used by the closed-form matrix elements (notes 7.5, 11.10).
#include <span>

namespace qgem {

/// e^{x^2} erfc(x), accurate to ~1e-16 relative on [-26, inf).
double erfcx(double x);

/// Int_0^inf r^p e^{-a r^2} dr = Gamma((p+1)/2) / (2 a^{(p+1)/2}),  a > 0, p >
/// -1. Evaluated with lgamma so large p / extreme a do not overflow.
double master_integral(double p, double a);

/// Q_k(a, mu) = Int_0^inf r^k e^{-a r^2 - mu r} dr for k = 0..kmax (kmax+1
/// values in Q). Upward recursion where it is stable, downward ratio recursion
/// otherwise (notes Rem. 7.9).
void screened_moments(int kmax, double a, double mu, std::span<double> Q);

/// Spherical Bessel function j_l(x) for integer l >= 0 and any x >= 0 (series
/// for small x, upward recurrence otherwise; std::sph_bessel refuses large
/// arguments).
double sph_bessel_j(int l, double x);

}  // namespace qgem
