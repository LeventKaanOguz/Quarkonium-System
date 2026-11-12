#pragma once
// Radial matrix elements <R_{n'l'}| O |R_{nl}> = Int r^2 R_{n'l'} O(r) R_{nl}
// dr in the real Gaussian basis. Bra and ket may have different L (tensor
// coupling, transitions).
#include <functional>

#include "basis/basis.hpp"
#include "core/types.hpp"

namespace qgem {

/// S_{n'n}; same L only.
Mat overlap(const Basis& b);

/// <p^2> radial matrix, -d^2 - (2/r) d + l(l+1)/r^2.  Kinetic energy is p2(b) /
/// (2 m_r).
Mat p2(const Basis& b);
Mat kinetic_nr(const Basis& b, double mr);

/// <r^k> = N' N Gamma(l_bar + (k+3)/2) / (2 a^{l_bar + (k+3)/2}), l_bar =
/// (l+l')/2. Requires l + l' + 3 + k > 0.
Mat r_power(const Basis& bra, const Basis& ket, double k);

/// <r^k exp(-beta r^2)> (default k = 0), beta >= 0, l + l' + 3 + k > 0.
Mat gaussian(const Basis& bra, const Basis& ket, double beta, double k = 0.0);

/// <erf(tau r)/r>, tau > 0.
Mat erf_over_r(const Basis& bra, const Basis& ket, double tau);

/// <r^k exp(-mu r)> = N' N Q_{l+l'+2+k}(a, mu), mu >= 0, l + l' + 2 + k >= 0.
Mat exp_power(const Basis& bra, const Basis& ket, double mu, int k);

/// Any V(r): Gauss-Legendre on [0, sqrt(40/a)] per pair, cached rule rescaled.
/// V must be finite and smooth on (0, inf); the integrand r^{l+l'+2} V(r) e^{-a
/// r^2} must be integrable.
Mat local(const Basis& bra, const Basis& ket,
          const std::function<double(double)>& V, int nq = 128);

}  // namespace qgem
