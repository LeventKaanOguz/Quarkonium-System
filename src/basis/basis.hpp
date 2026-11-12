#pragma once
// Gaussian GEM radial basis (notes 7.2, 7.10).
//   real range:     R_nl(r) = N_nl r^l exp(-nu_n r^2),  nu_n = 1/r_n^2, r_n =
//   r1 q^{n-1}, q = (rmax/r1)^{1/(n-1)} complex range:  u_n = r^l exp(-eta_n (1
//   + i omega) r^2), eta_n = 1/rho_n^2 (rho_n geometric); the real functions
//                   c_n = Re u_n and s_n = -Im u_n carry the oscillation that
//                   radial excitations need.
// A real set plus a complex-range set is one variational space
// (Hylleraas-Undheim-MacDonald bounds hold).
#include <vector>

#include "core/types.hpp"

namespace qgem {

/// Plain aggregate on purpose (designated initializers: {.n = 20, .r1 = 0.05,
/// .rmax = 30}).
struct BasisSpec {
  int n = 30;
  double r1 = 0.02;    ///< shortest range, GeV^-1
  double rmax = 30.0;  ///< longest range, GeV^-1
  double eta = 1e-12;  ///< orthogonalization threshold
  // optional complex-range set; nc = 0 means real Gaussians only
  int nc = 0;
  double rho1 = 1.0;
  double rhomax = 20.0;
  double omega = 0.5;  ///< phase at r = rho_n: how many oscillations fit inside
                       ///< the envelope
};

enum class Fn { Real, Cos, Sin };

struct Basis {
  int L = 0;
  BasisSpec spec;
  Vec nu;     ///< real exponents nu_n (real functions), eta_n (complex-range
              ///< functions)
  Vec omega;  ///< phase parameter of each function (0 for real Gaussians)
  std::vector<Fn> kind;  ///< Real, Cos or Sin
  Vec lognorm2;          ///< ln N^2: closed form for real Gaussians, from <f|f>
                         ///< otherwise

  int size() const { return static_cast<int>(nu.size()); }
  bool is_real() const { return spec.nc == 0; }
  double norm(int i) const;
};

Basis make_basis(int L, const BasisSpec& spec);

/// R(r) = sum_n c_n f_n(r).
double radial_value(const Basis& b, const Vec& c, double r);

/// Coefficient of r^L at the origin per function: R(r) ~ (c.v) r^L (S waves:
/// R(0) = c.v).
Vec origin_coefficients(const Basis& b);

}  // namespace qgem
