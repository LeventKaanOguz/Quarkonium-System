#pragma once
// Strong coupling alpha_s(Q^2) (notes 5) and the Coulomb-like potentials it
// generates.
//
// Momentum-space potential V(q) = -(16 pi/3) alpha_s(q^2)/q^2, hence
//   V(r)  = -(4/3)(2/pi) Int_0^inf dq sin(qr)/(qr) alpha_s(q^2) (3.5)
// and, after Gaussian smearing (multiply alpha_s by exp(-q^2/4 sigma^2)), every
// integrand is damped, so the same transform gives all the spin structures of
// the Breit-Fermi terms:
//   G(r)       = -(8/(3 pi r)) Int dq A(q) sin(qr)/q            A(q) =
//   alpha_s(q^2) exp(-q^2/4 sigma^2) G'(r)/r    = -(8/(3 pi r^2)) [ Int dq A
//   cos(qr) - (1/r) Int dq A sin(qr)/q ] lap G(r)   = (8/(3 pi r)) Int dq q A
//   sin(qr)                (the generalized contact term; = (16 pi/3) sum a_k
//   rho_k) T(r)       = G'/r - G'' = 3 G'/r - lap G                    (tensor)
// For the Godfrey-Isgur Gaussian sum all four are analytic (erf and Gaussians):
// the numerical transform is then a check, and the analytic form is what the
// models use.
#include <complex>
#include <functional>
#include <memory>
#include <string>
#include <vector>

#include "core/types.hpp"

namespace qgem {

double beta0(int nf);  ///< 11 - 2 nf / 3 (eq. 5.1)
/// The MS-bar QCD coupling at the scale mu (GeV): two-loop running from
/// alpha_s(M_Z) = 0.1180 (PDG), with n_f = 5 above m_b = 4.18 GeV, 4 down to
/// m_c = 1.27 GeV and 3 below (continuous at the thresholds, exact at two
/// loops). For the QCD corrections of the annihilation widths, not for the
/// potential.
double alpha_s_msbar(double mu);

class RunningCoupling {
 public:
  virtual ~RunningCoupling() = default;
  virtual std::string name() const = 0;
  virtual double alpha(double Q2) const = 0;  ///< Q^2 >= 0, GeV^2
  /// alpha_s at complex q^2 (for the contour-rotated transform of unsmeared
  /// potentials); false if unavailable.
  virtual bool alpha_complex(cplx /*Q2*/, cplx& /*out*/) const { return false; }
};

std::unique_ptr<RunningCoupling> make_constant_coupling(double alpha);
/// 4 pi / (beta0 ln((Q^2 + 4 Mg^2)/Lambda^2)) (eq. 5.4 with the 4 Mg^2
/// variant), finite at Q = 0 for 2 Mg > Lambda.
std::unique_ptr<RunningCoupling> make_frozen_coupling(double Lambda, double Mg,
                                                      int nf);
/// Shirkov-Solovtsov analytic perturbation theory (5.5); alpha(0) = 4 pi /
/// beta0.
std::unique_ptr<RunningCoupling> make_apt_coupling(double Lambda, int nf);
/// Richardson (5.6); confining: V(r) -> b r with b = 8 pi Lambda^2 / (33 - 2
/// nf).
std::unique_ptr<RunningCoupling> make_richardson_coupling(double Lambda,
                                                          int nf);
/// sum_k alpha_k exp(-Q^2 / (4 gamma_k^2)) (5.7).
std::unique_ptr<RunningCoupling> make_gaussian_sum_coupling(
    std::vector<double> alpha, std::vector<double> gamma);

/// Godfrey-Isgur's three Gaussians (notes 5.2.6): alpha = (0.25, 0.15, 0.20),
/// gamma = (1/2, sqrt(10)/2, sqrt(1000)/2) GeV.
std::vector<double> gi_alpha();
std::vector<double> gi_gamma();

/// The *shape* alpha(Q^2)/alpha(0) of a coupling that is finite at Q = 0,
/// approximated by sum_k w_k exp(-Q^2/(4 gamma_k^2)) on fixed ranges gamma_k:
/// non-negative least squares (Lawson-Hanson) on relative residuals over a
/// log grid Q in [0.02, qmax] GeV, with sum_k w_k = 1 enforced (so alpha(0) is
/// the model's saturation value alpha_s). Every Coulomb-like matrix element of
/// such a sum is closed-form (erf/r and Gaussians), which is how models A/B
/// take an arbitrary finite running coupling.
struct GaussianShape {
  std::vector<double> w, gamma;
  double max_rel_err = 0;  ///< max |fit/shape - 1| for 0.05 < Q < qmax GeV
                           ///< (below, r > 4 fm, confinement dominates)
};
GaussianShape fit_gaussian_shape(const RunningCoupling& a,
                                 std::vector<double> gamma, double qmax = 50.0);
/// Default ranges for fit_gaussian_shape: 0.25 * 2^k GeV, k = -4..8
/// (1/64 to 64 GeV: the APT shape falls steeply near Q = 0).
std::vector<double> default_shape_gamma();

/// V(r) = -(4/3)(2/pi) Int dq sin(qr)/(qr) alpha_s(q^2) for forms with an
/// analytic continuation (contour rotation: plain quadrature on [0, q0],
/// rotated Gauss-Laguerre beyond). Throws if the form has no continuation.
double static_potential(const RunningCoupling& a, double r);

/// The smeared Coulomb-like functions for one coupling and one smearing width
/// (all by numerical transform).
struct SmearedCoulombValues {
  double G, Gp_over_r, lap, T;
};
SmearedCoulombValues smeared_coulomb_transform(const RunningCoupling& a,
                                               double sigma, double r);

/// Cubic-spline tables of the four functions on a log-spaced grid (used by
/// `local` matrix elements). Regular at the origin, so values below r_min are
/// clamped to r_min.
class SmearedCoulombTable {
 public:
  SmearedCoulombTable(const RunningCoupling& a, double sigma,
                      double rmin = 2e-3, double rmax = 80.0, int npts = 500);
  double G(double r) const { return eval(0, r); }
  double Gp_over_r(double r) const { return eval(1, r); }
  double lap(double r) const { return eval(2, r); }
  double T(double r) const { return eval(3, r); }

 private:
  double eval(int which, double r) const;
  std::vector<double> x_;  // ln r
  std::vector<std::vector<double>> y_,
      m_;  // values and second derivatives of the spline, per function
};

}  // namespace qgem
