#pragma once
// Observables from eigenvectors and masses (notes 7.9, 8). Leading order in v^2
// and alpha_s unless a QCD correction factor is requested. Units: GeV, GeV^-1;
// widths are returned in GeV.
#include <string>
#include <vector>

#include "models/cache.hpp"
#include "solve/spectrum.hpp"

namespace qgem {

inline constexpr double kAlphaEM = 1.0 / 137.036;

// ---- wave-function quantities
// -------------------------------------------------------------------

/// <r^k> summed over the components of the state (diagonal blocks).
double moment(const SolvedChannel& s, int level, MatrixCache& cache, double k);

/// |R(0)|^2 of the S-wave component (zero if the state has none).
double R0_squared(const SolvedChannel& s, int level, MatrixCache& cache);

/// Value at the origin of the L-th derivative of the radial function of the
/// component with angular momentum L: R ~ A r^L, so R^{(L)}(0) = L! A (L = 1:
/// R'(0); L = 2: R''(0)).
double origin_derivative(const SolvedChannel& s, int level, MatrixCache& cache,
                         int L);

/// The wave function of one level: per component (L, S) the expansion
/// coefficients c on the normalised basis functions (R_L(r) = sum_n c_n N_nL
/// r^L exp(-nu_n r^2)) and its weight c^T S c; their sum (the norm, 1); and
/// the radial nodes of the dominant component: sign changes of R(r) for r < 4
/// r_rms where |R| > 1e-3 max|R| (the Gaussian tail is excluded).
struct WaveFunction {
  std::vector<int> L, S;
  std::vector<double> weight;
  std::vector<Vec> c;
  double norm = 0;
  int dominant = 0;  ///< index of the component with the largest weight
  int nodes = -1;
};
WaveFunction wave_function(const SolvedChannel& s, int level,
                           MatrixCache& cache);

/// <f| r^k |i> over the matching-spin components (the radial overlap of
/// eq. 7.13, l_bar = (L_f+L_i)/2).
double radial_overlap(const SolvedChannel& f, int lf, const SolvedChannel& i,
                      int li, MatrixCache& cache, double k);

/// <f| j0(k r / 2) |i>  by Gauss-Legendre quadrature of the radial functions (S
/// waves).
double radial_overlap_j0(const SolvedChannel& f, int lf, const SolvedChannel& i,
                         int li, MatrixCache& cache, double k_photon);

/// <p^2> = 2 m_r <T>  summed over components (needs the kinetic matrix: p2 / (2
/// m_r) -> p2 directly).
double p2_expectation(const SolvedChannel& s, int level, MatrixCache& cache);
/// <p^4> = |p^2 psi|^2 in the basis: c^T P2 S^-1 P2 c per partial wave (the
/// projection of p^2 psi onto the basis; converges with it).
double p4_expectation(const SolvedChannel& s, int level, MatrixCache& cache);

/// <p^4> of an exact eigenstate as (2 m_r)^2 <(E - V)^2> = (2 m_r)^2 (E^2 - 2 E
/// <V> + <V^2>) (notes 12.2): converges much faster than the direct p^4 matrix.
/// V and V^2 are given as matrices in the state's basis (single-component
/// states).
double p4_from_potential(const Vec& c, double E, const Mat& V, const Mat& V2,
                         double mr);

// ---- widths and decay constants (GeV)
// ------------------------------------------------------------

/// V -> e+e-:  4 alpha^2 e_Q^2 |R(0)|^2 / M^2 * (1 - 16 alpha_s/(3 pi) if qcd).
double width_ee(double R0sq, double M, double eQ, double alpha_s = 0.0,
                bool qcd = false);
/// 1S0 -> gamma gamma:  12 alpha^2 e_Q^4 |R(0)|^2 / M^2 * (1 - (20 - pi^2)
/// alpha_s/(3 pi) if qcd).
double width_gammagamma_S(double R0sq, double M, double eQ,
                          double alpha_s = 0.0, bool qcd = false);
/// 3PJ -> gamma gamma (J = 0, 2):  27 e_Q^4 alpha^2 |R'(0)|^2 / m^4 (J=0), 36/5
/// (J=2).
double width_gammagamma_P(int J, double Rp0sq, double mQ, double eQ);
/// 3D1 -> e+e-:  25 e_Q^2 alpha^2 |R''(0)|^2 / (2 m^4 M^2).
double width_ee_D(double Rpp0sq, double M, double mQ, double eQ);
/// V -> e+e- of a 3S1, 3D1 or mixed 3S1-3D1 state:  4 alpha^2 e_Q^2 A^2 / M^2
/// with the amplitude A = R_S(0) + 5 R''_D(0) / (2 sqrt2 m^2) (signed radial
/// values; the tensor coupling <3D1|S12|3S1> = +sqrt8 fixes the relative sign).
/// Pure 3S1: width_ee; pure 3D1: width_ee_D.
double width_ee_SD(double RS0, double RD2, double M, double mQ, double eQ);
/// One-loop QCD correction factors of the annihilation widths (Barbieri et
/// al.; Kwong et al., PRD 37, 3210 (1988)): 3S1 -> e+e- 1 - 16 a/(3 pi);
/// 1S0 -> gg 1 - (20 - pi^2) a/(3 pi); 3P0 -> gg 1 + (pi^2/3 - 28/9) a/pi;
/// 3P2 -> gg 1 - 16 a/(3 pi). `process`: "ee", "gg_S", "gg_P0", "gg_P2".
double qcd_factor(const std::string& process, double alpha_s);
/// f_V^2 = 3 |R(0)|^2 / (pi M) (Van Royen-Weisskopf).
double decay_constant(double R0sq, double M);

/// Photon energy in the rest frame of the initial state: k = (Mi^2 - Mf^2)/(2
/// Mi).
double photon_energy(double Mi, double Mf);
/// Wigner 6j symbol {j1 j2 j3; j4 j5 j6} (integer arguments), Racah formula.
double sixj(int j1, int j2, int j3, int j4, int j5, int j6);
/// Angular factor of E1: C_fi = max(L, L') (2J'+1) {L' J' S; J L 1}^2 (J, L
/// initial; J', L' final).
double e1_angular_factor(int Li, int Ji, int Lf, int Jf, int S);
/// E1: (4/3) C_fi <e>^2 alpha k^3 |<f|r|i>|^2 ;  <e> = (m2 e1 - m1 e2)/(m1+m2)
/// -> e_Q for equal masses.
double width_E1(double C_fi, double eQ_eff, double k, double r_overlap);
/// M1 (S-wave n3S1 -> n'1S0): (4/3)(e_Q^2 alpha/m^2) k^3 |<j0(kr/2)>|^2 ;
/// multiply by 3 for 1S0 -> 3S1.
double width_M1(double eQ, double mQ, double k, double j0_overlap);

}  // namespace qgem
