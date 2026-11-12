#pragma once
// Spin-angle algebra of the Breit-Fermi terms (notes 3.5): every rotational
// scalar acts as a number in a coupled state |(L S) J>. Derived in the notes;
// verified in tests against explicit construction of the operators with
// Clebsch-Gordan coefficients.
#include <string>
#include <vector>

namespace qgem {

/// <S1.S2>: -3/4 (S = 0), +1/4 (S = 1).
double s1s2(int S);

/// <L.S> = [J(J+1) - L(L+1) - S(S+1)]/2.
double ls(int L, int S, int J);

/// Diagonal tensor element <S12>, S = 1 only: J = L+1: -2L/(2L+3); J = L: 2; J
/// = L-1: -2(L+1)/(2L-1). Zero for L = 0.
double s12_diag(int L, int J);

/// <(J-1) 1 J| S12 |(J+1) 1 J> = 6 sqrt(J(J+1)) / (2J+1)   (3S1-3D1: 2 sqrt 2),
/// Condon-Shortley phases.
double s12_offdiag(int J);

/// <1L_L| L.(S1 - S2) |3L_L> = sqrt(L(L+1)): the antisymmetric spin-orbit term,
/// m1 != m2 only.
double ls_minus(int L);

/// One partial wave (L, S).
struct Component {
  int L, S;
};

/// A block H c = E S c: one component, or two coupled ones (tensor L = J -+ 1
/// with S = 1; 1L_L - 3L_L for m1 != m2).
struct Channel {
  int J = 0, P = 0, C = 0;  ///< J^{PC}; C = 0 when m1 != m2
  std::vector<Component> comps;
  std::string label() const;  ///< "1S0", "3S1-3D1", "3P0", ...
};

/// All channels with L <= Lmax (J up to Lmax+1 so that e.g. 3D1 appears with
/// 3S1).
std::vector<Channel> enumerate_channels(int Lmax, bool equal_masses);

}  // namespace qgem
