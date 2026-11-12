#pragma once
// Units: GeV and GeV^-1 everywhere inside the program, hbar = c = 1.
// Masses arrive from the PDG file in MeV and are converted once, at the
// boundary.

namespace qgem {

inline constexpr double hbarc_GeV_fm =
    0.1973269804;                      // GeV fm; for printing radii in fm only
inline constexpr double MeV = 1.0e-3;  // 1 MeV in GeV

}  // namespace qgem
