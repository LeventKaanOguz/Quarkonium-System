#pragma once
// A dataset: measured masses attached to model levels, with the data covariance
// (independent experimental errors, optional common systematics within a group,
// notes Example 9.1).
#include <string>
#include <vector>

#include "core/types.hpp"
#include "data/pdg.hpp"
#include "solve/spectrum.hpp"

namespace qgem {

struct Datum {
  std::string label;  ///< PDG name
  long mcid = 0;
  StateId state;          ///< which model level this datum measures
  double y = 0;           ///< GeV
  double err = 0;         ///< GeV, independent experimental error
  std::string sys_group;  ///< states sharing a common systematic (e.g.
                          ///< chi_bJ(1P)); empty: none
  double sys_common = 0;  ///< GeV
  bool in_fit = true;
  /// relative size of the theory floor for this state: its floor is
  /// floor_w * sigma_th (1 = the uniform floor)
  double floor_w = 1.0;
};

struct Dataset {
  std::vector<Datum> d;
  /// Experimental covariance of the fitted data (diagonal sigma_exp^2 + common
  /// systematics), order of in_fit data.
  Mat V_exp() const;
  std::vector<const Datum*> fitted() const;
  int n_fit() const;
};

/// Assignment file: one state per line, "mcid channel n [fit|predict] [group]",
/// '#' comments. (The assignment mcid -> (channel, radial index) is physics
/// input, kept out of the code.)
Dataset load_dataset(const std::string& states_file, const PdgDatabase& pdg);

}  // namespace qgem
