// Where the tests read the PDG data: the full PDG SQLite file when it is in
// data/ (it is not redistributed), otherwise the small CC BY 4.0 fixture
// tests/data/pdg-subset.sqlite (made by tests/data/make_pdg_subset.py), which
// holds every row the reader uses for the states of config/*.states. The two
// give identical numbers, so the fit tests run in the CI as well.
#pragma once
#include <filesystem>
#include <string>

namespace qgem_test {
inline std::string pdg_data_path() {
  const std::string full =
      std::string(QGEM_SOURCE_DIR) + "/data/pdg-2026.0.sqlite";
  if (std::filesystem::exists(full)) return full;
  return std::string(QGEM_SOURCE_DIR) + "/tests/data/pdg-subset.sqlite";
}
inline bool have_pdg_data() { return std::filesystem::exists(pdg_data_path()); }
}  // namespace qgem_test
