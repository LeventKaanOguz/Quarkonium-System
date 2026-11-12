#pragma once
// Read-only access to the PDG SQLite snapshot (notes 11.10). Masses only;
// everything is converted to GeV at this boundary.
#include <optional>
#include <string>

struct sqlite3;

namespace qgem {

struct PdgMass {
  std::string name;
  double value = 0;  ///< GeV
  double err = 0;    ///< GeV: mean of the positive and negative error
  double err_pos = 0, err_neg = 0;  ///< GeV
  double common_err =
      0;  ///< GeV: second ("systematic") error from the display text, 0 if none
  std::string value_type;  ///< PDG value type (AC, FC, V, ...)
};

struct PdgWidth {
  double value_keV = 0;
  double err_keV = 0;    ///< mean of the positive and negative errors
  bool derived = false;  ///< summary-table BR x total width (errors in
                         ///< quadrature, treated as uncorrelated)
};

class PdgDatabase {
 public:
  explicit PdgDatabase(const std::string& path);
  ~PdgDatabase();
  PdgDatabase(const PdgDatabase&) = delete;
  PdgDatabase& operator=(const PdgDatabase&) = delete;

  /// Summary-table mass of the particle with this Monte Carlo id; nullopt if
  /// not listed.
  std::optional<PdgMass> mass(long mcid) const;

  /// Partial width G(X --> mode) in keV (e.g. mode = "e+ e-", "gamma gamma"):
  /// the listed width, else the summary-table branching fraction times the
  /// summary-table total width (`derived`); nullopt if neither, or a limit.
  std::optional<PdgWidth> partial_width(long mcid,
                                        const std::string& mode) const;

 private:
  std::optional<PdgWidth> derived_width(long mcid,
                                        const std::string& mode) const;
  sqlite3* db_ = nullptr;
};

/// "9859.44+-0.42+-0.31" -> {0.42, 0.31}; one error -> {e, 0}. Returns the
/// second error in the unit of the text.
double second_error_from_display(const std::string& display);

}  // namespace qgem
