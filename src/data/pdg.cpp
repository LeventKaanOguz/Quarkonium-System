#include "data/pdg.hpp"

#include <sqlite3.h>

#include <cmath>
#include <stdexcept>
#include <vector>

namespace qgem {

double second_error_from_display(const std::string& s) {
  // split on "+-"; first token is the value, then the errors
  std::vector<std::string> tok;
  std::size_t pos = 0;
  while (true) {
    const std::size_t q = s.find("+-", pos);
    if (q == std::string::npos) {
      tok.push_back(s.substr(pos));
      break;
    }
    tok.push_back(s.substr(pos, q - pos));
    pos = q + 2;
  }
  if (tok.size() < 3) return 0.0;
  try {
    return std::stod(tok[2]);
  } catch (...) {
    return 0.0;
  }
}

PdgDatabase::PdgDatabase(const std::string& path) {
  if (sqlite3_open_v2(path.c_str(), &db_, SQLITE_OPEN_READONLY, nullptr) !=
      SQLITE_OK) {
    const std::string msg = db_ ? sqlite3_errmsg(db_) : "out of memory";
    sqlite3_close(db_);
    throw std::runtime_error("cannot open PDG database " + path + ": " + msg);
  }
}

PdgDatabase::~PdgDatabase() { sqlite3_close(db_); }

std::optional<PdgMass> PdgDatabase::mass(long mcid) const {
  // in_summary_table = 1 picks the fit value where PDG lists both a fit and an
  // average.
  static const char* sql =
      "SELECT p.name, d.value, d.error_positive, d.error_negative, "
      "d.unit_text, d.display_value_text, d.value_type "
      "FROM pdgparticle p "
      "JOIN pdgid i ON i.parent_pdgid = p.pdgid AND i.data_type = 'M' "
      "JOIN pdgdata d ON d.pdgid_id = i.id "
      "WHERE p.mcid = ?1 AND d.in_summary_table = 1;";
  sqlite3_stmt* st = nullptr;
  if (sqlite3_prepare_v2(db_, sql, -1, &st, nullptr) != SQLITE_OK)
    throw std::runtime_error(sqlite3_errmsg(db_));
  sqlite3_bind_int64(st, 1, mcid);
  std::optional<PdgMass> out;
  auto text = [&](int c) {
    const unsigned char* p = sqlite3_column_text(st, c);
    return p ? std::string(reinterpret_cast<const char*>(p)) : std::string();
  };
  while (sqlite3_step(st) == SQLITE_ROW) {
    if (text(4) != "MeV")
      continue;  // all masses of interest are stored in MeV; anything else is
                 // skipped
    PdgMass m;
    m.name = text(0);
    m.value = sqlite3_column_double(st, 1) * 1e-3;
    m.err_pos = sqlite3_column_double(st, 2) * 1e-3;
    m.err_neg = sqlite3_column_double(st, 3) * 1e-3;
    m.err = 0.5 * (m.err_pos + m.err_neg);
    m.common_err = second_error_from_display(text(5)) * 1e-3;
    m.value_type = text(6);
    out = m;
  }
  sqlite3_finalize(st);
  return out;
}

std::optional<PdgWidth> PdgDatabase::partial_width(
    long mcid, const std::string& mode) const {
  static const char* sql =
      "SELECT d.value, d.error_positive, d.error_negative, d.unit_text, "
      "d.display_value_text "
      "FROM pdgparticle p "
      "JOIN pdgid i ON i.parent_pdgid = p.pdgid "
      "JOIN pdgdata d ON d.pdgid_id = i.id "
      "WHERE p.mcid = ?1 AND i.description = 'G(' || p.name || ' --> ' || ?2 "
      "|| ')' ORDER BY d.in_summary_table DESC;";
  sqlite3_stmt* st = nullptr;
  if (sqlite3_prepare_v2(db_, sql, -1, &st, nullptr) != SQLITE_OK)
    throw std::runtime_error(sqlite3_errmsg(db_));
  sqlite3_bind_int64(st, 1, mcid);
  sqlite3_bind_text(st, 2, mode.c_str(), -1, SQLITE_TRANSIENT);
  std::optional<PdgWidth> out;
  auto text = [&](int c) {
    const unsigned char* p = sqlite3_column_text(st, c);
    return p ? std::string(reinterpret_cast<const char*>(p)) : std::string();
  };
  while (sqlite3_step(st) == SQLITE_ROW) {
    const std::string unit = text(3), disp = text(4);
    if (disp.find('<') != std::string::npos ||
        disp.find('>') != std::string::npos)
      continue;  // limits
    double f;
    if (unit == "keV")
      f = 1.0;
    else if (unit == "eV")
      f = 1e-3;
    else if (unit == "MeV")
      f = 1e3;
    else
      continue;
    PdgWidth w;
    w.value_keV = sqlite3_column_double(st, 0) * f;
    w.err_keV =
        0.5 * (sqlite3_column_double(st, 1) + sqlite3_column_double(st, 2)) * f;
    out = w;
    break;
  }
  sqlite3_finalize(st);
  if (!out) out = derived_width(mcid, mode);
  return out;
}

// Summary-table branching fraction times summary-table total width, when the
// PDG lists no partial width (chi_c0,2 -> gamma gamma, psi(3770) -> e+e-).
std::optional<PdgWidth> PdgDatabase::derived_width(
    long mcid, const std::string& mode) const {
  auto one = [&](const char* sql) -> std::optional<std::pair<double, double>> {
    sqlite3_stmt* st = nullptr;
    if (sqlite3_prepare_v2(db_, sql, -1, &st, nullptr) != SQLITE_OK)
      throw std::runtime_error(sqlite3_errmsg(db_));
    sqlite3_bind_int64(st, 1, mcid);
    sqlite3_bind_text(st, 2, mode.c_str(), -1, SQLITE_TRANSIENT);
    std::optional<std::pair<double, double>> r;
    if (sqlite3_step(st) == SQLITE_ROW &&
        sqlite3_column_type(st, 0) != SQLITE_NULL) {
      const double v = sqlite3_column_double(st, 0);
      const double e =
          0.5 * (sqlite3_column_double(st, 1) + sqlite3_column_double(st, 2));
      if (v > 0) r = std::make_pair(v, e);
    }
    sqlite3_finalize(st);
    return r;
  };
  const auto br = one(
      "SELECT d.value, d.error_positive, d.error_negative FROM pdgparticle p "
      "JOIN pdgid i ON i.parent_pdgid = p.pdgid JOIN pdgdata d ON d.pdgid_id = "
      "i.id WHERE p.mcid = ?1 AND i.data_type = 'BFX' AND d.in_summary_table "
      "= 1 AND d.limit_type IS NULL AND i.description = p.name || ' --> ' || "
      "?2;");
  // the total width, in MeV (?2 is bound but unused)
  const auto gt = one(
      "SELECT d.value, d.error_positive, d.error_negative FROM pdgparticle p "
      "JOIN pdgid i ON i.parent_pdgid = p.pdgid JOIN pdgdata d ON d.pdgid_id = "
      "i.id WHERE p.mcid = ?1 AND ?2 = ?2 AND i.data_type = 'G' AND "
      "d.in_summary_table = 1 AND d.unit_text = 'MeV' AND i.description = "
      "p.name || ' WIDTH';");
  if (!br || !gt) return std::nullopt;
  PdgWidth w;
  w.value_keV = br->first * gt->first * 1e3;
  w.err_keV =
      w.value_keV * std::hypot(br->second / br->first, gt->second / gt->first);
  w.derived = true;
  return w;
}

}  // namespace qgem
