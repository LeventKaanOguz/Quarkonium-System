#pragma once
// Minimal JSON writer (objects, arrays, numbers, strings): results carry their
// configuration and seed, so every run is reproducible (notes R9). No external
// dependency.
#include <cmath>
#include <iomanip>
#include <sstream>
#include <string>
#include <vector>

#include "core/types.hpp"

namespace qgem::json {

inline std::string quote(const std::string& s) {
  std::string o = "\"";
  for (char c : s) {
    if (c == '"' || c == '\\') o += '\\';
    if (c == '\n')
      o += "\\n";
    else
      o += c;
  }
  return o + "\"";
}
inline std::string num(double x) {
  if (!std::isfinite(x)) return "null";
  std::ostringstream ss;
  ss << std::setprecision(12) << x;
  return ss.str();
}
inline std::string array(const Vec& v) {
  std::string s = "[";
  for (int i = 0; i < v.size(); ++i) s += (i ? ", " : "") + num(v(i));
  return s + "]";
}
inline std::string array(const Mat& m) {
  std::string s = "[";
  for (int i = 0; i < m.rows(); ++i)
    s += (i ? ", " : "") + array(Vec(m.row(i).transpose()));
  return s + "]";
}
inline std::string array(const std::vector<std::string>& v) {
  std::string s = "[";
  for (std::size_t i = 0; i < v.size(); ++i) s += (i ? ", " : "") + quote(v[i]);
  return s + "]";
}

/// Object builder: pairs of already-serialized values, rendered with one key
/// per line.
class Object {
 public:
  Object& put(const std::string& key, const std::string& raw) {
    items_.emplace_back(key, raw);
    return *this;
  }
  Object& str(const std::string& k, const std::string& v) {
    return put(k, quote(v));
  }
  Object& number(const std::string& k, double v) { return put(k, num(v)); }
  Object& integer(const std::string& k, long v) {
    return put(k, std::to_string(v));
  }
  std::string render(int indent = 0) const {
    const std::string pad(indent + 2, ' '), end(indent, ' ');
    std::string s = "{\n";
    for (std::size_t i = 0; i < items_.size(); ++i)
      s += pad + quote(items_[i].first) + ": " + items_[i].second +
           (i + 1 < items_.size() ? ",\n" : "\n");
    return s + end + "}";
  }

 private:
  std::vector<std::pair<std::string, std::string>> items_;
};

}  // namespace qgem::json
