#pragma once
// Minimal forward-mode dual number over the (up to) 8 model parameters. Used
// only for the scalar prefactors of the Hamiltonian terms, so that dH/dtheta
// (Feynman-Hellmann, notes 7.11) follows from the very same expressions that
// build H, with no hand-copied derivative to drift out of sync.
#include <array>
#include <cmath>

namespace qgem {

inline constexpr int kMaxP = 24;

struct Dual {
  double v = 0.0;
  std::array<double, kMaxP> d{};
  Dual() = default;
  Dual(double x) : v(x) {}  // NOLINT: implicit on purpose (constants)
  static Dual var(double x, int i) {
    Dual r(x);
    r.d[i] = 1.0;
    return r;
  }
};

inline Dual operator+(const Dual& a, const Dual& b) {
  Dual r(a.v + b.v);
  for (int i = 0; i < kMaxP; ++i) r.d[i] = a.d[i] + b.d[i];
  return r;
}
inline Dual operator-(const Dual& a, const Dual& b) {
  Dual r(a.v - b.v);
  for (int i = 0; i < kMaxP; ++i) r.d[i] = a.d[i] - b.d[i];
  return r;
}
inline Dual operator-(const Dual& a) { return Dual(0.0) - a; }
inline Dual operator*(const Dual& a, const Dual& b) {
  Dual r(a.v * b.v);
  for (int i = 0; i < kMaxP; ++i) r.d[i] = a.d[i] * b.v + a.v * b.d[i];
  return r;
}
inline Dual operator/(const Dual& a, const Dual& b) {
  Dual r(a.v / b.v);
  for (int i = 0; i < kMaxP; ++i) r.d[i] = (a.d[i] - r.v * b.d[i]) / b.v;
  return r;
}
inline Dual pow(const Dual& a, double p) {
  Dual r(std::pow(a.v, p));
  const double f = p * std::pow(a.v, p - 1.0);
  for (int i = 0; i < kMaxP; ++i) r.d[i] = f * a.d[i];
  return r;
}
inline Dual sqrt(const Dual& a) { return pow(a, 0.5); }

}  // namespace qgem
