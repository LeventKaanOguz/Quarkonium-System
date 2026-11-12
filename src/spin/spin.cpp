#include "spin/spin.hpp"

#include <cmath>
#include <stdexcept>

namespace qgem {

double s1s2(int S) {
  if (S != 0 && S != 1) throw std::invalid_argument("s1s2: S in {0, 1}");
  return 0.5 * S * (S + 1) - 0.75;
}

double ls(int L, int S, int J) {
  return 0.5 * (J * (J + 1) - L * (L + 1) - S * (S + 1));
}

double s12_diag(int L, int J) {
  if (L == 0) return 0.0;
  if (J == L + 1) return -2.0 * L / (2.0 * L + 3.0);
  if (J == L) return 2.0;
  if (J == L - 1) return -2.0 * (L + 1) / (2.0 * L - 1.0);
  throw std::invalid_argument("s12_diag: |L - 1| <= J <= L + 1");
}

double s12_offdiag(int J) {
  if (J < 1) throw std::invalid_argument("s12_offdiag: J >= 1");
  return 6.0 * std::sqrt(static_cast<double>(J) * (J + 1)) / (2.0 * J + 1.0);
}

double ls_minus(int L) { return std::sqrt(static_cast<double>(L) * (L + 1)); }

namespace {
const char* kLetters = "SPDFGHIK";
std::string comp_label(const Component& c) {
  return std::to_string(2 * c.S + 1) + kLetters[c.L];
}
}  // namespace

std::string Channel::label() const {
  std::string s;
  for (std::size_t i = 0; i < comps.size(); ++i) {
    if (i) s += "-";
    s += comp_label(comps[i]);
    // J only once for a single component; for a coupled pair the J of the first
    // is appended at the end
  }
  if (comps.size() == 1)
    return std::to_string(2 * comps[0].S + 1) + kLetters[comps[0].L] +
           std::to_string(J);
  return comp_label(comps[0]) + std::to_string(J) + "-" + comp_label(comps[1]) +
         std::to_string(J);
}

std::vector<Channel> enumerate_channels(int Lmax, bool equal_masses) {
  if (Lmax < 0 || Lmax >= 8)
    throw std::invalid_argument("enumerate_channels: 0 <= Lmax < 8");
  std::vector<Channel> out;
  auto parity = [](int L) { return (L % 2 == 0) ? -1 : 1; };
  auto charge = [&](int L, int S) {
    return equal_masses ? ((L + S) % 2 == 0 ? 1 : -1) : 0;
  };
  for (int J = 0; J <= Lmax + 1; ++J) {
    // S = 0, L = J and S = 1, L = J (the latter needs J >= 1: |L - S| <= J);
    // coupled for unequal masses (C is not a symmetry there), and J = 0 then
    // has the singlet only.
    if (J <= Lmax) {
      if (equal_masses) {
        out.push_back({J, parity(J), charge(J, 0), {{J, 0}}});
        if (J >= 1) out.push_back({J, parity(J), charge(J, 1), {{J, 1}}});
      } else if (J == 0) {
        out.push_back({0, parity(0), 0, {{0, 0}}});
      } else {
        out.push_back({J, parity(J), 0, {{J, 0}, {J, 1}}});
      }
    }
    // S = 1, L = J-1, J+1 (same parity (-1)^J)
    const int Lm = J - 1, Lp = J + 1;
    const int P = parity(Lp);  // = parity(Lm)
    std::vector<Component> comps;
    if (Lm >= 0 && Lm <= Lmax) comps.push_back({Lm, 1});
    if (Lp <= Lmax) comps.push_back({Lp, 1});
    if (!comps.empty()) out.push_back({J, P, charge(comps[0].L, 1), comps});
  }
  return out;
}

}  // namespace qgem
