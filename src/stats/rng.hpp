#pragma once
// Reproducible random numbers: std::mt19937_64 is fully specified by the
// standard, but std::normal_distribution is not, so Gaussians are generated
// here (Box-Muller) to make bootstrap/MCMC bit-identical across standard
// libraries (notes 11.3).
#include <cmath>
#include <cstdint>
#include <numbers>
#include <random>

namespace qgem {

class Rng {
 public:
  explicit Rng(std::uint64_t seed, std::uint64_t stream = 0)
      : eng_(mix(seed, stream)) {}
  double uniform() {
    return (static_cast<double>(eng_() >> 11) + 0.5) *
           (1.0 / 9007199254740992.0);
  }  // (0,1)
  double normal() {
    if (have_) {
      have_ = false;
      return spare_;
    }
    const double u1 = uniform(), u2 = uniform();
    const double r = std::sqrt(-2.0 * std::log(u1));
    spare_ = r * std::sin(2.0 * std::numbers::pi * u2);
    have_ = true;
    return r * std::cos(2.0 * std::numbers::pi * u2);
  }

 private:
  static std::uint64_t mix(std::uint64_t a,
                           std::uint64_t b) {  // splitmix64 of (seed, stream)
    std::uint64_t x = a + 0x9E3779B97F4A7C15ULL * (b + 1);
    x = (x ^ (x >> 30)) * 0xBF58476D1CE4E5B9ULL;
    x = (x ^ (x >> 27)) * 0x94D049BB133111EBULL;
    return x ^ (x >> 31);
  }
  std::mt19937_64 eng_;
  bool have_ = false;
  double spare_ = 0;
};

}  // namespace qgem
