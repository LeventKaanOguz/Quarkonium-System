#pragma once
// Independent reference integrator for tests: composite Gauss-Legendre whose
// nodes come from Newton iteration on the Legendre recurrence (deliberately NOT
// the Golub-Welsch route the library uses), in long double.
#include <cmath>
#include <functional>
#include <vector>

namespace qgem::test {

inline void legendre_nodes(int n, std::vector<long double>& x,
                           std::vector<long double>& w) {
  x.assign(n, 0);
  w.assign(n, 0);
  const long double pi = 3.141592653589793238462643383279502884L;
  for (int i = 0; i < n; ++i) {
    long double z = std::cos(pi * (i + 0.75L) / (n + 0.5L)), dp = 1;
    for (int it = 0; it < 100; ++it) {
      long double p0 = 1, p1 = z;
      for (int k = 2; k <= n; ++k) {
        long double p2 = ((2 * k - 1) * z * p1 - (k - 1) * p0) / k;
        p0 = p1;
        p1 = p2;
      }
      dp = n * (z * p1 - p0) / (z * z - 1);
      long double dz = p1 / dp;
      z -= dz;
      if (std::fabs(dz) < 1e-19L) break;
    }
    x[i] = z;
    w[i] = 2 / ((1 - z * z) * dp * dp);
  }
}

/// Int_a^b f(r) dr with `panels` panels of an `n`-point rule.
inline double integrate(const std::function<double(double)>& f, double a,
                        double b, int panels = 400, int n = 24) {
  static std::vector<long double> x, w;
  static int cached = 0;
  if (cached != n) {
    legendre_nodes(n, x, w);
    cached = n;
  }
  long double sum = 0, h = (static_cast<long double>(b) - a) / panels;
  for (int p = 0; p < panels; ++p) {
    long double lo = a + p * h, mid = lo + h / 2;
    for (int i = 0; i < n; ++i)
      sum += w[i] * f(static_cast<double>(mid + x[i] * h / 2));
  }
  return static_cast<double>(sum * h / 2);
}

}  // namespace qgem::test
