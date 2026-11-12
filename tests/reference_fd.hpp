#pragma once
// Independent finite-difference solver of the radial problem for tests:
//   -u''/(2 mr) + [V(r) + l(l+1)/(2 mr r^2)] u = E u,  u(0) = u(rmax) = 0
// on a uniform grid; eigenvalues by Sturm-sequence bisection (no GEM code
// involved), then Richardson extrapolation in the step h.
#include <algorithm>
#include <cmath>
#include <functional>
#include <vector>

namespace qgem::test {

inline int sturm_count(const std::vector<double>& d, double e,
                       double x) {  // eigenvalues < x
  int cnt = 0;
  double q = d[0] - x;
  if (q < 0) ++cnt;
  for (std::size_t i = 1; i < d.size(); ++i) {
    if (q == 0) q = 1e-300;
    q = d[i] - x - e * e / q;
    if (q < 0) ++cnt;
  }
  return cnt;
}

inline double fd_level_grid(const std::function<double(double)>& V, int l,
                            double mr, double rmax, int N, int k) {
  const double h = rmax / N;
  std::vector<double> d(N - 1);
  double lo = 1e300, hi = -1e300;
  const double e = -1.0 / (2 * mr * h * h);
  for (int i = 1; i < N; ++i) {
    const double r = i * h;
    d[i - 1] = 1.0 / (mr * h * h) + V(r) + l * (l + 1) / (2 * mr * r * r);
    lo = std::min(lo, d[i - 1] + 2 * e);  // Gershgorin: |e| off-diagonals
    hi = std::max(hi, d[i - 1] - 2 * e);
  }
  lo -= 1;
  hi += 1;
  for (int it = 0; it < 200; ++it) {
    const double mid = 0.5 * (lo + hi);
    if (sturm_count(d, e, mid) >= k + 1)
      hi = mid;
    else
      lo = mid;
    if (hi - lo < 1e-14 * std::max(1.0, std::fabs(mid))) break;
  }
  return 0.5 * (lo + hi);
}

/// Richardson extrapolation h -> 0 from N, 2N, 4N assuming E(h) = E0 + c1 h +
/// c2 h^2.
inline double fd_level(const std::function<double(double)>& V, int l, double mr,
                       double rmax, int N, int k) {
  const double e1 = fd_level_grid(V, l, mr, rmax, N, k),
               e2 = fd_level_grid(V, l, mr, rmax, 2 * N, k),
               e3 = fd_level_grid(V, l, mr, rmax, 4 * N, k);
  // E0 = (8 e3 - 6 e2 + e1)/3 eliminates h and h^2 terms for step ratios 1:2:4.
  return (8 * e3 - 6 * e2 + e1) / 3.0;
}

}  // namespace qgem::test
