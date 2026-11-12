#include "math/optimize.hpp"

#include <algorithm>
#include <cmath>
#include <numeric>
#include <vector>

namespace qgem {

NelderMeadResult nelder_mead(const std::function<double(const Vec&)>& f,
                             const Vec& x0, const NelderMeadOptions& opt) {
  const int n = static_cast<int>(x0.size());
  std::vector<Vec> X(n + 1, x0);
  for (int i = 0; i < n; ++i) X[i + 1](i) += opt.step;
  std::vector<double> F(n + 1);
  for (int i = 0; i <= n; ++i) F[i] = f(X[i]);
  NelderMeadResult R;
  std::vector<int> idx(n + 1);
  for (int it = 0; it < opt.max_iter; ++it) {
    R.iterations = it + 1;
    std::iota(idx.begin(), idx.end(), 0);
    std::sort(idx.begin(), idx.end(),
              [&](int a, int b) { return F[a] < F[b]; });
    std::vector<Vec> Xs;
    std::vector<double> Fs;
    for (int i : idx) {
      Xs.push_back(X[i]);
      Fs.push_back(F[i]);
    }
    X = Xs;
    F = Fs;
    double size = 0;
    for (int i = 1; i <= n; ++i)
      size = std::max(size, (X[i] - X[0]).lpNorm<Eigen::Infinity>());
    if (size < opt.xtol &&
        std::fabs(F[n] - F[0]) < opt.ftol * (1 + std::fabs(F[0]))) {
      R.converged = true;
      break;
    }
    Vec c = Vec::Zero(n);
    for (int i = 0; i < n; ++i) c += X[i] / n;
    const Vec xr = c + (c - X[n]);
    const double fr = f(xr);
    if (fr < F[0]) {
      const Vec xe = c + 2.0 * (c - X[n]);
      const double fe = f(xe);
      if (fe < fr) {
        X[n] = xe;
        F[n] = fe;
      } else {
        X[n] = xr;
        F[n] = fr;
      }
    } else if (fr < F[n - 1]) {
      X[n] = xr;
      F[n] = fr;
    } else {
      const Vec xc =
          (fr < F[n]) ? Vec(c + 0.5 * (xr - c)) : Vec(c + 0.5 * (X[n] - c));
      const double fc = f(xc);
      if (fc < std::min(fr, F[n])) {
        X[n] = xc;
        F[n] = fc;
      } else {
        for (int i = 1; i <= n; ++i) {
          X[i] = X[0] + 0.5 * (X[i] - X[0]);
          F[i] = f(X[i]);
        }
      }
    }
  }
  int best = static_cast<int>(std::min_element(F.begin(), F.end()) - F.begin());
  R.x = X[best];
  R.f = F[best];
  return R;
}

}  // namespace qgem
