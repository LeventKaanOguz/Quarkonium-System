#include "obs/observables.hpp"

#include <algorithm>
#include <cmath>
#include <numbers>
#include <stdexcept>

#include "math/quadrature.hpp"

namespace qgem {
namespace {

int block_offset(const SolvedChannel& s, int comp) {
  int o = 0;
  for (int i = 0; i < comp; ++i) o += s.mats.block_size[i];
  return o;
}
Vec block_vec(const SolvedChannel& s, int level, int comp) {
  return s.levels.at(level).c.segment(block_offset(s, comp),
                                      s.mats.block_size[comp]);
}
double fact(int n) { return std::tgamma(n + 1.0); }

}  // namespace

double moment(const SolvedChannel& s, int level, MatrixCache& cache, double k) {
  double sum = 0;
  for (std::size_t i = 0; i < s.channel.comps.size(); ++i) {
    const int L = s.channel.comps[i].L;
    const Vec c = block_vec(s, level, static_cast<int>(i));
    sum += c.dot(qgem::r_power(cache.basis(L), cache.basis(L), k) * c);
  }
  return sum;
}

double R0_squared(const SolvedChannel& s, int level, MatrixCache& cache) {
  for (std::size_t i = 0; i < s.channel.comps.size(); ++i)
    if (s.channel.comps[i].L == 0) {
      const double R0 = block_vec(s, level, static_cast<int>(i))
                            .dot(origin_coefficients(cache.basis(0)));
      return R0 * R0;
    }
  return 0.0;
}

double origin_derivative(const SolvedChannel& s, int level, MatrixCache& cache,
                         int L) {
  for (std::size_t i = 0; i < s.channel.comps.size(); ++i)
    if (s.channel.comps[i].L == L)
      return fact(L) * block_vec(s, level, static_cast<int>(i))
                           .dot(origin_coefficients(cache.basis(L)));
  return 0.0;
}

WaveFunction wave_function(const SolvedChannel& s, int level,
                           MatrixCache& cache) {
  WaveFunction w;
  for (std::size_t i = 0; i < s.channel.comps.size(); ++i) {
    const int L = s.channel.comps[i].L;
    const Vec c = block_vec(s, level, static_cast<int>(i));
    w.L.push_back(L);
    w.S.push_back(s.channel.comps[i].S);
    w.c.push_back(c);
    w.weight.push_back(c.dot(cache.overlap(L) * c));
    w.norm += w.weight.back();
  }
  w.dominant = static_cast<int>(
      std::max_element(w.weight.begin(), w.weight.end()) - w.weight.begin());
  const int L = w.L[w.dominant];
  const Basis& b = cache.basis(L);
  const Vec& c = w.c[w.dominant];
  const double rmax = 4 * std::sqrt(moment(s, level, cache, 2));
  const int n = 4000;
  std::vector<double> R(n);
  double big = 0;
  for (int k = 0; k < n; ++k) {
    R[k] = radial_value(b, c, rmax * (k + 1) / n);
    big = std::max(big, std::fabs(R[k]));
  }
  w.nodes = 0;
  int last = 0;
  for (int k = 0; k < n; ++k) {
    if (std::fabs(R[k]) <= 1e-3 * big) continue;
    const int sg = R[k] > 0 ? 1 : -1;
    if (last != 0 && sg != last) ++w.nodes;
    last = sg;
  }
  return w;
}

double radial_overlap(const SolvedChannel& f, int lf, const SolvedChannel& i,
                      int li, MatrixCache& cache, double k) {
  double sum = 0;
  for (std::size_t a = 0; a < f.channel.comps.size(); ++a)
    for (std::size_t b = 0; b < i.channel.comps.size(); ++b) {
      if (f.channel.comps[a].S != i.channel.comps[b].S) continue;
      const int Lf = f.channel.comps[a].L, Li = i.channel.comps[b].L;
      const Vec cf = block_vec(f, lf, static_cast<int>(a)),
                ci = block_vec(i, li, static_cast<int>(b));
      sum += cf.dot(qgem::r_power(cache.basis(Lf), cache.basis(Li), k) * ci);
    }
  return sum;
}

double radial_overlap_j0(const SolvedChannel& f, int lf, const SolvedChannel& i,
                         int li, MatrixCache& cache, double kp) {
  // <f| j0(kr/2) |i> = Int r^2 R_f R_i j0(k r/2) dr, S-wave components,
  // composite Gauss-Legendre.
  double sum = 0;
  for (std::size_t a = 0; a < f.channel.comps.size(); ++a)
    for (std::size_t b = 0; b < i.channel.comps.size(); ++b) {
      const auto &ca = f.channel.comps[a], &cb = i.channel.comps[b];
      if (ca.S != cb.S || ca.L != cb.L) continue;
      const Basis& B = cache.basis(ca.L);
      const Vec cf = block_vec(f, lf, static_cast<int>(a)),
                ci = block_vec(i, li, static_cast<int>(b));
      const double nu_min = B.nu.minCoeff();
      const double rc = std::sqrt(
          45.0 /
          (2 * nu_min));  // beyond this the broadest Gaussian squared < e^-45
      const int panels = 64;
      const QuadRule& q = gauss_legendre01(24);
      double acc = 0;
      for (int p = 0; p < panels; ++p) {
        const double lo = rc * p / panels, w = rc / panels;
        for (int m = 0; m < q.x.size(); ++m) {
          const double r = lo + w * q.x(m);
          const double x = 0.5 * kp * r;
          const double j0 = (x < 1e-8) ? 1.0 : std::sin(x) / x;
          acc += q.w(m) * w * r * r * radial_value(B, cf, r) *
                 radial_value(B, ci, r) * j0;
        }
      }
      sum += acc;
    }
  return sum;
}

double p2_expectation(const SolvedChannel& s, int level, MatrixCache& cache) {
  double sum = 0;
  for (std::size_t i = 0; i < s.channel.comps.size(); ++i) {
    const int L = s.channel.comps[i].L;
    const Vec c = block_vec(s, level, static_cast<int>(i));
    sum += c.dot(cache.p2(L) * c);
  }
  return sum;
}

double p4_expectation(const SolvedChannel& s, int level, MatrixCache& cache) {
  double sum = 0;
  for (std::size_t i = 0; i < s.channel.comps.size(); ++i) {
    const int L = s.channel.comps[i].L;
    const Vec c = block_vec(s, level, static_cast<int>(i));
    const Vec pc = cache.p2(L) * c;
    sum += pc.dot(cache.overlap(L).ldlt().solve(pc));
  }
  return sum;
}

double p4_from_potential(const Vec& c, double E, const Mat& V, const Mat& V2,
                         double mr) {
  return 4 * mr * mr * (E * E - 2 * E * c.dot(V * c) + c.dot(V2 * c));
}

double width_ee(double R0sq, double M, double eQ, double als, bool qcd) {
  const double w = 4 * kAlphaEM * kAlphaEM * eQ * eQ * R0sq / (M * M);
  return qcd ? w * (1 - 16 * als / (3 * std::numbers::pi)) : w;
}
double width_gammagamma_S(double R0sq, double M, double eQ, double als,
                          bool qcd) {
  const double w = 12 * kAlphaEM * kAlphaEM * std::pow(eQ, 4) * R0sq / (M * M);
  return qcd ? w * (1 - (20 - std::numbers::pi * std::numbers::pi) * als /
                            (3 * std::numbers::pi))
             : w;
}
double width_gammagamma_P(int J, double Rp0sq, double mQ, double eQ) {
  if (J != 0 && J != 2)
    throw std::invalid_argument("width_gammagamma_P: J in {0, 2}");
  const double f = (J == 0) ? 27.0 : 36.0 / 5.0;
  return f * std::pow(eQ, 4) * kAlphaEM * kAlphaEM * Rp0sq / std::pow(mQ, 4);
}
double width_ee_D(double Rpp0sq, double M, double mQ, double eQ) {
  return 25 * eQ * eQ * kAlphaEM * kAlphaEM * Rpp0sq /
         (2 * std::pow(mQ, 4) * M * M);
}
double width_ee_SD(double RS0, double RD2, double M, double mQ, double eQ) {
  const double A = RS0 + 5 * RD2 / (2 * std::numbers::sqrt2 * mQ * mQ);
  return 4 * kAlphaEM * kAlphaEM * eQ * eQ * A * A / (M * M);
}
double qcd_factor(const std::string& process, double a) {
  const double pi = std::numbers::pi;
  if (process == "ee" || process == "gg_P2") return 1 - 16 * a / (3 * pi);
  if (process == "gg_S") return 1 - (20 - pi * pi) * a / (3 * pi);
  if (process == "gg_P0") return 1 + (pi * pi / 3 - 28.0 / 9) * a / pi;
  throw std::invalid_argument("qcd_factor: unknown process " + process);
}
double decay_constant(double R0sq, double M) {
  return std::sqrt(3 * R0sq / (std::numbers::pi * M));
}

double photon_energy(double Mi, double Mf) {
  return (Mi * Mi - Mf * Mf) / (2 * Mi);
}

namespace {
double tri(int a, int b, int c) {  // Delta(abc)
  return std::sqrt(fact(a + b - c) * fact(a - b + c) * fact(-a + b + c) /
                   fact(a + b + c + 1));
}
bool triad(int a, int b, int c) { return c >= std::abs(a - b) && c <= a + b; }
}  // namespace

double sixj(int j1, int j2, int j3, int j4, int j5, int j6) {
  if (!triad(j1, j2, j3) || !triad(j1, j5, j6) || !triad(j4, j2, j6) ||
      !triad(j4, j5, j3))
    return 0.0;
  const double pre =
      tri(j1, j2, j3) * tri(j1, j5, j6) * tri(j4, j2, j6) * tri(j4, j5, j3);
  const int a1 = j1 + j2 + j3, a2 = j1 + j5 + j6, a3 = j4 + j2 + j6,
            a4 = j4 + j5 + j3;
  const int b1 = j1 + j2 + j4 + j5, b2 = j2 + j3 + j5 + j6,
            b3 = j3 + j1 + j6 + j4;
  double sum = 0;
  for (int t = std::max({a1, a2, a3, a4}); t <= std::min({b1, b2, b3}); ++t)
    sum += ((t % 2) ? -1.0 : 1.0) * fact(t + 1) /
           (fact(t - a1) * fact(t - a2) * fact(t - a3) * fact(t - a4) *
            fact(b1 - t) * fact(b2 - t) * fact(b3 - t));
  return pre * sum;
}

double e1_angular_factor(int Li, int Ji, int Lf, int Jf, int S) {
  const double s6 = sixj(Lf, Jf, S, Ji, Li, 1);
  return std::max(Li, Lf) * (2 * Jf + 1) * s6 * s6;
}
double width_E1(double C, double eQ, double k, double r) {
  return 4.0 / 3.0 * C * eQ * eQ * kAlphaEM * k * k * k * r * r;
}
double width_M1(double eQ, double mQ, double k, double j0) {
  return 4.0 / 3.0 * eQ * eQ * kAlphaEM / (mQ * mQ) * k * k * k * j0 * j0;
}

}  // namespace qgem
