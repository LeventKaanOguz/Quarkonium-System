#include "models/alpha_s.hpp"

#include <Eigen/Dense>
#include <algorithm>
#include <cmath>
#include <numbers>
#include <stdexcept>

#include "math/quadrature.hpp"

namespace qgem {
namespace {
constexpr double kPi = std::numbers::pi;

class Constant final : public RunningCoupling {
 public:
  explicit Constant(double a) : a_(a) {}
  std::string name() const override { return "constant"; }
  double alpha(double) const override { return a_; }
  bool alpha_complex(cplx, cplx& out) const override {
    out = a_;
    return true;
  }

 private:
  double a_;
};

class Frozen final : public RunningCoupling {
 public:
  Frozen(double L, double Mg, int nf)
      : L2_(L * L), M2_(4 * Mg * Mg), b0_(beta0(nf)) {
    if (!(M2_ > L2_))
      throw std::invalid_argument("frozen coupling needs 2 Mg > Lambda");
  }
  std::string name() const override { return "frozen"; }
  double alpha(double Q2) const override {
    return 4 * kPi / (b0_ * std::log((Q2 + M2_) / L2_));
  }
  bool alpha_complex(cplx Q2, cplx& out) const override {
    out = 4.0 * kPi / (b0_ * std::log((Q2 + M2_) / L2_));
    return true;
  }

 private:
  double L2_, M2_, b0_;
};

class Apt final : public RunningCoupling {
 public:
  Apt(double L, int nf) : L2_(L * L), b0_(beta0(nf)) {}
  std::string name() const override { return "APT"; }
  double alpha(double Q2) const override {
    const double x = Q2 / L2_;
    if (std::fabs(x - 1.0) < 1e-6)
      return 4 * kPi / b0_ *
             0.5;  // limit of the bracket at Q^2 = Lambda^2 is 1/2
    return 4 * kPi / b0_ * (1.0 / std::log(x) + 1.0 / (1.0 - x));
  }
  bool alpha_complex(cplx Q2, cplx& out) const override {
    const cplx x = Q2 / L2_;
    out = 4.0 * kPi / b0_ * (1.0 / std::log(x) + 1.0 / (1.0 - x));
    return true;
  }

 private:
  double L2_, b0_;
};

class Richardson final : public RunningCoupling {
 public:
  Richardson(double L, int nf) : L2_(L * L), b0_(beta0(nf)) {}
  std::string name() const override { return "Richardson"; }
  double alpha(double Q2) const override {
    return 4 * kPi / (b0_ * std::log1p(Q2 / L2_));
  }
  bool alpha_complex(cplx Q2, cplx& out) const override {
    out = 4.0 * kPi / (b0_ * std::log(1.0 + Q2 / L2_));
    return true;
  }

 private:
  double L2_, b0_;
};

class GaussianSum final : public RunningCoupling {
 public:
  GaussianSum(std::vector<double> a, std::vector<double> g)
      : a_(std::move(a)), g_(std::move(g)) {
    if (a_.size() != g_.size() || a_.empty())
      throw std::invalid_argument(
          "GaussianSum: matching, non-empty (alpha, gamma)");
  }
  std::string name() const override { return "Gaussian sum"; }
  double alpha(double Q2) const override {
    double s = 0;
    for (std::size_t k = 0; k < a_.size(); ++k)
      s += a_[k] * std::exp(-Q2 / (4 * g_[k] * g_[k]));
    return s;
  }

 private:
  std::vector<double> a_, g_;
};

}  // namespace

double beta0(int nf) { return 11.0 - 2.0 * nf / 3.0; }

double alpha_s_msbar(double mu) {
  // d a / d ln mu^2 = -b0 a^2/(4 pi) - b1 a^3/(16 pi^2), RK4 in t = ln mu^2
  auto run = [](double a, double t0, double t1, int nf) {
    const double b0 = beta0(nf), b1 = 102.0 - 38.0 * nf / 3.0;
    auto f = [&](double x) {
      return -b0 * x * x / (4 * std::numbers::pi) -
             b1 * x * x * x / (16 * std::numbers::pi * std::numbers::pi);
    };
    const int n = 400;
    const double h = (t1 - t0) / n;
    for (int i = 0; i < n; ++i) {
      const double k1 = f(a), k2 = f(a + 0.5 * h * k1),
                   k3 = f(a + 0.5 * h * k2), k4 = f(a + h * k3);
      a += h * (k1 + 2 * k2 + 2 * k3 + k4) / 6;
    }
    return a;
  };
  const double mz = 91.1876, mb = 4.18, mc = 1.27;
  const double t = std::log(mu * mu);
  double a = 0.1180;
  if (mu >= mb) return run(a, std::log(mz * mz), t, 5);
  a = run(a, std::log(mz * mz), std::log(mb * mb), 5);
  if (mu >= mc) return run(a, std::log(mb * mb), t, 4);
  a = run(a, std::log(mb * mb), std::log(mc * mc), 4);
  return run(a, std::log(mc * mc), t, 3);
}
std::unique_ptr<RunningCoupling> make_constant_coupling(double a) {
  return std::make_unique<Constant>(a);
}
std::unique_ptr<RunningCoupling> make_frozen_coupling(double L, double Mg,
                                                      int nf) {
  return std::make_unique<Frozen>(L, Mg, nf);
}
std::unique_ptr<RunningCoupling> make_apt_coupling(double L, int nf) {
  return std::make_unique<Apt>(L, nf);
}
std::unique_ptr<RunningCoupling> make_richardson_coupling(double L, int nf) {
  return std::make_unique<Richardson>(L, nf);
}
std::unique_ptr<RunningCoupling> make_gaussian_sum_coupling(
    std::vector<double> a, std::vector<double> g) {
  return std::make_unique<GaussianSum>(std::move(a), std::move(g));
}
std::vector<double> gi_alpha() { return {0.25, 0.15, 0.20}; }
std::vector<double> gi_gamma() {
  return {0.5, std::sqrt(10.0) / 2, std::sqrt(1000.0) / 2};
}

double static_potential(const RunningCoupling& a, double r) {
  if (!(r > 0)) throw std::invalid_argument("static_potential: r > 0");
  cplx probe;
  if (!a.alpha_complex(cplx(10.0, 0.0), probe))
    throw std::invalid_argument("static_potential: " + a.name() +
                                " has no analytic continuation");
  // V(r) = -(8/(3 pi r)) Int_0^inf dq sin(qr) alpha(q^2)/q
  const double q0 =
      10.0;  // GeV; far above Lambda, alpha is analytic and smooth beyond it
  const QuadRule& gl = gauss_legendre01(200);
  double low = 0;
  for (int k = 0; k < gl.x.size(); ++k) {
    const double q = q0 * gl.x(k);
    low += gl.w(k) * q0 * ((q * r < 1e-9) ? r : std::sin(q * r) / q) *
           a.alpha(q * q);
  }
  // Int_{q0}^inf sin(qr) g(q) dq = Im[ i e^{i q0 r} Int_0^inf e^{-t r} g(q0 + i
  // t) dt ],  g = alpha(q^2)/q
  const QuadRule& lag = gauss_laguerre(120, 0.0);
  cplx acc = 0;
  for (int k = 0; k < lag.x.size(); ++k) {
    const cplx q(q0, lag.x(k) / r);
    cplx al;
    a.alpha_complex(q * q, al);
    acc += lag.w(k) / r * al / q;
  }
  const double high = (cplx(0, 1) * std::exp(cplx(0, q0 * r)) * acc).imag();
  return -(8.0 / (3.0 * kPi * r)) * (low + high);
}

SmearedCoulombValues smeared_coulomb_transform(const RunningCoupling& a,
                                               double sigma, double r) {
  // A(q) = alpha(q^2) exp(-q^2/4 sigma^2) is Gaussian-damped: composite
  // Gauss-Legendre to q_max, panels fine enough for sin(q r) (>= ~4 points per
  // half-wavelength).
  const double qmax = sigma * std::sqrt(4.0 * 42.0);
  const int panels =
      std::max(40, static_cast<int>(std::ceil(qmax * r / kPi * 2.0)));
  const QuadRule& g = gauss_legendre01(12);
  double I_sin_over_q = 0, I_cos = 0, I_q_sin = 0;
  const double h = qmax / panels;
  for (int p = 0; p < panels; ++p)
    for (int k = 0; k < g.x.size(); ++k) {
      const double q = (p + g.x(k)) * h, w = g.w(k) * h;
      const double A = a.alpha(q * q) * std::exp(-q * q / (4 * sigma * sigma));
      const double s = std::sin(q * r);
      I_sin_over_q += w * A * ((q * r < 1e-9) ? r : s / q);
      I_cos += w * A * std::cos(q * r);
      I_q_sin += w * A * q * s;
    }
  SmearedCoulombValues v;
  const double c = 8.0 / (3.0 * kPi);
  v.G = -c * I_sin_over_q / r;
  const double Gp = -c * (I_cos - I_sin_over_q / r) / r;
  v.Gp_over_r = Gp / r;
  v.lap = c * I_q_sin / r;
  v.T = 3.0 * v.Gp_over_r - v.lap;
  return v;
}

SmearedCoulombTable::SmearedCoulombTable(const RunningCoupling& a, double sigma,
                                         double rmin, double rmax, int npts) {
  x_.resize(npts);
  y_.assign(4, std::vector<double>(npts));
  m_.assign(4, std::vector<double>(npts, 0.0));
  const double l0 = std::log(rmin), l1 = std::log(rmax);
  for (int i = 0; i < npts; ++i) {
    x_[i] = l0 + (l1 - l0) * i / (npts - 1);
    const SmearedCoulombValues v =
        smeared_coulomb_transform(a, sigma, std::exp(x_[i]));
    y_[0][i] = v.G;
    y_[1][i] = v.Gp_over_r;
    y_[2][i] = v.lap;
    y_[3][i] = v.T;
  }
  // natural cubic spline second derivatives (tridiagonal solve) for each
  // function
  for (int f = 0; f < 4; ++f) {
    const int n = npts;
    std::vector<double> u(n, 0.0), m(n, 0.0);
    for (int i = 1; i < n - 1; ++i) {
      const double sig = (x_[i] - x_[i - 1]) / (x_[i + 1] - x_[i - 1]);
      const double p = sig * m[i - 1] + 2.0;
      m[i] = (sig - 1.0) / p;
      const double d = (y_[f][i + 1] - y_[f][i]) / (x_[i + 1] - x_[i]) -
                       (y_[f][i] - y_[f][i - 1]) / (x_[i] - x_[i - 1]);
      u[i] = (6.0 * d / (x_[i + 1] - x_[i - 1]) - sig * u[i - 1]) / p;
    }
    m[n - 1] = 0.0;
    for (int k = n - 2; k >= 0; --k) m[k] = m[k] * m[k + 1] + u[k];
    m_[f] = m;
  }
}

double SmearedCoulombTable::eval(int f, double r) const {
  const double x = std::log(std::max(r, std::exp(x_.front())));
  if (x >= x_.back()) {
    // beyond the table: pure Coulomb tails; G ~ -(4/3) alpha_0 / r is not
    // assumed here, clamp the cubic extrapolation by linear extrapolation of
    // the last interval
    const int n = static_cast<int>(x_.size());
    const double h = x_[n - 1] - x_[n - 2];
    const double slope = (y_[f][n - 1] - y_[f][n - 2]) / h +
                         h * (m_[f][n - 2] + 2 * m_[f][n - 1]) / 6.0;
    return y_[f][n - 1] + slope * (x - x_[n - 1]);
  }
  int lo =
      static_cast<int>(std::upper_bound(x_.begin(), x_.end(), x) - x_.begin()) -
      1;
  lo = std::max(0, std::min(lo, static_cast<int>(x_.size()) - 2));
  const double h = x_[lo + 1] - x_[lo], A = (x_[lo + 1] - x) / h,
               B = (x - x_[lo]) / h;
  return A * y_[f][lo] + B * y_[f][lo + 1] +
         ((A * A * A - A) * m_[f][lo] + (B * B * B - B) * m_[f][lo + 1]) * h *
             h / 6.0;
}

namespace {
// Lawson-Hanson non-negative least squares: min |A x - b|, x >= 0.
Vec nnls(const Mat& A, const Vec& b) {
  const int n = static_cast<int>(A.cols());
  Vec x = Vec::Zero(n);
  std::vector<bool> P(n, false);
  for (int outer = 0; outer < 3 * n; ++outer) {
    const Vec w = A.transpose() * (b - A * x);
    int jmax = -1;
    double wmax = 1e-12 * (1.0 + b.norm());
    for (int j = 0; j < n; ++j)
      if (!P[j] && w(j) > wmax) {
        wmax = w(j);
        jmax = j;
      }
    if (jmax < 0) break;
    P[jmax] = true;
    for (int inner = 0; inner < 3 * n; ++inner) {
      std::vector<int> idx;
      for (int j = 0; j < n; ++j)
        if (P[j]) idx.push_back(j);
      Mat Ap(A.rows(), static_cast<int>(idx.size()));
      for (std::size_t k = 0; k < idx.size(); ++k)
        Ap.col(static_cast<int>(k)) = A.col(idx[k]);
      const Vec zp = Ap.colPivHouseholderQr().solve(b);
      Vec z = Vec::Zero(n);
      for (std::size_t k = 0; k < idx.size(); ++k)
        z(idx[k]) = zp(static_cast<int>(k));
      bool ok = true;
      for (int j : idx) ok = ok && z(j) > 0;
      if (ok) {
        x = z;
        break;
      }
      double alpha = 1.0;
      for (int j : idx)
        if (z(j) <= 0) alpha = std::min(alpha, x(j) / (x(j) - z(j)));
      x += alpha * (z - x);
      for (int j : idx)
        if (x(j) <= 1e-14) {
          x(j) = 0;
          P[j] = false;
        }
    }
  }
  return x;
}
}  // namespace

std::vector<double> default_shape_gamma() {
  std::vector<double> g;
  for (int k = -4; k <= 8; ++k) g.push_back(0.25 * std::pow(2.0, k));
  return g;
}

GaussianShape fit_gaussian_shape(const RunningCoupling& a,
                                 std::vector<double> gamma, double qmax) {
  const double a0 = a.alpha(0.0);
  if (!(std::isfinite(a0) && a0 > 0))
    throw std::invalid_argument(
        "fit_gaussian_shape: the coupling must be finite and positive at Q = 0 "
        "(a confining form such as Richardson's has no Gaussian-sum shape)");
  const int nq = 400, nk = static_cast<int>(gamma.size());
  Mat A(nq + 1, nk);
  Vec b(nq + 1);
  std::vector<double> Q(nq);
  for (int i = 0; i < nq; ++i)
    Q[i] = 0.02 * std::pow(qmax / 0.02, static_cast<double>(i) / (nq - 1));
  for (int i = 0; i < nq; ++i) {
    const double s = a.alpha(Q[i] * Q[i]) / a0;
    for (int k = 0; k < nk; ++k)
      A(i, k) = std::exp(-Q[i] * Q[i] / (4 * gamma[k] * gamma[k])) / s;
    b(i) = 1.0;  // relative residual: fit/shape - 1
  }
  const double big = 1e4;  // sum w = 1 (alpha(0) = alpha_s exactly)
  A.row(nq).setConstant(big);
  b(nq) = big;
  const Vec x = nnls(A, b);
  GaussianShape out;
  out.gamma = gamma;
  const double tot = x.sum();
  for (int k = 0; k < nk; ++k) out.w.push_back(x(k) / tot);
  for (int i = 0; i < nq; ++i) {
    if (Q[i] < 0.05) continue;  // r > 4 fm: linear confinement dominates
    double f = 0;
    for (int k = 0; k < nk; ++k)
      f += out.w[k] * std::exp(-Q[i] * Q[i] / (4 * gamma[k] * gamma[k]));
    out.max_rel_err = std::max(
        out.max_rel_err, std::fabs(f / (a.alpha(Q[i] * Q[i]) / a0) - 1.0));
  }
  return out;
}

}  // namespace qgem
