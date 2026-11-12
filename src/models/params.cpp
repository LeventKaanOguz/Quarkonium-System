#include "models/params.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace qgem {

ParamSet::ParamSet(std::vector<ParamSpec> specs) : specs_(std::move(specs)) {
  tie_src_.assign(specs_.size(), -1);
  for (std::size_t i = 0; i < specs_.size(); ++i) {
    if (!specs_[i].tied_to.empty()) {
      specs_[i].free = false;
      tie_src_[i] = index(specs_[i].tied_to);
      if (tie_src_[i] == static_cast<int>(i))
        throw std::invalid_argument("ParamSet: parameter tied to itself");
    }
    if (specs_[i].tr == Transform::Log && specs_[i].free &&
        !(specs_[i].value > 0))
      throw std::invalid_argument("ParamSet: log-transformed parameter '" +
                                  specs_[i].name + "' must be > 0");
    if (specs_[i].free) free_idx_.push_back(static_cast<int>(i));
  }
  for (std::size_t i = 0; i < specs_.size(); ++i)
    if (tie_src_[i] >= 0 && tie_src_[tie_src_[i]] >= 0)
      throw std::invalid_argument("ParamSet: chained ties are not supported");
}

int ParamSet::index(const std::string& name) const {
  for (std::size_t i = 0; i < specs_.size(); ++i)
    if (specs_[i].name == name) return static_cast<int>(i);
  throw std::invalid_argument("ParamSet: unknown parameter '" + name + "'");
}

std::vector<std::string> ParamSet::free_names() const {
  std::vector<std::string> out;
  for (int i : free_idx_) out.push_back(specs_[i].name);
  return out;
}

Vec ParamSet::theta0() const {
  Vec t(n_total());
  for (int i = 0; i < n_total(); ++i) t(i) = specs_[i].value;
  for (int i = 0; i < n_total(); ++i)
    if (tie_src_[i] >= 0) t(i) = t(tie_src_[i]);
  return t;
}

Vec ParamSet::u0() const { return u(theta0()); }

Vec ParamSet::u(const Vec& th) const {
  Vec out(n_free());
  for (int k = 0; k < n_free(); ++k) {
    const ParamSpec& s = specs_[free_idx_[k]];
    out(k) = (s.tr == Transform::Log) ? std::log(th(free_idx_[k]))
                                      : th(free_idx_[k]);
  }
  return out;
}

Vec ParamSet::theta(const Vec& u) const {
  Vec t = theta0();
  for (int k = 0; k < n_free(); ++k) {
    const ParamSpec& s = specs_[free_idx_[k]];
    t(free_idx_[k]) = (s.tr == Transform::Log) ? std::exp(u(k)) : u(k);
  }
  for (int i = 0; i < n_total(); ++i)
    if (tie_src_[i] >= 0) t(i) = t(tie_src_[i]);
  return t;
}

Mat ParamSet::dtheta_du(const Vec& u) const {
  Mat K = Mat::Zero(n_total(), n_free());
  const Vec th = theta(u);
  for (int k = 0; k < n_free(); ++k) {
    const int i = free_idx_[k];
    K(i, k) = (specs_[i].tr == Transform::Log) ? th(i) : 1.0;
  }
  for (int i = 0; i < n_total(); ++i)
    if (tie_src_[i] >= 0) K.row(i) = K.row(tie_src_[i]);
  return K;
}

bool ParamSet::within_bounds(const Vec& th) const {
  for (int i = 0; i < n_total(); ++i)
    if (!(th(i) >= specs_[i].lo && th(i) <= specs_[i].hi)) return false;
  return true;
}

ParamSet ParamSet::with_values(const Vec& th) const {
  auto s = specs_;
  for (int i = 0; i < n_total(); ++i) s[i].value = th(i);
  return ParamSet(std::move(s));
}

ParamSet ParamSet::with_free(const std::vector<std::string>& names) const {
  auto s = specs_;
  for (auto& p : s)
    if (p.tied_to.empty())
      p.free = std::find(names.begin(), names.end(), p.name) != names.end();
  return ParamSet(std::move(s));
}

}  // namespace qgem
