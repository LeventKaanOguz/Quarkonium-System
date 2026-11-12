#include "data/dataset.hpp"

#include <fstream>
#include <sstream>
#include <stdexcept>

namespace qgem {

std::vector<const Datum*> Dataset::fitted() const {
  std::vector<const Datum*> out;
  for (const Datum& x : d)
    if (x.in_fit) out.push_back(&x);
  return out;
}

int Dataset::n_fit() const { return static_cast<int>(fitted().size()); }

Mat Dataset::V_exp() const {
  const auto f = fitted();
  const int n = static_cast<int>(f.size());
  Mat V = Mat::Zero(n, n);
  for (int i = 0; i < n; ++i) {
    V(i, i) = f[i]->err * f[i]->err;
    for (int j = 0; j < n; ++j)
      if (!f[i]->sys_group.empty() && f[i]->sys_group == f[j]->sys_group)
        V(i, j) += f[i]->sys_common * f[j]->sys_common;
  }
  return V;
}

Dataset load_dataset(const std::string& path, const PdgDatabase& pdg) {
  std::ifstream in(path);
  if (!in) throw std::runtime_error("cannot open states file " + path);
  Dataset ds;
  std::string line;
  int ln = 0;
  while (std::getline(in, line)) {
    ++ln;
    if (const auto h = line.find('#'); h != std::string::npos) line.erase(h);
    std::istringstream ss(line);
    long mcid;
    std::string channel, mode, group;
    int n;
    if (!(ss >> mcid)) continue;
    if (!(ss >> channel >> n))
      throw std::runtime_error(
          path + ":" + std::to_string(ln) +
          ": expected 'mcid channel n [fit|predict] [group]'");
    ss >> mode >> group;
    const auto m = pdg.mass(mcid);
    if (!m)
      throw std::runtime_error(path + ":" + std::to_string(ln) + ": mcid " +
                               std::to_string(mcid) + " not in the PDG file");
    Datum x;
    x.label = m->name;
    x.mcid = mcid;
    x.state = {channel, n};
    x.y = m->value;
    x.err = m->err;
    x.in_fit = (mode != "predict");
    if (!group.empty() && m->common_err > 0) {
      x.sys_group = group;
      x.sys_common = m->common_err;
    }
    ds.d.push_back(x);
  }
  return ds;
}

}  // namespace qgem
