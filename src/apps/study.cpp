// qgem_study: small reproducible studies that are not fits (they evaluate a
// model at given parameters). Every number of the sections of the paper that
// refer to it is printed by one subcommand.
//   qgem_study gi-levels [--sandwich closure|exact] [--n N --r1 R --rmax R]
//       Godfrey-Isgur Hamiltonian with the published parameters (BGS for
//       charmonium, Godfrey-Moats for bottomonium): the lowest levels (MeV) of
//       every channel up to L = 3.
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

#include "models/gi_model.hpp"
#include "models/nr_model.hpp"
#include "obs/observables.hpp"
#include "solve/spectrum.hpp"

using namespace qgem;

namespace {

int gi_levels(const std::vector<std::string>& args) {
  bool exact = false;
  BasisSpec spec{.n = 30, .r1 = 0.02, .rmax = 30.0, .eta = 1e-12};
  int nlev = 8;
  double exact_rmax = 40.0;
  for (std::size_t i = 0; i < args.size(); ++i) {
    const std::string& k = args[i];
    auto val = [&]() -> const std::string& {
      if (++i >= args.size()) throw std::runtime_error(k + ": value missing");
      return args[i];
    };
    if (k == "--sandwich")
      exact = val() == "exact";
    else if (k == "--n")
      spec.n = std::stoi(val());
    else if (k == "--r1")
      spec.r1 = std::stod(val());
    else if (k == "--rmax")
      spec.rmax = std::stod(val());
    else if (k == "--levels")
      nlev = std::stoi(val());
    else if (k == "--exact-rmax")
      exact_rmax = std::stod(val());
    else
      throw std::runtime_error("unknown option " + k);
  }
  GIModel model({.exact_sandwich = exact, .exact_rmax = exact_rmax});
  std::cout << std::fixed << std::setprecision(2);
  for (double m : {1.628, 4.977}) {
    std::cout << std::defaultfloat << "m=" << m << std::fixed << "\n";
    const Vec th = model.godfrey_isgur_params(m, m).theta0();
    MatrixCache cache(spec);
    for (const Channel& ch : model.channels(3, true)) {
      if (ch.J > 4) continue;
      const SolvedChannel sc = solve_channel(model, ch, th, cache);
      std::cout << std::left << std::setw(10) << ch.label() << std::right;
      for (int k = 0; k < nlev && k < static_cast<int>(sc.levels.size()); ++k)
        std::cout << std::setw(10) << sc.level(k).M * 1000.0;
      std::cout << "\n";
    }
  }
  return 0;
}

// qgem_study contact-r0 --m M --alpha A --b B --sigma S
//   Model A (constant alpha_s): |R(0)|^2 of the 1S and 2S vector states with
//   the smeared contact hyperfine term (a) absent from H (SpinMode::None: the
//   Radford-Repko "H0") and (b) present in H (SpinMode::FirstOrder), at the
//   same parameters. Answers whether our nearly spin-independent |R(0)|^2 is a
//   consequence of keeping the contact term in H at all orders.
int contact_r0(const std::vector<std::string>& args) {
  double m = 1.467, a = 0.534, b = 0.1485, c = 0.0, sg = 1.130;
  BasisSpec spec{.n = 36, .r1 = 0.02, .rmax = 60.0, .eta = 1e-12};
  for (std::size_t i = 0; i + 1 < args.size(); i += 2) {
    const std::string& k = args[i];
    const double v = std::stod(args[i + 1]);
    if (k == "--m")
      m = v;
    else if (k == "--alpha")
      a = v;
    else if (k == "--b")
      b = v;
    else if (k == "--sigma")
      sg = v;
    else if (k == "--c")
      c = v;
    else
      throw std::runtime_error("unknown option " + k);
  }
  std::cout << std::fixed << std::setprecision(4);
  for (SpinMode sm : {SpinMode::None, SpinMode::FirstOrder}) {
    NRModel model({.spin = sm});
    Vec th(kNRParams);
    th << m, m, a, b, c, sg, 0.1;
    MatrixCache cache(spec);
    std::cout << (sm == SpinMode::None ? "H0 (no contact)   "
                                       : "contact in H      ");
    for (const char* lab : {"1S0", "3S1"}) {
      Channel ch;
      for (const Channel& cc : model.channels(0, true))
        if (cc.label() == lab) ch = cc;
      const SolvedChannel sc = solve_channel(model, ch, th, cache);
      for (int k = 0; k < 2; ++k)
        std::cout << lab << "(" << k + 1 << "S) M=" << std::setw(7)
                  << sc.level(k).M << " R0^2=" << std::setw(7)
                  << R0_squared(sc, k, cache) << "  ";
    }
    std::cout << "\n";
  }
  return 0;
}

// qgem_study complex-levels [--levels N]
//   Model A (constant coupling, first-order spin) at the fitted c c-bar
//   parameters of the paper: the lowest N levels of the 3S1 and 3P2 channels on
//   real-range bases and on mixed real + complex-range bases of the same total
//   number of functions (each complex range is two real functions, Re and Im of
//   r^l exp(-eta (1 + i w) r^2)), and the deviation of each from the lowest
//   value of each level over a large real basis (n = 160) and a large mixed
//   one. Masses in MeV.
int complex_levels(const std::vector<std::string>& args) {
  int nlev = 12;
  for (std::size_t i = 0; i + 1 < args.size(); i += 2) {
    if (args[i] == "--levels")
      nlev = std::stoi(args[i + 1]);
    else
      throw std::runtime_error("unknown option " + args[i]);
  }
  NRModel model({.spin = SpinMode::FirstOrder});
  Vec th(kNRParams);
  th << 1.46700880384, 1.46700880384, 0.533952502333, 0.148548883611, 0.0,
      1.13013331312, 0.1;
  auto levels = [&](const BasisSpec& spec, const std::string& lab) {
    MatrixCache cache(spec);
    Channel c32;
    std::vector<double> out;
    for (const char* name : {"3S1", "3P2"}) {
      Channel ch;
      for (const Channel& cc : model.channels(1, true))
        if (cc.label() == name) ch = cc;
      const SolvedChannel sc = solve_channel(model, ch, th, cache);
      for (int k = 0; k < nlev; ++k)
        out.push_back(k < static_cast<int>(sc.levels.size())
                          ? 1e3 * sc.level(k).M
                          : std::nan(""));
    }
    return out;
  };
  // every finite basis is an upper bound, so the reference is the lowest value
  // of each level over a large real basis and a large mixed one
  auto ref =
      levels({.n = 160, .r1 = 0.005, .rmax = 120.0, .eta = 1e-13}, "ref");
  const auto ref2 = levels({.n = 70,
                            .r1 = 0.01,
                            .rmax = 60.0,
                            .eta = 1e-13,
                            .nc = 40,
                            .rho1 = 0.5,
                            .rhomax = 60.0,
                            .omega = 0.5},
                           "ref2");
  for (std::size_t i = 0; i < ref.size(); ++i)
    ref[i] = std::min(ref[i], ref2[i]);
  struct Cfg {
    const char* name;
    BasisSpec spec;
  };
  const Cfg cfgs[] = {
      {"real 24", {.n = 24, .r1 = 0.05, .rmax = 25.0}},
      {"real 12 + cx 6 (24)",
       {.n = 12,
        .r1 = 0.05,
        .rmax = 25.0,
        .nc = 6,
        .rho1 = 0.8,
        .rhomax = 25.0,
        .omega = 0.5}},
      {"real 36", {.n = 36, .r1 = 0.05, .rmax = 25.0}},
      {"real 18 + cx 9 (36)",
       {.n = 18,
        .r1 = 0.05,
        .rmax = 25.0,
        .nc = 9,
        .rho1 = 0.8,
        .rhomax = 25.0,
        .omega = 0.5}},
      {"real 48", {.n = 48, .r1 = 0.05, .rmax = 25.0}},
      {"real 24 + cx 12 (48)",
       {.n = 24,
        .r1 = 0.05,
        .rmax = 25.0,
        .nc = 12,
        .rho1 = 0.8,
        .rhomax = 25.0,
        .omega = 0.5}},
  };
  std::cout
      << std::fixed << std::setprecision(3)
      << "reference (lowest of real n = 160 and mixed 70 + 40): masses (MeV)\n";
  for (int ch = 0; ch < 2; ++ch) {
    std::cout << (ch ? "3P2 " : "3S1 ");
    for (int k = 0; k < nlev; ++k)
      std::cout << std::setw(10) << ref[ch * nlev + k];
    std::cout << "\n";
  }
  std::cout << "deviation from the reference (MeV), basis size = number of "
               "real functions\n";
  for (const Cfg& c : cfgs) {
    const auto v = levels(c.spec, c.name);
    for (int ch = 0; ch < 2; ++ch) {
      std::cout << std::left << std::setw(22) << c.name
                << (ch ? "3P2 " : "3S1 ") << std::right;
      for (int k = 0; k < nlev; ++k)
        std::cout << std::setw(10) << v[ch * nlev + k] - ref[ch * nlev + k];
      std::cout << "\n";
    }
  }
  return 0;
}

}  // namespace

int main(int argc, char** argv) {
  try {
    if (argc < 2) {
      std::cerr << "usage: qgem_study gi-levels [options]\n";
      return 2;
    }
    const std::string cmd = argv[1];
    std::vector<std::string> args(argv + 2, argv + argc);
    if (cmd == "gi-levels") return gi_levels(args);
    if (cmd == "complex-levels") return complex_levels(args);
    if (cmd == "contact-r0") return contact_r0(args);
    std::cerr << "unknown study " << cmd << "\n";
    return 2;
  } catch (const std::exception& e) {
    std::cerr << "qgem_study: " << e.what() << "\n";
    return 1;
  }
}
