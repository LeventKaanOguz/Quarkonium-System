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
//   Model A (constant alpha_s): |R(0)|^2 of the 1S and 2S vector states with the
//   smeared contact hyperfine term (a) absent from H (SpinMode::None: the
//   Radford-Repko "H0") and (b) present in H (SpinMode::FirstOrder), at the same
//   parameters. Answers whether our nearly spin-independent |R(0)|^2 is a
//   consequence of keeping the contact term in H at all orders.
int contact_r0(const std::vector<std::string>& args) {
  double m = 1.467, a = 0.534, b = 0.1485, c = 0.0, sg = 1.130;
  BasisSpec spec{.n = 36, .r1 = 0.02, .rmax = 60.0, .eta = 1e-12};
  for (std::size_t i = 0; i + 1 < args.size(); i += 2) {
    const std::string& k = args[i];
    const double v = std::stod(args[i + 1]);
    if (k == "--m") m = v; else if (k == "--alpha") a = v;
    else if (k == "--b") b = v; else if (k == "--sigma") sg = v;
    else if (k == "--c") c = v; else throw std::runtime_error("unknown option " + k);
  }
  std::cout << std::fixed << std::setprecision(4);
  for (SpinMode sm : {SpinMode::None, SpinMode::FirstOrder}) {
    NRModel model({.spin = sm});
    Vec th(kNRParams); th << m, m, a, b, c, sg, 0.1;
    MatrixCache cache(spec);
    std::cout << (sm == SpinMode::None ? "H0 (no contact)   " : "contact in H      ");
    for (const char* lab : {"1S0", "3S1"}) {
      Channel ch;
      for (const Channel& cc : model.channels(0, true)) if (cc.label() == lab) ch = cc;
      const SolvedChannel sc = solve_channel(model, ch, th, cache);
      for (int k = 0; k < 2; ++k)
        std::cout << lab << "(" << k + 1 << "S) M=" << std::setw(7) << sc.level(k).M
                  << " R0^2=" << std::setw(7) << R0_squared(sc, k, cache) << "  ";
    }
    std::cout << "\n";
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
    if (cmd == "contact-r0") return contact_r0(args);
    std::cerr << "unknown study " << cmd << "\n";
    return 2;
  } catch (const std::exception& e) {
    std::cerr << "qgem_study: " << e.what() << "\n";
    return 1;
  }
}
