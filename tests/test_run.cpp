// The end-to-end run: default configuration prints every PDG state found and
// reproduces the worked fit.
#include <gtest/gtest.h>

#include <filesystem>
#include <sstream>

#include "report/run.hpp"

using namespace qgem;

namespace {
bool have_data() {
  return std::filesystem::exists(std::string(QGEM_SOURCE_DIR) +
                                 "/data/pdg-2026.0.sqlite");
}
std::string output(RunOptions o) {
  o.root = QGEM_SOURCE_DIR;
  std::ostringstream s;
  EXPECT_EQ(run(o, s), 0);
  return s.str();
}
}  // namespace

TEST(Run, DefaultConfigurationPrintsEveryStateAndReproducesTheWorkedFit) {
  if (!have_data()) GTEST_SKIP();
  RunOptions
      o;  // defaults: model A, constant alpha_s, c = 0, ccbar + bbar + bc
  const std::string t = output(o);
  // every PDG state in the three state lists appears
  for (const char* name :
       {"eta_c(1S)", "J/psi(1S)", "chi_c0(1P)", "psi(3770)", "chi_c2(3930)",
        "psi(4040)", "psi(4160)", "psi(4415)", "eta_b(1S)", "Upsilon(1S)",
        "chi_b2(3P)", "Upsilon(4S)", "Upsilon(10860)", "Upsilon(11020)",
        "B_c+"})
    EXPECT_NE(t.find(name), std::string::npos) << name;
  EXPECT_NE(t.find("c fixed to 0"), std::string::npos);
  EXPECT_NE(t.find("theory floor sigma_th = 8.344 MeV"),
            std::string::npos);  // notes 10.10: 8.34 MeV
  // widths at LO and NLO, compared with the PDG (also BR x total width)
  EXPECT_NE(t.find("one-loop QCD factor"), std::string::npos);
  EXPECT_NE(t.find("e+ e-"), std::string::npos);
  EXPECT_NE(t.find("+- 0.1123*"), std::string::npos);  // chi_c0 -> gamma gamma
  // the truncation test of the non-relativistic model
  EXPECT_NE(t.find("a-priori floor sigma_a"), std::string::npos);
}

TEST(Run, ModelRunningAndConstantChoices) {
  if (!have_data()) GTEST_SKIP();
  RunOptions o;
  o.systems = {"ccbar"};
  o.with_c = true;
  std::string t = output(o);
  EXPECT_NE(t.find("c FREE"), std::string::npos);
  EXPECT_NE(t.find("  c  "),
            std::string::npos);  // c appears among the free parameters
  o.with_c = false;
  o.model = "B";
  EXPECT_NE(output(o).find("model B"), std::string::npos);
  o.model = "A";
  o.running = "gi";
  EXPECT_NE(output(o).find("running alpha_s"), std::string::npos);
  for (const char* r : {"frozen", "apt"}) {  // finite running forms, A/B only
    o.running = r;
    const std::string t = output(o);
    EXPECT_NE(t.find(std::string(r) == "apt" ? "APT shape" : "frozen shape"),
              std::string::npos)
        << r;
  }
  o.running =
      "free";  // free amplitudes exist only for the Godfrey-Isgur models
  std::ostringstream s;
  EXPECT_THROW(run(o, s), std::runtime_error);
  o.model = "X";
  o.running = "const";
  EXPECT_THROW(run(o, s), std::runtime_error);
}

TEST(Run, SeedFromAnotherMinimiserIsKeptWhenBetter) {
  // Model A, APT coupling, smeared spin terms, bottomonium: from the default
  // start Levenberg-Marquardt runs the smearing sigma into its flat sigma -> 0
  // corner (sigma_th ~ 107 MeV); Minuit2 finds sigma = 3.54 (2.27 MeV). A seed
  // at the Minuit2 point must be fitted and kept.
  if (!have_data()) GTEST_SKIP();
  RunOptions o;
  o.systems = {"bbar"};
  o.running = "apt";
  o.spin = "smeared";
  o.starts["bbar"].push_back(
      {{"m1", 4.969}, {"alpha_s", 1.636}, {"b", 0.1385}, {"sigma", 3.5409}});
  o.starts["ccbar"].push_back({{"m1", 1.0}});  // other systems: ignored
  const std::string t = output(o);
  EXPECT_NE(t.find("1 seed,"), std::string::npos);
  EXPECT_NE(t.find("best = start 1"), std::string::npos);
  EXPECT_NE(t.find("theory floor sigma_th = 2.27"), std::string::npos);
}

TEST(Run, FreeFixAndLeaveMultipletOut) {
  if (!have_data()) GTEST_SKIP();
  RunOptions o;
  o.systems = {"ccbar"};
  o.fixed = {{"alpha_s", 0.5}};
  o.freed = {"c"};  // the model keeps c fixed at 0 unless --c on
  o.lmo = true;
  const std::string t = output(o);
  EXPECT_NE(t.find("\n  c  "), std::string::npos);        // c is free
  EXPECT_EQ(t.find("\n  alpha_s  "), std::string::npos);  // alpha_s is not
  EXPECT_NE(t.find("alpha_s=0.5"), std::string::npos);    // listed as fixed
  // every (n, L) multiplet of the charmonium fit is left out once
  EXPECT_NE(t.find("leave-multiplet-out"), std::string::npos);
  for (const char* g : {"  1S ", "  2S ", "  1P "})
    EXPECT_NE(t.find(g), std::string::npos) << g;
  o.freed = {"no_such_parameter"};
  std::ostringstream s;
  EXPECT_THROW(run(o, s), std::exception);
}
