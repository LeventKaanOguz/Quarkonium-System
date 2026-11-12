// Milestone M0: the toolchain works end to end (compiler, Eigen, GoogleTest,
// CMake).
#include <gtest/gtest.h>

#include "core/constants.hpp"
#include "core/types.hpp"

TEST(Toolchain, EigenSolvesASymmetricEigenproblem) {
  qgem::Mat a(2, 2);
  a << 2, 1, 1, 2;  // eigenvalues 1 and 3
  Eigen::SelfAdjointEigenSolver<qgem::Mat> es(a);
  ASSERT_EQ(es.info(), Eigen::Success);
  EXPECT_NEAR(es.eigenvalues()[0], 1.0, 1e-14);
  EXPECT_NEAR(es.eigenvalues()[1], 3.0, 1e-14);
}

TEST(Toolchain, UnitsAreSane) {
  EXPECT_NEAR(qgem::hbarc_GeV_fm, 0.19733, 1e-5);
  EXPECT_DOUBLE_EQ(1000.0 * qgem::MeV, 1.0);
}
