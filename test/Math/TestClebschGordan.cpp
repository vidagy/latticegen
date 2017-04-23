#include <TestUtils/base.h>

#include <Math/ClebschGordan.h>

using namespace Math;

namespace
{
  static const double tol = std::numeric_limits<double>::epsilon();
}

TEST(TestClebschGordan, Zero)
{
  // m out of bounds
  EXPECT_EQ(ClebschGordan::calculate(0, 1, 0, 0, 0, 0), 0.0);
  EXPECT_EQ(ClebschGordan::calculate(0, -1, 0, 0, 0, 0), 0.0);
  EXPECT_EQ(ClebschGordan::calculate(0, 0, 0, 1, 0, 0), 0.0);
  EXPECT_EQ(ClebschGordan::calculate(0, 0, 0, -1, 0, 0), 0.0);
  EXPECT_EQ(ClebschGordan::calculate(0, 0, 0, 0, 0, 1), 0.0);
  EXPECT_EQ(ClebschGordan::calculate(0, 0, 0, 0, 0, -1), 0.0);

  // L out of bounds
  EXPECT_EQ(ClebschGordan::calculate(-1, 0, 0, 0, 0, 0), 0.0);
  EXPECT_EQ(ClebschGordan::calculate(0, 0, -1, 0, 0, 0), 0.0);
  EXPECT_EQ(ClebschGordan::calculate(0, 0, 0, 0, -1, 0), 0.0);
}

TEST(TestClebschGordan, j0)
{
  EXPECT_EQ(ClebschGordan::calculate(0, 0, 0, 0, 0, 0), 1.0);

  auto j1 = 0.5;
  EXPECT_NEAR(ClebschGordan::calculate(j1, 0.5, 0, 0, j1, 0.5), 1.0, tol);
  EXPECT_NEAR(ClebschGordan::calculate(j1, 0.5, 0, 0, j1, -0.5), 0.0, tol);
  EXPECT_NEAR(ClebschGordan::calculate(j1, -0.5, 0, 0, j1, 0.5), 0.0, tol);
  EXPECT_NEAR(ClebschGordan::calculate(j1, -0.5, 0, 0, j1, -0.5), 1.0, tol);
  EXPECT_NEAR(ClebschGordan::calculate(0, 0, j1, 0.5, j1, 0.5), 1.0, tol);
  EXPECT_NEAR(ClebschGordan::calculate(0, 0, j1, 0.5, j1, -0.5), 0.0, tol);
  EXPECT_NEAR(ClebschGordan::calculate(0, 0, j1, -0.5, j1, 0.5), 0.0, tol);
  EXPECT_NEAR(ClebschGordan::calculate(0, 0, j1, -0.5, j1, -0.5), 1.0, tol);

  j1 = 1.0;
  EXPECT_NEAR(ClebschGordan::calculate(j1, 1, 0, 0, j1, 1), 1.0, tol);
  EXPECT_NEAR(ClebschGordan::calculate(j1, 1, 0, 0, j1, 0), 0.0, tol);
  EXPECT_NEAR(ClebschGordan::calculate(j1, 1, 0, 0, j1, -1), 0.0, tol);
  EXPECT_NEAR(ClebschGordan::calculate(j1, 0, 0, 0, j1, 1), 0.0, tol);
  EXPECT_NEAR(ClebschGordan::calculate(j1, 0, 0, 0, j1, 0), 1.0, tol);
  EXPECT_NEAR(ClebschGordan::calculate(j1, 0, 0, 0, j1, -1), 0.0, tol);
  EXPECT_NEAR(ClebschGordan::calculate(j1, -1, 0, 0, j1, 1), 0.0, tol);
  EXPECT_NEAR(ClebschGordan::calculate(j1, -1, 0, 0, j1, 0), 0.0, tol);
  EXPECT_NEAR(ClebschGordan::calculate(j1, -1, 0, 0, j1, -1), 1.0, tol);
  EXPECT_NEAR(ClebschGordan::calculate(0, 0, j1, 1, j1, 1), 1.0, tol);
  EXPECT_NEAR(ClebschGordan::calculate(0, 0, j1, 1, j1, 0), 0.0, tol);
  EXPECT_NEAR(ClebschGordan::calculate(0, 0, j1, 1, j1, -1), 0.0, tol);
  EXPECT_NEAR(ClebschGordan::calculate(0, 0, j1, 0, j1, 1), 0.0, tol);
  EXPECT_NEAR(ClebschGordan::calculate(0, 0, j1, 0, j1, 0), 1.0, tol);
  EXPECT_NEAR(ClebschGordan::calculate(0, 0, j1, 0, j1, -1), 0.0, tol);
  EXPECT_NEAR(ClebschGordan::calculate(0, 0, j1, -1, j1, 1), 0.0, tol);
  EXPECT_NEAR(ClebschGordan::calculate(0, 0, j1, -1, j1, 0), 0.0, tol);
  EXPECT_NEAR(ClebschGordan::calculate(0, 0, j1, -1, j1, -1), 1.0, tol);
}

TEST(TestClebschGordan, j0505)
{
  // L is off
  EXPECT_EQ(ClebschGordan::calculate(0.5, 0.5, 0.5, 0.5, -0.5, 1.0), 0.0);
  EXPECT_EQ(ClebschGordan::calculate(0.5, 0.5, 0.5, -0.5, 2.0, 0.0), 0.0);
  EXPECT_EQ(ClebschGordan::calculate(0.5, -0.5, 0.5, 0.5, 0.5, 0.0), 0.0);
  EXPECT_EQ(ClebschGordan::calculate(0.5, -0.5, 0.5, -0.5, -1.0, 1.0), 0.0);
  // M is off
  EXPECT_EQ(ClebschGordan::calculate(0.5, 0.5, 0.5, 0.5, 1.0, 0.0), 0.0);
  EXPECT_EQ(ClebschGordan::calculate(0.5, 0.5, 0.5, -0.5, 1.0, 3.0), 0.0);
  EXPECT_EQ(ClebschGordan::calculate(0.5, -0.5, 0.5, 0.5, 1.0, 2.0), 0.0);
  EXPECT_EQ(ClebschGordan::calculate(0.5, -0.5, 0.5, -0.5, 1.0, 0.0), 0.0);

  EXPECT_NEAR(ClebschGordan::calculate(0.5, 0.5, 0.5, 0.5, 1.0, 1.0), 1.0, tol);
  EXPECT_NEAR(ClebschGordan::calculate(0.5, -0.5, 0.5, -0.5, 1.0, -1.0), 1.0, tol);

  EXPECT_NEAR(ClebschGordan::calculate(0.5, -0.5, 0.5, 0.5, 1.0, 0.0), sqrt(1.0 / 2.0), tol);
  EXPECT_NEAR(ClebschGordan::calculate(0.5, -0.5, 0.5, 0.5, 0.0, 0.0), -sqrt(1.0 / 2.0), tol);
  EXPECT_NEAR(ClebschGordan::calculate(0.5, 0.5, 0.5, -0.5, 1.0, 0.0), sqrt(1.0 / 2.0), tol);
  EXPECT_NEAR(ClebschGordan::calculate(0.5, 0.5, 0.5, -0.5, 0.0, 0.0), sqrt(1.0 / 2.0), tol);
}

TEST(TestClebschGordan, j105)
{
  EXPECT_NEAR(ClebschGordan::calculate(0.5, 0.5, 1.0, 1.0, 0.5, 0.5), 0.0, tol);
  EXPECT_NEAR(ClebschGordan::calculate(0.5, 0.5, 1.0, 1.0, 0.5, -0.5), 0.0, tol);
  EXPECT_NEAR(ClebschGordan::calculate(0.5, -0.5, 1.0, 1.0, 0.5, 0.5), -sqrt(2.0 / 3.0), tol);
  EXPECT_NEAR(ClebschGordan::calculate(0.5, -0.5, 1.0, 1.0, 0.5, -0.5), 0.0, tol);
  EXPECT_NEAR(ClebschGordan::calculate(0.5, 0.5, 1.0, 0.0, 0.5, 0.5), sqrt(1.0 / 3.0), tol);
  EXPECT_NEAR(ClebschGordan::calculate(0.5, 0.5, 1.0, 0.0, 0.5, -0.5), 0.0, tol);
  EXPECT_NEAR(ClebschGordan::calculate(0.5, -0.5, 1.0, 0.0, 0.5, 0.5), 0.0, tol);
  EXPECT_NEAR(ClebschGordan::calculate(0.5, -0.5, 1.0, 0.0, 0.5, -0.5), -sqrt(1.0 / 3.0), tol);
  EXPECT_NEAR(ClebschGordan::calculate(0.5, 0.5, 1.0, -1.0, 0.5, 0.5), 0.0, tol);
  EXPECT_NEAR(ClebschGordan::calculate(0.5, 0.5, 1.0, -1.0, 0.5, -0.5), sqrt(2.0 / 3.0), tol);
  EXPECT_NEAR(ClebschGordan::calculate(0.5, -0.5, 1.0, -1.0, 0.5, 0.5), 0.0, tol);
  EXPECT_NEAR(ClebschGordan::calculate(0.5, -0.5, 1.0, -1.0, 0.5, -0.5), 0.0, tol);

  EXPECT_NEAR(ClebschGordan::calculate(0.5, 0.5, 1.0, 1.0, 1.5, 0.5), 0.0, tol);
  EXPECT_NEAR(ClebschGordan::calculate(0.5, 0.5, 1.0, 1.0, 1.5, -0.5), 0.0, tol);
  EXPECT_NEAR(ClebschGordan::calculate(0.5, -0.5, 1.0, 1.0, 1.5, 0.5), sqrt(1.0 / 3.0), tol);
  EXPECT_NEAR(ClebschGordan::calculate(0.5, -0.5, 1.0, 1.0, 1.5, -0.5), 0.0, tol);
  EXPECT_NEAR(ClebschGordan::calculate(0.5, 0.5, 1.0, 0.0, 1.5, 0.5), sqrt(2.0 / 3.0), tol);
  EXPECT_NEAR(ClebschGordan::calculate(0.5, 0.5, 1.0, 0.0, 1.5, -0.5), 0.0, tol);
  EXPECT_NEAR(ClebschGordan::calculate(0.5, -0.5, 1.0, 0.0, 1.5, 0.5), 0.0, tol);
  EXPECT_NEAR(ClebschGordan::calculate(0.5, -0.5, 1.0, 0.0, 1.5, -0.5), sqrt(2.0 / 3.0), tol);
  EXPECT_NEAR(ClebschGordan::calculate(0.5, 0.5, 1.0, -1.0, 1.5, 0.5), 0.0, tol);
  EXPECT_NEAR(ClebschGordan::calculate(0.5, 0.5, 1.0, -1.0, 1.5, -0.5), sqrt(1.0 / 3.0), tol);
  EXPECT_NEAR(ClebschGordan::calculate(0.5, -0.5, 1.0, -1.0, 1.5, 0.5), 0.0, tol);
  EXPECT_NEAR(ClebschGordan::calculate(0.5, -0.5, 1.0, -1.0, 1.5, -0.5), 0.0, tol);

  EXPECT_NEAR(ClebschGordan::calculate(0.5, 0.5, 1.0, 1.0, 1.5, 1.5), 1.0, tol);
  EXPECT_NEAR(ClebschGordan::calculate(0.5, -0.5, 1.0, 1.0, 1.5, 1.5), 0.0, tol);
  EXPECT_NEAR(ClebschGordan::calculate(0.5, 0.5, 1.0, 0.0, 1.5, 1.5), 0.0, tol);
  EXPECT_NEAR(ClebschGordan::calculate(0.5, -0.5, 1.0, 0.0, 1.5, 1.5), 0.0, tol);
  EXPECT_NEAR(ClebschGordan::calculate(0.5, 0.5, 1.0, -1.0, 1.5, 1.5), 0.0, tol);
  EXPECT_NEAR(ClebschGordan::calculate(0.5, -0.5, 1.0, -1.0, 1.5, 1.5), 0.0, tol);
}