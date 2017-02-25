#include <gtest/gtest.h>
#include <gmock/gmock.h>

#include <Math/Integrator.h>

using namespace Math;

namespace
{
  static const double tolerance = std::numeric_limits<double>::epsilon();
}

TEST(TestIntegrator, Constant)
{
  auto const_1_1 = std::vector<double>(2, 1.0);
  EXPECT_NEAR(IntegratorEquidistant::simpson(const_1_1, 1.0 / 1.0), 1.0, tolerance);
  auto const_1_2 = std::vector<double>(3, 1.0);
  EXPECT_NEAR(IntegratorEquidistant::simpson(const_1_2, 1.0 / 2.0), 1.0, tolerance);
  auto const_1_3 = std::vector<double>(4, 1.0);
  EXPECT_NEAR(IntegratorEquidistant::simpson(const_1_3, 1.0 / 3.0), 1.0, tolerance);
  auto const_1_4 = std::vector<double>(5, 1.0);
  EXPECT_NEAR(IntegratorEquidistant::simpson(const_1_4, 1.0 / 4.0), 1.0, tolerance);
  auto const_1_5 = std::vector<double>(6, 1.0);
  EXPECT_NEAR(IntegratorEquidistant::simpson(const_1_5, 1.0 / 5.0), 1.0, tolerance);
  auto const_1_6 = std::vector<double>(7, 1.0);
  EXPECT_NEAR(IntegratorEquidistant::simpson(const_1_6, 1.0 / 6.0), 1.0, tolerance);

  auto const_1_100 = std::vector<double>(101, 1.0);
  EXPECT_NEAR(IntegratorEquidistant::simpson(const_1_100, 1.0 / 100.0), 1.0, tolerance);
}

TEST(TestIntegrator, Polinom)
{
  const int n = 201;
  auto p = [](double x)
  {
    return 3.0 * x * x - 5.0 * x + 2.0;
  };
  auto P = [](double x)
  {
    return x * x * x - 5.0 / 2.0 * x * x + 2.0 * x;
  };
  auto a = -2.0;
  auto b = 4.0;
  auto dx = (b - a) / (n - 1);
  auto x = [a, b, dx](unsigned int i)
  {
    return a + dx * i;
  };

  std::vector<double> polynomial;
  polynomial.reserve(n);
  for (unsigned int i = 0; i < n; ++i) {
    polynomial.push_back(p(x(i)));
  }

  EXPECT_NEAR(IntegratorEquidistant::simpson(polynomial, dx), P(b) - P(a), tolerance);
}