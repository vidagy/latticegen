#include <gtest/gtest.h>
#include <gmock/gmock.h>

#include <Math/Integrator.h>

using namespace Math;

namespace
{
  static const double tolerance = 8.0 * std::numeric_limits<double>::epsilon();

  void test_simpson_const_up_to(double f, size_t n)
  {
    for (auto i = 2u; i <= n; ++i) {
      auto const_n_f = std::vector<double>(i, f);
      EXPECT_NEAR(IntegratorEquidistant::simpson(const_n_f, 1.0 / (i - 1)), f, tolerance)
              << " with n= " << n << " i= " << i << " f= " << f;
    }
  }
}

TEST(TestIntegrator, SimpsonConstant)
{
  test_simpson_const_up_to(8.4, 15);

  auto const_1_100 = std::vector<double>(101, 1.0);
  EXPECT_NEAR(IntegratorEquidistant::simpson(const_1_100, 1.0 / 100.0), 1.0, tolerance);
}

TEST(TestIntegrator, SimpsonPolinom)
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

namespace
{
  void test_simpson_3_8_const_up_to(double f, size_t n)
  {
    for (auto i = 2u; i <= n; ++i) {
      auto const_n_f = std::vector<double>(i, f);
      EXPECT_NEAR(IntegratorEquidistant::simpson_3_8(const_n_f, 1.0 / (i - 1)), f, tolerance)
              << " with n= " << n << " i= " << i << " f= " << f;
    }
  }
}

TEST(TestIntegrator, Simpson38Constant)
{
  test_simpson_3_8_const_up_to(1.0, 15);

  auto const_1_100 = std::vector<double>(101, 1.0);
  EXPECT_NEAR(IntegratorEquidistant::simpson_3_8(const_1_100, 1.0 / 100.0), 1.0, tolerance);
}

TEST(TestIntegrator, Simpson38Polinom)
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

  EXPECT_NEAR(IntegratorEquidistant::simpson_3_8(polynomial, dx), P(b) - P(a), 8.0 * tolerance);
}