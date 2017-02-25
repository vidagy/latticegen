#include <gtest/gtest.h>
#include <gmock/gmock.h>

#include <Math/Legendre.h>

using namespace Math;

TEST(TestLegendre, L0M0)
{
  auto legendre_l0_m0 = legendre_polynomial(0, 0);

  EXPECT_DOUBLE_EQ(legendre_l0_m0(0.0), 1.0);
}

namespace
{
  void compare_equals(const std::function<double(double)> &lhs, const std::function<double(double)> &rhs)
  {
    for (double x = -0.95; x < 1.0; x += 0.05) {
      double lhs_res = lhs(x);
      double rhs_res = rhs(x);

      EXPECT_NEAR(lhs_res, rhs_res, 64 * std::numeric_limits<double>::epsilon());
    }
  }

  void compare_fast_vs_slow(int l, int m)
  {
    auto fast = legendre_polynomial(l, m);
    auto slow = [l, m](double x)
    {
      return legendre_polynomial_slow(l, m, x);
    };
    compare_equals(fast, slow);
  }
}

TEST(TestLegendre, NegativeL)
{
  auto negative = legendre_polynomial(-2, 1);
  auto positive = legendre_polynomial(1, 1);

  compare_equals(negative, positive);
}

TEST(TestLegendre, CompareL0)
{
  compare_fast_vs_slow(0, 0);
}

TEST(TestLegendre, CompareL1)
{
  compare_fast_vs_slow(1, 1);
  compare_fast_vs_slow(1, 0);
  compare_fast_vs_slow(1, -1);
}

TEST(TestLegendre, CompareL2)
{
  compare_fast_vs_slow(2, 2);
  compare_fast_vs_slow(2, 1);
  compare_fast_vs_slow(2, 0);
  compare_fast_vs_slow(2, -1);
  compare_fast_vs_slow(2, -2);
}

TEST(TestLegendre, CompareL3)
{
  compare_fast_vs_slow(3, 3);
  compare_fast_vs_slow(3, 2);
  compare_fast_vs_slow(3, 1);
  compare_fast_vs_slow(3, 0);
  compare_fast_vs_slow(3, -1);
  compare_fast_vs_slow(3, -2);
  compare_fast_vs_slow(3, -3);
}

TEST(TestLegendre, CompareL4)
{
  compare_fast_vs_slow(4, 4);
  compare_fast_vs_slow(4, 3);
  compare_fast_vs_slow(4, 2);
  compare_fast_vs_slow(4, 1);
  compare_fast_vs_slow(4, 0);
  compare_fast_vs_slow(4, -1);
  compare_fast_vs_slow(4, -2);
  compare_fast_vs_slow(4, -3);
  compare_fast_vs_slow(4, -4);
}

/*
 * For performance profiling
const int n = 1e7;

TEST(TestPerformance, inlined)
{
  double a = 0.0;
  for (double x = -1.0; x < 1.0; x += 2.0 / n)
  {
    a += legendre_polynomial_43(x);
  }
}

TEST(TestPerformance, intoAuto)
{
  auto f = &legendre_polynomial_43;
  double a = 0.0;
  for (double x = -1.0; x < 1.0; x += 2.0 / n)
  {
    a += f(x);
  }
}

TEST(TestPerformance, functionPointer)
{
  double (*f)(double) = legendre_polynomial_43;
  double a = 0.0;
  for (double x = -1.0; x < 1.0; x += 2.0 / n)
  {
    a += f(x);
  }
}

TEST(TestPerformance, slow)
{
  double a = 0.0;
  for (double x = -1.0; x < 1.0; x += 2.0 / n)
  {
    a += legendre_polynomial_slow(4,3,x);
  }
}

TEST(TestPerformance, trueInline)
{
  double a = 0.0;
  for (double x = -1.0; x < 1.0; x += 2.0 / n)
  {
    a += -105.0 * x * sqrt(1.0 - x * x) * (1.0 - x * x);
  }
}

TEST(TestPerformance, stdFunction)
{
  std::function<double(double)> f(legendre_polynomial_43);
  double a = 0.0;
  for (double x = -1.0; x < 1.0; x += 2.0 / n)
  {
    a += f(x);
  }
}
*/
