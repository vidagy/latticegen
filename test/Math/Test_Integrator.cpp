#include <gtest/gtest.h>
#include <gmock/gmock.h>

#include <Math/Integrator.h>
#include <Core/ExponentialMesh.h>

using namespace Core;
using namespace Math;

namespace
{
  static const double epsilon = std::numeric_limits<double>::epsilon();

  void test_const_up_to(double (*integrator)(const std::vector<double> &, double), double f, size_t n,
                        double tol_multiplier = 1.0)
  {
    for (auto i = 2u; i <= n; ++i) {
      auto const_n_f = std::vector<double>(i, f);
      EXPECT_NEAR(integrator(const_n_f, 1.0 / (i - 1)), f, tol_multiplier * epsilon)
              << " with n= " << n << " i= " << i << " f= " << f;
    }
  }
}

TEST(TestIntegratorEquidistant, SimpsonConstant)
{
  test_const_up_to(IntegratorEquidistant::simpson, 8.4, 15, 8);

  auto const_1_100 = std::vector<double>(501, 1.0);
  EXPECT_NEAR(IntegratorEquidistant::simpson(const_1_100, 1.0 / 500.0), 1.0, 4 * epsilon);
}

TEST(TestIntegratorEquidistant, Simpson38Constant)
{
  test_const_up_to(IntegratorEquidistant::simpson_3_8, 4.76, 15, 4);

  auto const_1_100 = std::vector<double>(501, 1.0);
  EXPECT_NEAR(IntegratorEquidistant::simpson_3_8(const_1_100, 1.0 / 500.0), 1.0, 64 * epsilon);
}

TEST(TestIntegratorEquidistant, SimpsonAltConstant)
{
  test_const_up_to(IntegratorEquidistant::simpson_alt, 1.73, 15);

  auto const_1_100 = std::vector<double>(501, 1.0);
  EXPECT_NEAR(IntegratorEquidistant::simpson_alt(const_1_100, 1.0 / 500.0), 1.0, 4 * epsilon);
}

namespace
{
  static const auto p = [](double x)
  {
    return 3.0 * x * x - 5.0 * x + 2.0;
  };
  static const auto P = [](double x)
  {
    return x * x * x - 5.0 / 2.0 * x * x + 2.0 * x;
  };

  constexpr double get_dx(int n, double a, double b)
  {
    return (b - a) / n;
  }

  std::vector<double> get_polynomial(unsigned int n = 201, double a = -2.0, double b = 4.0)
  {
    auto x = [a, b, n](unsigned int i)
    {
      return a + get_dx(n, a, b) * i;
    };

    std::vector<double> polynomial;
    polynomial.reserve(n + 1);
    for (unsigned int i = 0; i <= n; ++i) {
      polynomial.push_back(p(x(i)));
    }
    return polynomial;
  }
}

static const unsigned int n = 200;
static const double a = -2.0;
static const double b = 4.0;
static const double dx = get_dx(n, a, b);

TEST(TestIntegratorEquidistant, SimpsonPolinom)
{
  auto polynomial = get_polynomial(n, a, b);
  EXPECT_NEAR(IntegratorEquidistant::simpson(polynomial, dx), P(b) - P(a), epsilon);
}

TEST(TestIntegratorEquidistant, Simpson38Polinom)
{
  auto polynomial = get_polynomial(n, a, b);
  EXPECT_NEAR(IntegratorEquidistant::simpson_3_8(polynomial, dx), P(b) - P(a), 32 * epsilon);
}

TEST(TestIntegratorEquidistant, SimpsonAltPolinom)
{
  auto polynomial = get_polynomial(n, a, b);
  EXPECT_NEAR(IntegratorEquidistant::simpson_alt(polynomial, dx), P(b) - P(a), epsilon);
}

TEST(TestIntegratorGeneric, Constant)
{
  auto const_n_f = std::vector<double>(501, 1.0);
  auto points = ExponentialMesh(1.0, 0.1, 501).points;
  EXPECT_NEAR(IntegratorGeneric::integrate(const_n_f, points), 1.0, epsilon);
}

TEST(TestIntegratorGeneric, Exponential)
{
  auto n = 10000u;
  auto b = 20.0;

  std::vector<double> points;
  points.reserve(n);
  for (auto i = 0u; i < n; ++i) {
    points.push_back((b) / (n - 1) * i);
  }

  std::vector<double> polynomial;
  polynomial.reserve(points.size());
  for (auto x: points) {
    polynomial.push_back(exp(-x));
  }

  EXPECT_NEAR(IntegratorGeneric::integrate(polynomial, points), 1.0 - exp(-b), 4e-7);
}