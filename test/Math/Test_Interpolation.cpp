#include <TestUtils/base.h>

#include <Math/CubicSplineInterpolation.h>

using namespace Core;
using namespace Math;

TEST(InterpolationCubicSpline, AllZero)
{
  const auto x = std::vector<double>{0.0, 1.0};
  const auto y = std::vector<double>{0.0, 0.0};
  auto spline = CubicSplineInterpolation<double, double>(x, y);

  for (auto xi = 0.0; xi <= 1.0; xi += 0.1)
    EXPECT_NEAR(spline.interpolate(xi), 0.0, std::numeric_limits<double>::epsilon());
}

TEST(InterpolationCubicSpline, Linear)
{
  const auto x = std::vector<double>{0.0, 1.0, 2.0, 3.0};
  const auto y = std::vector<double>{0.1, 0.2, 0.3, 0.4};
  auto spline = CubicSplineInterpolation<double, double>(x, y);

  for (auto xi = 0.0; xi <= 3.0; xi += 0.1)
    EXPECT_NEAR(spline.interpolate(xi), 0.1 * (xi + 1.0), std::numeric_limits<double>::epsilon());
}

TEST(InterpolationCubicSpline, LinearComplex)
{
  const auto x = std::vector<double>{0.0, 1.0, 2.0, 3.0};
  const auto y = std::vector<std::complex<double>>{{0.1, -1.1},
                                                   {0.2, -1.2},
                                                   {0.3, -1.3},
                                                   {0.4, -1.4}};
  auto spline = CubicSplineInterpolation<double, std::complex<double>>(x, y);

  for (auto xi = 0.0; xi <= 3.0; xi += 0.1) {
    EXPECT_NEAR(spline.interpolate(xi).real(), 0.1 * (xi + 1.0), std::numeric_limits<double>::epsilon());
    EXPECT_NEAR(spline.interpolate(xi).imag(), -0.1 * (xi + 1.0) - 1.0, std::numeric_limits<double>::epsilon());
  }
}
