#include <gtest/gtest.h>
#include <gmock/gmock.h>

#include <Core/Point3D.h>
#include <Core/SphericalHarmonics.h>

using namespace Core;

namespace
{
  void
  compare_equals(const std::function<double(double)> &lhs, const std::function<double(double)> &rhs, unsigned int l,
                 int m)
  {
    for (double x = -0.95; x < 1.0; x += 0.05) {
      double lhs_res = lhs(x);
      double rhs_res = rhs(x);

      EXPECT_NEAR(lhs_res, rhs_res, 32 * std::numeric_limits<double>::epsilon()) << "l: " << l << " m: " << m;
    }
  }

  void compare_fast_vs_slow(unsigned int l, int m)
  {
    auto fast = Core::spherical_harmonic(l, m);
    auto slow = [l, m](double x)
    {
      return spherical_harmonic_slow(l, m, x);
    };
    compare_equals(fast, slow, l, m);
  }
}

TEST(TestSphericalHarmonic, CompareL0)
{
  compare_fast_vs_slow(0, 0);
}

TEST(TestSphericalHarmonic, CompareL1)
{
  compare_fast_vs_slow(1, 1);
  compare_fast_vs_slow(1, 0);
  compare_fast_vs_slow(1, -1);
}

TEST(TestSphericalHarmonic, CompareL2)
{
  compare_fast_vs_slow(2, 2);
  compare_fast_vs_slow(2, 1);
  compare_fast_vs_slow(2, 0);
  compare_fast_vs_slow(2, -1);
  compare_fast_vs_slow(2, -2);
}

TEST(TestSphericalHarmonic, CompareL3)
{
  compare_fast_vs_slow(3, 3);
  compare_fast_vs_slow(3, 2);
  compare_fast_vs_slow(3, 1);
  compare_fast_vs_slow(3, 0);
  compare_fast_vs_slow(3, -1);
  compare_fast_vs_slow(3, -2);
  compare_fast_vs_slow(3, -3);
}

TEST(TestSphericalHarmonic, CompareL4)
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
