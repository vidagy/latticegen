#include <gtest/gtest.h>
#include <gmock/gmock.h>

#include <Core/Point3D.h>
#include <Math/SphericalHarmonics.h>

using namespace Math;

namespace
{
  void compare_equals(
    const std::function<double(const Vector3D &)> &lhs,
    const std::function<double(const Vector3D &)> &rhs, unsigned int l, int m)
  {
    for (double x = -1.95; x < 2.0; x += 0.5) {
      for (double y = -1.95; y < 2.0; y += 0.5) {
        for (double z = -1.95; z < 2.0; z += 0.5) {
          double lhs_res = lhs({x, y, z});
          double rhs_res = rhs({x, y, z});

          EXPECT_NEAR(lhs_res, rhs_res, 32 * std::numeric_limits<double>::epsilon()) << "l: " << l << " m: " << m;
        }
      }
    }
  }

  void compare_fast_vs_slow(unsigned int l, int m)
  {
    auto fast = spherical_harmonic(l, m);
    auto slow = [l, m](const Vector3D &x)
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
