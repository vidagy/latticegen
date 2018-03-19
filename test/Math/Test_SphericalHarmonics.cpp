#include <TestUtils/base.h>

#include <Math/SphericalHarmonics.h>

using namespace Math;

namespace
{
  template<class T>
  void compare(T lhs, T rhs, double tol, unsigned int l, int m);

  template<>
  void compare<double>(double lhs, double rhs, double tol, unsigned int l, int m)
  {
    EXPECT_NEAR(lhs, rhs, tol) << "l: " << l << " m: " << m;
  }

  template<>
  void
  compare<std::complex<double>>(std::complex<double> lhs, std::complex<double> rhs, double tol, unsigned int l, int m)
  {
    EXPECT_NEAR(lhs.real(), rhs.real(), tol) << "l: " << l << " m: " << m;
    EXPECT_NEAR(lhs.imag(), rhs.imag(), tol) << "l: " << l << " m: " << m;
  }

  template<class T>
  void compare_equals(
    const std::function<T(const Point3D &)> &lhs,
    const std::function<T(const Point3D &)> &rhs, unsigned int l, int m)
  {
    for (double x = -1.95; x < 2.0; x += 0.5) {
      for (double y = -1.95; y < 2.0; y += 0.5) {
        for (double z = -1.95; z < 2.0; z += 0.5) {
          T lhs_res = lhs({x, y, z});
          T rhs_res = rhs({x, y, z});

          compare<T>(lhs_res, rhs_res, 32 * std::numeric_limits<double>::epsilon(), l, m);
        }
      }
    }
  }

  void compare_fast_vs_slow_real(unsigned int l, int m)
  {
    using std::placeholders::_1;
    auto fast = std::bind(Real::spherical_harmonic, l, m, _1);
    auto slow = std::bind(Real::spherical_harmonic_slow, l, m, _1);

    compare_equals<double>(fast, slow, l, m);
  }

  void compare_fast_vs_slow_complex(unsigned int l, int m)
  {
    using std::placeholders::_1;
    auto fast = std::bind(Complex::spherical_harmonic, l, m, _1);
    auto slow = std::bind(Complex::spherical_harmonic_slow, l, m, _1);

    compare_equals<std::complex<double>>(fast, slow, l, m);
  }
}

TEST(TestSphericalHarmonic, CompareRealL0)
{
  compare_fast_vs_slow_real(0, 0);
}

TEST(TestSphericalHarmonic, CompareRealL1)
{
  compare_fast_vs_slow_real(1, 1);
  compare_fast_vs_slow_real(1, 0);
  compare_fast_vs_slow_real(1, -1);
}

TEST(TestSphericalHarmonic, CompareRealL2)
{
  compare_fast_vs_slow_real(2, 2);
  compare_fast_vs_slow_real(2, 1);
  compare_fast_vs_slow_real(2, 0);
  compare_fast_vs_slow_real(2, -1);
  compare_fast_vs_slow_real(2, -2);
}

TEST(TestSphericalHarmonic, CompareRealL3)
{
  compare_fast_vs_slow_real(3, 3);
  compare_fast_vs_slow_real(3, 2);
  compare_fast_vs_slow_real(3, 1);
  compare_fast_vs_slow_real(3, 0);
  compare_fast_vs_slow_real(3, -1);
  compare_fast_vs_slow_real(3, -2);
  compare_fast_vs_slow_real(3, -3);
}

TEST(TestSphericalHarmonic, CompareRealL4)
{
  compare_fast_vs_slow_real(4, 4);
  compare_fast_vs_slow_real(4, 3);
  compare_fast_vs_slow_real(4, 2);
  compare_fast_vs_slow_real(4, 1);
  compare_fast_vs_slow_real(4, 0);
  compare_fast_vs_slow_real(4, -1);
  compare_fast_vs_slow_real(4, -2);
  compare_fast_vs_slow_real(4, -3);
  compare_fast_vs_slow_real(4, -4);
}


TEST(TestSphericalHarmonic, CompareComplexL0)
{
  compare_fast_vs_slow_complex(0, 0);
}

TEST(TestSphericalHarmonic, CompareComplexL1)
{
  compare_fast_vs_slow_complex(1, 1);
  compare_fast_vs_slow_complex(1, 0);
  compare_fast_vs_slow_complex(1, -1);
}

TEST(TestSphericalHarmonic, CompareComplexL2)
{
  compare_fast_vs_slow_complex(2, 2);
  compare_fast_vs_slow_complex(2, 1);
  compare_fast_vs_slow_complex(2, 0);
  compare_fast_vs_slow_complex(2, -1);
  compare_fast_vs_slow_complex(2, -2);
}

TEST(TestSphericalHarmonic, CompareComplexL3)
{
  compare_fast_vs_slow_complex(3, 3);
  compare_fast_vs_slow_complex(3, 2);
  compare_fast_vs_slow_complex(3, 1);
  compare_fast_vs_slow_complex(3, 0);
  compare_fast_vs_slow_complex(3, -1);
  compare_fast_vs_slow_complex(3, -2);
  compare_fast_vs_slow_complex(3, -3);
}

TEST(TestSphericalHarmonic, CompareComplexL4)
{
  compare_fast_vs_slow_complex(4, 4);
  compare_fast_vs_slow_complex(4, 3);
  compare_fast_vs_slow_complex(4, 2);
  compare_fast_vs_slow_complex(4, 1);
  compare_fast_vs_slow_complex(4, 0);
  compare_fast_vs_slow_complex(4, -1);
  compare_fast_vs_slow_complex(4, -2);
  compare_fast_vs_slow_complex(4, -3);
  compare_fast_vs_slow_complex(4, -4);
}
