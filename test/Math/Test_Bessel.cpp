#include <TestUtils/base.h>

#include <Math/Bessel.h>
#include <boost/math/special_functions/bessel.hpp>
#include <boost/math/special_functions/hankel.hpp>

using namespace Math;
using namespace std::complex_literals;

namespace
{
  void compare_bessel_to_boost(unsigned int n)
  {
    for (double x = 0.1; x < 10.0; x += 0.005) {
      double explicit_res = Math::bessel<double>(n, x);
      double boost_res = boost::math::sph_bessel(n, x);
      double slow = Math::bessel_slow<double>(n, x);
      EXPECT_NEAR(explicit_res, boost_res, 3e-12)
              << "for n = " << n << " x = " << x;
      EXPECT_NEAR(slow, boost_res, 3e-12)
              << "for n = " << n << " x = " << x;
      // std::cout << std::setw(20) << std::fixed << std::setprecision(15) << std::scientific
      //   << x << '\t'
      //   << std::abs(slow - boost_res) << '\t'
      //   << std::abs(explicit_res - boost_res) << std::endl;
    }
  }
}

TEST(TestBessel, CompareToBoost)
{
  compare_bessel_to_boost(0);
  compare_bessel_to_boost(1);
  compare_bessel_to_boost(2);
  compare_bessel_to_boost(3);
}

TEST(TestBessel, Complex)
{
  auto explicit_res = Math::bessel<std::complex<double>>(3, 5.0 + 1i);
  auto slow = Math::bessel_slow<std::complex<double>>(3, 5.0 + 1i);
  auto reference = 0.283003548719885467 - 0.0520917970290797114i;
  EXPECT_NEAR(explicit_res.real(), reference.real(), 1e-15);
  EXPECT_NEAR(explicit_res.imag(), reference.imag(), 1e-15);
  EXPECT_NEAR(slow.real(), reference.real(), 1e-15);
  EXPECT_NEAR(slow.imag(), reference.imag(), 1e-15);
}

namespace
{
  void compare_neumann_to_boost(unsigned int n)
  {
    for (double x = 0.1; x < 10.0; x += 0.005) {
      double explicit_res = Math::neumann<double>(n, x);
      double boost_res = boost::math::sph_neumann(n, x);
      double slow = Math::neumann_slow<double>(n, x);
      EXPECT_NEAR(explicit_res, boost_res, 3e-11)
              << "for n = " << n << " x = " << x;
      EXPECT_NEAR(slow, boost_res, 3e-11)
              << "for n = " << n << " x = " << x;
    }
  }
}

TEST(TestNeumann, CompareToBoost)
{
  compare_neumann_to_boost(0);
  compare_neumann_to_boost(1);
  compare_neumann_to_boost(2);
  compare_neumann_to_boost(3);
}

TEST(TestNeumann, Complex)
{
  auto explicit_res = Math::neumann<std::complex<double>>(3, 5.0 + 1i);
  auto slow = Math::neumann_slow<std::complex<double>>(3, 5.0 + 1i);
  auto reference = 0.014791676815047422 + 0.18676911997446946i;
  EXPECT_NEAR(explicit_res.real(), reference.real(), 1e-15);
  EXPECT_NEAR(explicit_res.imag(), reference.imag(), 1e-15);
  EXPECT_NEAR(slow.real(), reference.real(), 1e-15);
  EXPECT_NEAR(slow.imag(), reference.imag(), 1e-15);
}

namespace
{
  void compare_hankel_1_to_boost(unsigned int n)
  {
    for (double x = 0.1; x < 10.0; x += 0.005) {
      std::complex<double> explicit_res = Math::hankel_1(n, x);
      auto boost_res = boost::math::sph_hankel_1(n, x);
      EXPECT_NEAR(explicit_res.real(), boost_res.real(), 3e-11)
              << "real part for n = " << n << " x = " << x;
      EXPECT_NEAR(explicit_res.imag(), boost_res.imag(), 3e-11)
              << "imag part for n = " << n << " x = " << x;
    }
  }

  void compare_hankel_2_to_boost(unsigned int n)
  {
    for (double x = 0.1; x < 10.0; x += 0.005) {
      std::complex<double> explicit_res = Math::hankel_2(n, x);
      auto boost_res = boost::math::sph_hankel_2(n, x);
      EXPECT_NEAR(explicit_res.real(), boost_res.real(), 3e-11)
              << "real part for n = " << n << " x = " << x;
      EXPECT_NEAR(explicit_res.imag(), boost_res.imag(), 3e-11)
              << "imag part for n = " << n << " x = " << x;
    }
  }
}

TEST(TestHankel, CompareToBoost)
{
  compare_hankel_1_to_boost(0);
  compare_hankel_1_to_boost(1);
  compare_hankel_1_to_boost(2);
  compare_hankel_1_to_boost(3);

  compare_hankel_2_to_boost(0);
  compare_hankel_2_to_boost(1);
  compare_hankel_2_to_boost(2);
  compare_hankel_2_to_boost(3);
}


