#include <Math/SphericalHarmonics.h>
#include <Math/Legendre.h>

using namespace Math::Real;

double Math::Real::spherical_harmonic_slow(unsigned int l, int m, const Point3D &v) {
  double cos_theta = v(2) / v.norm();
  double phi = atan2(v(1), v(0));

  if (m < 0) {
    int abs_m = abs(m);
    return sqrt(2) * sqrt((2 * l + 1) / (4.0 * pi) * factorial(l - abs_m) / factorial(l + abs_m)) *
           legendre_polynomial_slow(l, abs_m, cos_theta) * sin(abs_m * phi);
  } else if (m == 0) {
    return sqrt((2 * l + 1) / (4.0 * pi)) * legendre_polynomial_slow(l, 0, cos_theta);
  } else // m > 0
  {
    return sqrt(2) * sqrt((2 * l + 1) / (4.0 * pi) * factorial(l - m) / factorial(l + m)) *
           legendre_polynomial_slow(l, m, cos_theta) * cos(m * phi);
  }
}

double Math::Real::spherical_harmonic(unsigned int l, int m, const Point3D &v)
{
  if ((unsigned int) abs(m) > l)
    THROW_INVALID_ARGUMENT("m is not valid: " + std::to_string(m) + " while l is " + std::to_string(l));

  switch (l) {
    case 0:
      return spherical_harmonic_00(v);
    case 1:
      switch (m) {
        case -1:
          return spherical_harmonic_1m1(v);
        case 0:
          return spherical_harmonic_10(v);
        case 1:
          return spherical_harmonic_11(v);
      }
    case 2:
      switch (m) {
        case -2:
          return spherical_harmonic_2m2(v);
        case -1:
          return spherical_harmonic_2m1(v);
        case 0:
          return spherical_harmonic_20(v);
        case 1:
          return spherical_harmonic_21(v);
        case 2:
          return spherical_harmonic_22(v);
      }
    case 3:
      switch (m) {
        case -3:
          return spherical_harmonic_3m3(v);
        case -2:
          return spherical_harmonic_3m2(v);
        case -1:
          return spherical_harmonic_3m1(v);
        case 0:
          return spherical_harmonic_30(v);
        case 1:
          return spherical_harmonic_31(v);
        case 2:
          return spherical_harmonic_32(v);
        case 3:
          return spherical_harmonic_33(v);
      }
    case 4:
      switch (m) {
        case -4:
          return spherical_harmonic_4m4(v);
        case -3:
          return spherical_harmonic_4m3(v);
        case -2:
          return spherical_harmonic_4m2(v);
        case -1:
          return spherical_harmonic_4m1(v);
        case 0:
          return spherical_harmonic_40(v);
        case 1:
          return spherical_harmonic_41(v);
        case 2:
          return spherical_harmonic_42(v);
        case 3:
          return spherical_harmonic_43(v);
        case 4:
          return spherical_harmonic_44(v);
      }
  }

  return spherical_harmonic_slow(l, m, v);
}

std::complex<double> Math::Complex::spherical_harmonic_slow(unsigned int l, int m, const Point3D &v)
{
  if (m == 0) {
    return std::complex<double>(Math::Real::spherical_harmonic_slow(l, m, v), 0.0);
  }

  auto abs_m = abs(m);
  auto Y_lpm = Math::Real::spherical_harmonic_slow(l, abs_m, v);
  auto Y_lmm = Math::Real::spherical_harmonic_slow(l, -abs_m, v);

  if (m < 0) {
    return 1.0 / sqrt(2.0) * std::complex<double>(Y_lpm, -Y_lmm);
  } else {
    return Math::sign(m) / sqrt(2.0) * std::complex<double>(Y_lpm, Y_lmm);
  }

}

std::complex<double> Math::Complex::spherical_harmonic(unsigned int l, int m, const Point3D &v)
{

  switch (l) {
    case 0:
      return Real::spherical_harmonic_00(v);
    case 1:
      switch (m) {
        case -1:
          return 1.0 / sqrt(2.0) * std::complex<double>(spherical_harmonic_11(v), -spherical_harmonic_1m1(v));
        case 0:
          return Real::spherical_harmonic_10(v);
        case 1:
          return -1.0 / sqrt(2.0) * std::complex<double>(spherical_harmonic_11(v), spherical_harmonic_1m1(v));
      }
    case 2:
      switch (m) {
        case -2:
          return 1.0 / sqrt(2.0) * std::complex<double>(spherical_harmonic_22(v), -spherical_harmonic_2m2(v));
        case -1:
          return 1.0 / sqrt(2.0) * std::complex<double>(spherical_harmonic_21(v), -spherical_harmonic_2m1(v));
        case 0:
          return Real::spherical_harmonic_20(v);
        case 1:
          return -1.0 / sqrt(2.0) * std::complex<double>(spherical_harmonic_21(v), spherical_harmonic_2m1(v));
        case 2:
          return 1.0 / sqrt(2.0) * std::complex<double>(spherical_harmonic_22(v), spherical_harmonic_2m2(v));
      }
    case 3:
      switch (m) {
        case -3:
          return 1.0 / sqrt(2.0) * std::complex<double>(spherical_harmonic_33(v), -spherical_harmonic_3m3(v));
        case -2:
          return 1.0 / sqrt(2.0) * std::complex<double>(spherical_harmonic_32(v), -spherical_harmonic_3m2(v));
        case -1:
          return 1.0 / sqrt(2.0) * std::complex<double>(spherical_harmonic_31(v), -spherical_harmonic_3m1(v));
        case 0:
          return Real::spherical_harmonic_30(v);
        case 1:
          return -1.0 / sqrt(2.0) * std::complex<double>(spherical_harmonic_31(v), spherical_harmonic_3m1(v));
        case 2:
          return 1.0 / sqrt(2.0) * std::complex<double>(spherical_harmonic_32(v), spherical_harmonic_3m2(v));
        case 3:
          return -1.0 / sqrt(2.0) * std::complex<double>(spherical_harmonic_33(v), spherical_harmonic_3m3(v));
      }
    case 4:
      switch (m) {
        case -4:
          return 1.0 / sqrt(2.0) * std::complex<double>(spherical_harmonic_44(v), -spherical_harmonic_4m4(v));
        case -3:
          return 1.0 / sqrt(2.0) * std::complex<double>(spherical_harmonic_43(v), -spherical_harmonic_4m3(v));
        case -2:
          return 1.0 / sqrt(2.0) * std::complex<double>(spherical_harmonic_42(v), -spherical_harmonic_4m2(v));
        case -1:
          return 1.0 / sqrt(2.0) * std::complex<double>(spherical_harmonic_41(v), -spherical_harmonic_4m1(v));
        case 0:
          return Real::spherical_harmonic_40(v);
        case 1:
          return -1.0 / sqrt(2.0) * std::complex<double>(spherical_harmonic_41(v), spherical_harmonic_4m1(v));
        case 2:
          return 1.0 / sqrt(2.0) * std::complex<double>(spherical_harmonic_42(v), spherical_harmonic_4m2(v));
        case 3:
          return -1.0 / sqrt(2.0) * std::complex<double>(spherical_harmonic_43(v), spherical_harmonic_4m3(v));
        case 4:
          return 1.0 / sqrt(2.0) * std::complex<double>(spherical_harmonic_44(v), spherical_harmonic_4m4(v));
      }
  }

  return Complex::spherical_harmonic_slow(l, m, v);
}