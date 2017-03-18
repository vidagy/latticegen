#include <Math/SphericalHarmonics.h>
#include <Math/Legendre.h>
#include <Core/Exceptions.h>

using namespace Math;

double Math::spherical_harmonic_slow(unsigned int l, int m, const Vector3D &v)
{
  double cos_theta = v.z / v.length();
  double phi = atan2(v.y, v.x);

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

std::function<double(const Core::Vector3D &)> Math::spherical_harmonic(unsigned int l, int m)
{
  if ((unsigned int) abs(m) > l)
    THROW_INVALID_ARGUMENT("m is not valid: " + std::to_string(m) + " while l is " + std::to_string(l));

  switch (l) {
    case 0:
      return spherical_harmonic_00;
    case 1:
      switch (m) {
        case -1:
          return spherical_harmonic_1m1;
        case 0:
          return spherical_harmonic_10;
        case 1:
          return spherical_harmonic_11;
      }
    case 2:
      switch (m) {
        case -2:
          return spherical_harmonic_2m2;
        case -1:
          return spherical_harmonic_2m1;
        case 0:
          return spherical_harmonic_20;
        case 1:
          return spherical_harmonic_21;
        case 2:
          return spherical_harmonic_22;
      }
    case 3:
      switch (m) {
        case -3:
          return spherical_harmonic_3m3;
        case -2:
          return spherical_harmonic_3m2;
        case -1:
          return spherical_harmonic_3m1;
        case 0:
          return spherical_harmonic_30;
        case 1:
          return spherical_harmonic_31;
        case 2:
          return spherical_harmonic_32;
        case 3:
          return spherical_harmonic_33;
      }
    case 4:
      switch (m) {
        case -4:
          return spherical_harmonic_4m4;
        case -3:
          return spherical_harmonic_4m3;
        case -2:
          return spherical_harmonic_4m2;
        case -1:
          return spherical_harmonic_4m1;
        case 0:
          return spherical_harmonic_40;
        case 1:
          return spherical_harmonic_41;
        case 2:
          return spherical_harmonic_42;
        case 3:
          return spherical_harmonic_43;
        case 4:
          return spherical_harmonic_44;
      }
  }

  using std::placeholders::_1;
  return std::bind(spherical_harmonic_slow, l, m, _1);
}