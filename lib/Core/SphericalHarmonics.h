#ifndef LATTICEGEN_SPHERICALHARMONICS_H
#define LATTICEGEN_SPHERICALHARMONICS_H

#include <functional>
#include <Core/Point3D.h>

namespace
{
  static const double pi = 3.14159265358979323846;
}

namespace Core
{
  // https://en.wikipedia.org/wiki/Table_of_spherical_harmonics#Real_spherical_harmonics
  // l = 0
  double spherical_harmonic_00(const Vector3D &v)
  {
    return 0.5 / sqrt(pi);
  }

  // l = 1
  double spherical_harmonic_1m1(const Vector3D &v)
  {
    return sqrt(3.0 / (4.0 * pi)) * v.y / v.length();
  }

  double spherical_harmonic_10(const Vector3D &v)
  {
    return sqrt(3.0 / (4.0 * pi)) * v.z / v.length();
  }

  double spherical_harmonic_11(const Vector3D &v)
  {
    return -sqrt(3.0 / (4.0 * pi)) * v.x / v.length();
  }

  // l = 2
  double spherical_harmonic_2m2(const Vector3D &v)
  {
    double length = v.length();
    return sqrt(15.0 / pi) / 2.0 * v.x * v.y / (length * length);
  }

  double spherical_harmonic_2m1(const Vector3D &v)
  {
    double length = v.length();
    return sqrt(15.0 / pi) / 2.0 * v.y * v.z / (length * length);
  }

  double spherical_harmonic_20(const Vector3D &v)
  {
    double length = v.length();
    return sqrt(5.0 / pi) / 4.0 * (-v.x * v.x - v.y * v.y + 2.0 * v.z * v.z) / (length * length);
  }

  double spherical_harmonic_21(const Vector3D &v)
  {
    double length = v.length();
    return -sqrt(15.0 / pi) / 2.0 * v.x * v.z / (length * length);
  }

  double spherical_harmonic_22(const Vector3D &v)
  {
    double length = v.length();
    return sqrt(15.0 / pi) / 4.0 * (v.x * v.x - v.y * v.y) / (length * length);
  }

  // l = 3
  double spherical_harmonic_3m3(const Vector3D &v)
  {
    return sqrt(35.0 / (2.0 * pi)) / 4.0 * v.y * (3.0 * v.x * v.x - v.y * v.y) / pow(v.length(), 3);
  }

  double spherical_harmonic_3m2(const Vector3D &v)
  {
    return sqrt(105.0 / pi) / 2.0 * v.x * v.y * v.z / pow(v.length(), 3);
  }

  double spherical_harmonic_3m1(const Vector3D &v)
  {
    return sqrt(21.0 / (2.0 * pi)) / 4.0 * v.y * (4.0 * v.z * v.z - v.x * v.x - v.y * v.y) / pow(v.length(), 3);
  }

  double spherical_harmonic_30(const Vector3D &v)
  {
    return sqrt(7.0 / pi) / 4.0 * v.z * (2.0 * v.z * v.z - 3.0 * v.x * v.x - 3.0 * v.y * v.y) / pow(v.length(), 3);
  }

  double spherical_harmonic_31(const Vector3D &v)
  {
    return -sqrt(21.0 / (2.0 * pi)) / 4.0 * v.x * (4.0 * v.z * v.z - v.x * v.x - v.y * v.y) / pow(v.length(), 3);
  }

  double spherical_harmonic_32(const Vector3D &v)
  {
    return sqrt(105.0 / pi) / 4.0 * v.z * (v.x * v.x - v.y * v.y) / pow(v.length(), 3);
  }

  double spherical_harmonic_33(const Vector3D &v)
  {
    return -sqrt(35.0 / (2.0 * pi)) / 4.0 * v.x * (v.x * v.x - 3.0 * v.y * v.y) / pow(v.length(), 3);
  }

  // l = 4
  double spherical_harmonic_4m4(const Vector3D &v)
  {
    return sqrt(35.0 / (2.0 * pi)) * 0.75 * v.x * v.y * (v.x * v.x - v.y * v.y) / pow(v.length(), 4);
  }

  double spherical_harmonic_4m3(const Vector3D &v)
  {
    return sqrt(35.0 / (2.0 * pi)) * 0.75 * v.y * v.z * (3.0 * v.x * v.x - v.y * v.y) / pow(v.length(), 4);
  }

  double spherical_harmonic_4m2(const Vector3D &v)
  {
    double length2 = pow(v.length(), 2);
    return sqrt(5.0 / pi) * 0.75 * v.x * v.y * (7.0 * v.z * v.z - length2) / pow(length2, 2);
  }

  double spherical_harmonic_4m1(const Vector3D &v)
  {
    double length2 = pow(v.length(), 2);
    return sqrt(5.0 / (2.0 * pi)) * 0.75 * v.y * v.z * (7.0 * v.z * v.z - 3.0 * length2) / pow(length2, 2);
  }

  double spherical_harmonic_40(const Vector3D &v)
  {
    double z2 = v.z * v.z;
    double length2 = pow(v.length(), 2);
    return sqrt(1.0 / pi) * 3.0 / 16.0 * (35.0 * z2 * z2 - 30.0 * z2 * length2 + 3.0 * length2 * length2) /
           pow(length2, 2);
  }

  double spherical_harmonic_41(const Vector3D &v)
  {
    double length2 = pow(v.length(), 2);
    return -sqrt(5.0 / (2.0 * pi)) * 0.75 * v.x * v.z * (7.0 * v.z * v.z - 3.0 * length2) / pow(length2, 2);
  }

  double spherical_harmonic_42(const Vector3D &v)
  {
    double length2 = pow(v.length(), 2);
    return sqrt(5.0 / pi) * 3.0 / 8.0 * (v.x * v.x - v.y * v.y) * (7.0 * v.z * v.z - 1.0 * length2) / pow(length2, 2);
  }

  double spherical_harmonic_43(const Vector3D &v)
  {
    return -sqrt(35.0 / (2.0 * pi)) * 0.75 * v.x * v.z * (v.x * v.x - 3.0 * v.y * v.y) / pow(v.length(), 4);
  }

  double spherical_harmonic_44(const Vector3D &v)
  {
    double x2 = v.x * v.x;
    double y2 = v.y * v.y;
    return sqrt(35.0 / pi) * 3.0 / 16.0 * (x2 * (x2 - 3.0 * y2) - y2 * (3.0 * x2 - y2)) / pow(v.length(), 4);
  }

  double spherical_harmonic_slow(unsigned int l, int m, const Vector3D &v);

  std::function<double( const Vector3D
  &)>

  spherical_harmonic(unsigned int l, int m);
}

#endif //LATTICEGEN_SPHERICALHARMONICS_H
