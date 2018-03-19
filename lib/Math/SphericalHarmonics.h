#ifndef LATTICEGEN_SPHERICALHARMONICS_H
#define LATTICEGEN_SPHERICALHARMONICS_H

#include <Core/Point3D.h>
#include <Math/CommonFunctions.h>

#include <functional>
#include <complex>
#include <boost/config/no_tr1/complex.hpp>

using namespace Core;

namespace Math
{
  // https://en.wikipedia.org/wiki/Table_of_spherical_harmonics#Real_spherical_harmonics
  namespace Real
  {
    // l = 0
    LATTICEGEN_MUTE_BEGIN
    LATTICEGEN_MUTE_UNUSED_VAR

    inline double spherical_harmonic_00(const Point3D &v)
    {
      LATTICEGEN_MUTE_END
      return 0.5 / sqrt(pi);
    }

    // l = 1
    inline double spherical_harmonic_1m1(const Point3D &v)
    {
      return -sqrt(3.0 / (4.0 * pi)) * v(1) / v.norm();
    }

    inline double spherical_harmonic_10(const Point3D &v)
    {
      return sqrt(3.0 / (4.0 * pi)) * v(2) / v.norm();
    }

    inline double spherical_harmonic_11(const Point3D &v)
    {
      return -sqrt(3.0 / (4.0 * pi)) * v(0) / v.norm();
    }

    // l = 2
    inline double spherical_harmonic_2m2(const Point3D &v)
    {
      double length = v.norm();
      return sqrt(15.0 / pi) / 2.0 * v(0) * v(1) / (length * length);
    }

    inline double spherical_harmonic_2m1(const Point3D &v)
    {
      double length = v.norm();
      return -sqrt(15.0 / pi) / 2.0 * v(1) * v(2) / (length * length);
    }

    inline double spherical_harmonic_20(const Point3D &v)
    {
      double length = v.norm();
      return sqrt(5.0 / pi) / 4.0 * (-v(0) * v(0) - v(1) * v(1) + 2.0 * v(2) * v(2)) / (length * length);
    }

    inline double spherical_harmonic_21(const Point3D &v)
    {
      double length = v.norm();
      return -sqrt(15.0 / pi) / 2.0 * v(0) * v(2) / (length * length);
    }

    inline double spherical_harmonic_22(const Point3D &v)
    {
      double length = v.norm();
      return sqrt(15.0 / pi) / 4.0 * (v(0) * v(0) - v(1) * v(1)) / (length * length);
    }

    // l = 3
    inline double spherical_harmonic_3m3(const Point3D &v)
    {
      double length = v.norm();
      return -sqrt(35.0 / (2.0 * pi)) / 4.0 * v(1) * (3.0 * v(0) * v(0) - v(1) * v(1)) / (length * length * length);
    }

    inline double spherical_harmonic_3m2(const Point3D &v)
    {
      double length = v.norm();
      return sqrt(105.0 / pi) / 2.0 * v(0) * v(1) * v(2) / (length * length * length);
    }

    inline double spherical_harmonic_3m1(const Point3D &v)
    {
      double length = v.norm();
      return -sqrt(21.0 / (2.0 * pi)) / 4.0 * v(1) * (4.0 * v(2) * v(2) - v(0) * v(0) - v(1) * v(1)) /
             (length * length * length);
    }

    inline double spherical_harmonic_30(const Point3D &v)
    {
      double length = v.norm();
      return sqrt(7.0 / pi) / 4.0 * v(2) * (2.0 * v(2) * v(2) - 3.0 * v(0) * v(0) - 3.0 * v(1) * v(1)) /
             (length * length * length);
    }

    inline double spherical_harmonic_31(const Point3D &v)
    {
      double length = v.norm();
      return -sqrt(21.0 / (2.0 * pi)) / 4.0 * v(0) * (4.0 * v(2) * v(2) - v(0) * v(0) - v(1) * v(1)) /
             (length * length * length);
    }

    inline double spherical_harmonic_32(const Point3D &v)
    {
      double length = v.norm();
      return sqrt(105.0 / pi) / 4.0 * v(2) * (v(0) * v(0) - v(1) * v(1)) / (length * length * length);
    }

    inline double spherical_harmonic_33(const Point3D &v)
    {
      double length = v.norm();
      return -sqrt(35.0 / (2.0 * pi)) / 4.0 * v(0) * (v(0) * v(0) - 3.0 * v(1) * v(1)) / (length * length * length);
    }

    // l = 4
    inline double spherical_harmonic_4m4(const Point3D &v)
    {
      double length = v.norm();
      return sqrt(35.0 / pi) * 0.75 * v(0) * v(1) * (v(0) * v(0) - v(1) * v(1)) / (length * length * length * length);
    }

    inline double spherical_harmonic_4m3(const Point3D &v)
    {
      double length = v.norm();
      return -sqrt(35.0 / (2.0 * pi)) * 0.75 * v(1) * v(2) * (3.0 * v(0) * v(0) - v(1) * v(1)) /
             (length * length * length * length);
    }

    inline double spherical_harmonic_4m2(const Point3D &v)
    {
      double length2 = v.squaredNorm();
      return sqrt(5.0 / pi) * 0.75 * v(0) * v(1) * (7.0 * v(2) * v(2) - length2) / (length2 * length2);
    }

    inline double spherical_harmonic_4m1(const Point3D &v)
    {
      double length2 = v.squaredNorm();
      return -sqrt(5.0 / (2.0 * pi)) * 0.75 * v(1) * v(2) * (7.0 * v(2) * v(2) - 3.0 * length2) / (length2 * length2);
    }

    inline double spherical_harmonic_40(const Point3D &v)
    {
      double z2 = v(2) * v(2);
      double length2 = v.squaredNorm();
      return sqrt(1.0 / pi) * 3.0 / 16.0 * (35.0 * z2 * z2 - 30.0 * z2 * length2 + 3.0 * length2 * length2) /
             (length2 * length2);
    }

    inline double spherical_harmonic_41(const Point3D &v)
    {
      double length2 = v.squaredNorm();
      return -sqrt(5.0 / (2.0 * pi)) * 0.75 * v(0) * v(2) * (7.0 * v(2) * v(2) - 3.0 * length2) / (length2 * length2);
    }

    inline double spherical_harmonic_42(const Point3D &v)
    {
      double length2 = v.squaredNorm();
      return sqrt(5.0 / pi) * 3.0 / 8.0 * (v(0) * v(0) - v(1) * v(1)) * (7.0 * v(2) * v(2) - 1.0 * length2) /
             (length2 * length2);
    }

    inline double spherical_harmonic_43(const Point3D &v)
    {
      double length = v.norm();
      return -sqrt(35.0 / (2.0 * pi)) * 0.75 * v(0) * v(2) * (v(0) * v(0) - 3.0 * v(1) * v(1)) /
             (length * length * length * length);
    }

    inline double spherical_harmonic_44(const Point3D &v)
    {
      double x2 = v(0) * v(0);
      double y2 = v(1) * v(1);
      double length = v.norm();
      return sqrt(35.0 / pi) * 3.0 / 16.0 * (x2 * (x2 - 3.0 * y2) - y2 * (3.0 * x2 - y2)) /
             (length * length * length * length);
    }

    double spherical_harmonic_slow(unsigned int l, int m, const Point3D &v);

    double spherical_harmonic(unsigned int l, int m, const Point3D &v);
  }

  namespace Complex
  {
    std::complex<double> spherical_harmonic_slow(unsigned int l, int m, const Point3D &v);

    std::complex<double> spherical_harmonic(unsigned int l, int m, const Point3D &v);
  }
}

#endif //LATTICEGEN_SPHERICALHARMONICS_H
