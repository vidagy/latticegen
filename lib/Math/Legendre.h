#ifndef LATTICEGEN_LEGENDRE_H
#define LATTICEGEN_LEGENDRE_H

#include <cmath>
#include <functional>
#include <Math/Factorial.h>

namespace Math
{
  inline double zero(double)
  {
    return 0.0;
  }

  // l = 0
  inline double legendre_polynomial_00(double)
  {
    return 1.0;
  };

  // l = 1
  inline double legendre_polynomial_1m1(double x)
  {
    return sqrt(1 - x * x) / 2.0;
  };

  inline double legendre_polynomial_10(double x)
  {
    return x;
  };

  inline double legendre_polynomial_11(double x)
  {
    return -1 * sqrt(1 - x * x);
  };

  // l = 2
  inline double legendre_polynomial_2m2(double x)
  {
    return (1 - x * x) / 8.0;
  };

  inline double legendre_polynomial_2m1(double x)
  {
    return x * sqrt(1 - x * x) / 2.0;
  };

  inline double legendre_polynomial_20(double x)
  {
    return (3.0 * x * x - 1.0) / 2.0;
  };

  inline double legendre_polynomial_21(double x)
  {
    return -3.0 * x * sqrt(1 - x * x);
  };

  inline double legendre_polynomial_22(double x)
  {
    return 3.0 * (1 - x * x);
  };

  // l = 3
  inline double legendre_polynomial_3m3(double x)
  {
    return (1.0 - x * x) * sqrt(1.0 - x * x) / 48.0;
  };

  inline double legendre_polynomial_3m2(double x)
  {
    return x * (1 - x * x) / 8.0;
  };

  inline double legendre_polynomial_3m1(double x)
  {
    double xx = x * x;
    return (5.0 * xx - 1.0) * sqrt(1 - xx) / 8.0;
  };

  inline double legendre_polynomial_30(double x)
  {
    return (5.0 * x * x - 3.0) * x / 2.0;
  };

  inline double legendre_polynomial_31(double x)
  {
    double xx = x * x;
    return -1.5 * (5.0 * xx - 1.0) * sqrt(1 - xx);
  };

  inline double legendre_polynomial_32(double x)
  {
    return 15.0 * x * (1 - x * x);
  };

  inline double legendre_polynomial_33(double x)
  {
    return -15.0 * (1.0 - x * x) * sqrt(1.0 - x * x);
  };

  // l = 4
  inline double legendre_polynomial_4m4(double x)
  {
    return (1.0 - x * x) * (1.0 - x * x) / 384.0;
  };

  inline double legendre_polynomial_4m3(double x)
  {
    return x * (1.0 - x * x) * sqrt(1.0 - x * x) / 48.0;
  };

  inline double legendre_polynomial_4m2(double x)
  {
    double xx = x * x;
    return (7.0 * xx - 1.0) * (1 - xx) / 48.0;
  };

  inline double legendre_polynomial_4m1(double x)
  {
    double xx = x * x;
    return (7.0 * xx * x - 3.0 * x) * sqrt(1 - xx) / 8.0;
  };

  inline double legendre_polynomial_40(double x)
  {
    double xx = x * x;
    return (35.0 * xx * xx - 30.0 * xx + 3.0) / 8.0;
  };

  inline double legendre_polynomial_41(double x)
  {
    double xx = x * x;
    return -2.5 * (7.0 * xx * x - 3.0 * x) * sqrt(1 - xx);
  };

  inline double legendre_polynomial_42(double x)
  {
    double xx = x * x;
    return 7.5 * (7.0 * xx - 1.0) * (1 - xx);
  };

  inline double legendre_polynomial_43(double x)
  {
    return -105.0 * x * (1.0 - x * x) * sqrt(1.0 - x * x);
  };

  inline double legendre_polynomial_44(double x)
  {
    return 105.0 * (1.0 - x * x) * (1.0 - x * x);
  };

  double legendre_polynomial_slow(int l, int m, double x);

  std::function<double(double)> legendre_polynomial(int l, int m);
}

#endif //LATTICEGEN_LEGENDRE_H
