#ifndef LATTICEGEN_LEGENDRE_H
#define LATTICEGEN_LEGENDRE_H

#include <functional>
#include "Factorial.h"

namespace Core
{
  double legendre_polynomial_slow(int l, int m, double x)
  {
    // l >= 0, -l <= m <= l, , -1 <= x <= 1
    double pmm = 1.0;
    double sign = ((m & 1) == 0 ? 1 : -1);
    if (m < 0) {
      m = -m;
      pmm = sign * factorial(l - m) / factorial(l + m);
    }

    if (m > 0) {
      pmm *= sign * double_factorial(2 * m - 1) * pow(1 - x * x, m / 2.0);
    }
    if (l == m) {
      return pmm;
    }

    // P^m_{m+1}(x) = x*(2m + 1)*P^m_m(x)
    double next = x * (2 * m + 1) * pmm;
    if (l == m + 1) {
      return next;
    }

    // P^m_l(x) = ( x*(2l - 1)*P^m_{l-1} - (l + m - 1)*P^m_{l-2} ) / (l - m)
    for (int n = m + 2; n <= l; n++) {
      double nnext = (x * (2 * n - 1) * next - (n + m - 1) * pmm) / (n - m);
      pmm = next;
      next = nnext;
    }
    return next;
  }

  double zero(double x)
  {
    return 0.0;
  }

  double legendre_polynomial_00(double x)
  {
    return 1.0;
  };

  double legendre_polynomial_1m1(double x)
  {
    return sqrt(1 - x * x) / 2.0;
  };

  double legendre_polynomial_10(double x)
  {
    return x;
  };

  double legendre_polynomial_11(double x)
  {
    return -1 * sqrt(1 - x * x);
  };

  double legendre_polynomial_2m2(double x)
  {
    return (1 - x * x) / 8.0;
  };

  double legendre_polynomial_2m1(double x)
  {
    return x * sqrt(1 - x * x) / 2.0;
  };

  double legendre_polynomial_20(double x)
  {
    return (3.0 * x * x - 1.0) / 2.0;
  };

  double legendre_polynomial_21(double x)
  {
    return -3.0 * x * sqrt(1 - x * x);
  };

  double legendre_polynomial_22(double x)
  {
    return 3.0 * (1 - x * x);
  };

  double legendre_polynomial_3m3(double x)
  {
    return pow(1 - x * x, 1.5) / 48.0;
  };

  double legendre_polynomial_3m2(double x)
  {
    return x * (1 - x * x) / 8.0;
  };

  double legendre_polynomial_3m1(double x)
  {
    double xx = x * x;
    return (5.0 * xx - 1.0) * sqrt(1 - xx) / 8.0;
  };

  double legendre_polynomial_30(double x)
  {
    return (5.0 * x * x - 3.0) * x / 2.0;
  };

  double legendre_polynomial_31(double x)
  {
    double xx = x * x;
    return -1.5 * (5.0 * xx - 1.0) * sqrt(1 - xx);
  };

  double legendre_polynomial_32(double x)
  {
    return 15.0 * x * (1 - x * x);
  };

  double legendre_polynomial_33(double x)
  {
    return -15.0 * pow(1 - x * x, 1.5);
  };

  double legendre_polynomial_4m4(double x)
  {
    return pow(1 - x * x, 2) / 384.0;
  };

  double legendre_polynomial_4m3(double x)
  {
    return x * pow(1 - x * x, 1.5) / 48.0;
  };

  double legendre_polynomial_4m2(double x)
  {
    double xx = x * x;
    return (7.0 * xx - 1.0) * (1 - xx) / 48.0;
  };

  double legendre_polynomial_4m1(double x)
  {
    double xx = x * x;
    return (7.0 * xx * x - 3.0 * x) * sqrt(1 - xx) / 8.0;
  };

  double legendre_polynomial_40(double x)
  {
    double xx = x * x;
    return (35.0 * xx * xx - 30.0 * xx + 3.0) / 8.0;
  };

  double legendre_polynomial_41(double x)
  {
    double xx = x * x;
    return -2.5 * (7.0 * xx * x - 3.0 * x) * sqrt(1 - xx);
  };

  double legendre_polynomial_42(double x)
  {
    double xx = x * x;
    return 7.5 * (7.0 * xx - 1.0) * (1 - xx);
  };

  double legendre_polynomial_43(double x)
  {
    return -105.0 * x * pow(1 - x * x, 1.5);
  };

  double legendre_polynomial_44(double x)
  {
    return 105.0 * pow(1 - x * x, 2);
  };

  std::function<double(double)> legendre_polynomial(int l, int m)
  {
    if (l < 0)
      return legendre_polynomial(-l - 1, m);

    if (abs(m) > l)
      return zero;

    switch (l) {
      case 0:
        return legendre_polynomial_00;
      case 1:
        switch (m) {
          case -1:
            return legendre_polynomial_1m1;
          case 0:
            return legendre_polynomial_10;
          case 1:
            return legendre_polynomial_11;
        }
      case 2:
        switch (m) {
          case -2:
            return legendre_polynomial_2m2;
          case -1:
            return legendre_polynomial_2m1;
          case 0:
            return legendre_polynomial_20;
          case 1:
            return legendre_polynomial_21;
          case 2:
            return legendre_polynomial_22;
        }
      case 3:
        switch (m) {
          case -3:
            return legendre_polynomial_3m3;
          case -2:
            return legendre_polynomial_3m2;
          case -1:
            return legendre_polynomial_3m1;
          case 0:
            return legendre_polynomial_30;
          case 1:
            return legendre_polynomial_31;
          case 2:
            return legendre_polynomial_32;
          case 3:
            return legendre_polynomial_33;
        }
      case 4:
        switch (m) {
          case -4:
            return legendre_polynomial_4m4;
          case -3:
            return legendre_polynomial_4m3;
          case -2:
            return legendre_polynomial_4m2;
          case -1:
            return legendre_polynomial_4m1;
          case 0:
            return legendre_polynomial_40;
          case 1:
            return legendre_polynomial_41;
          case 2:
            return legendre_polynomial_42;
          case 3:
            return legendre_polynomial_43;
          case 4:
            return legendre_polynomial_44;
        }
    }

    using std::placeholders::_1;
    return std::bind(legendre_polynomial_slow, l, m, _1);
  }
}

#endif //LATTICEGEN_LEGENDRE_H
