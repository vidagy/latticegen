#ifndef LATTICEGEN_FACTORIAL_H
#define LATTICEGEN_FACTORIAL_H

#include <Core/Exceptions.h>
#include <complex>

namespace Math
{
  inline double factorial(int n)
  {
    static const double factorial_cache[16] =
      {1, 1, 2, 6, 24, 120, 720, 5040,
       40320, 362880, 3628800, 39916800, 479001600, 6227020800, 87178291200, 1307674368000};

    if (n < 0)
      THROW_INVALID_ARGUMENT("factorial invoked with n = " + std::to_string(n));

    if (n < 16)
      return factorial_cache[n];
    else {
      double result = 1.0;
      while (n > 15) {
        result *= (double) n;
        --n;
      }
      return result * factorial_cache[n];
    }
  }

  inline double double_factorial(int n)
  {
    static const double double_factorial_cache[16] =
      {1, 1, 2, 3, 8, 15, 48, 105, 384, 945, 3840, 10395, 46080, 135135, 645120, 2027025};

    if (n < 0)
      THROW_INVALID_ARGUMENT("double factorial invoked with n = " + std::to_string(n));

    if (n < 16)
      return double_factorial_cache[n];
    else {
      double result = 1.0;
      while (n > 15) {
        result *= (double) n;
        n -= 2;
      }
      return result * double_factorial_cache[n];
    }
  }

  inline double pow(double base, int exponent)
  {
    if (exponent < 0)
      THROW_INVALID_ARGUMENT("Negative exponent = " + std::to_string(exponent));

    double result = 1.0;

    while (exponent) {
      if (exponent & 1)
        result *= base;
      exponent >>= 1;
      base *= base;
    }

    return result;
  }

  inline unsigned int modulo(int num, unsigned int base)
  {
    if (num < 0)
      return (num % base) + base;
    else
      return num % base;
  }

  inline std::complex<double> ipow(int exponent)
  {
    using namespace std::complex_literals;

    auto mod = modulo(exponent, 4);
    switch (mod) {
      case 0:
        return 1.0;
      case 1:
        return 1i;
      case 2:
        return -1.0;
      case 3:
        return -1i;
      default:
        THROW_LOGIC_ERROR("modulo 4 returned" + std::to_string(mod));
    }
  }
}
#endif //LATTICEGEN_FACTORIAL_H
