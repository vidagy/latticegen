#ifndef LATTICEGEN_FACTORIAL_H
#define LATTICEGEN_FACTORIAL_H

#include <Core/Exceptions.h>
#include <complex>

namespace
{
  static const double pi = M_PI;
}

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

  /// @brief this function calculates \Gamma(n+1/2)
  inline double gamma_plus_half(unsigned int n)
  {
    static const double gamma_cache[16] =
      {1.0 * sqrt(pi), 0.5 * sqrt(pi), 0.75 * sqrt(pi), 1.875 * sqrt(pi), 6.5625 * sqrt(pi), 29.53125 * sqrt(pi),
       162.421875 * sqrt(pi), 1055.7421875 * sqrt(pi), 7918.06640625 * sqrt(pi), 67303.564453125 * sqrt(pi),
       639383.8623046875 * sqrt(pi), 6.71353055419921875e6 * sqrt(pi), 7.7205601373291015625e7 * sqrt(pi),
       9.650700171661376953125e8 * sqrt(pi), 1.302844523174285888671875e10 * sqrt(pi),
       1.889124558602714538574219e11 * sqrt(pi)
      };
    if (n < 16)
      return gamma_cache[n];
    else {
      return double_factorial(2 * n - 1) / static_cast<double>(1ull << n) * sqrt(pi);
    }
  }

  template<typename T>
  T pow(T base, int exponent)
  {
    if (exponent < 0)
      THROW_INVALID_ARGUMENT("Negative exponent = " + std::to_string(exponent));

    T result = 1.0;

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
        return 1.0i;
      case 2:
        return -1.0;
      case 3:
        return -1.0i;
      default:
        THROW_LOGIC_ERROR("modulo 4 returned" + std::to_string(mod));
    }
  }
}
#endif //LATTICEGEN_FACTORIAL_H
