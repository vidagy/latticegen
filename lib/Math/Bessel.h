#ifndef LATTICEGEN_BESSEL_H
#define LATTICEGEN_BESSEL_H

#include <cmath>
#include <utility>
#include <complex>

namespace Math
{

  // these are the spherical neumann, bessel and hankel+ functions as defined on wikipedia
  template<typename T>
  T bessel_0(const T &x) { return sin(x) / x; }

  template<typename T>
  T bessel_1(const T &x) { return (sin(x) / x - cos(x)) / x; }

  template<typename T>
  T bessel_2(const T &x) { return ((3.0 / x / x - 1.0) * sin(x) - 3.0 * cos(x) / x) / x; }

  template<typename T>
  T bessel_3(const T &x) { return ((15.0 / x / x - 6.0) * sin(x) / x - (15.0 / x / x - 1.0) * cos(x)) / x; }
  // write some more if necessary

  namespace
  {
    template<typename T>
    T slow(unsigned int n, const T &x, T(*zero)(const T &), T(*first)(const T &))
    {
      if (n == 0)
        return zero(x);
      if (n == 1)
        return first(x);

      auto im2 = zero(x);
      auto im1 = first(x);
      for (auto i = 2u; i <= n; ++i) {
        auto res = (2.0 * i - 1.0) / x * im1 - im2;
        std::swap(im1, im2);
        std::swap(res, im1);
      }
      return im1;
    }
  }

  template<typename T>
  T bessel_slow(unsigned int n, const T &x)
  {
    return slow(n, x, bessel_0<T>, bessel_1<T>);
  }

  template<typename T>
  T bessel(unsigned int n, const T &x)
  {
    switch (n) {
      case 0:
        return bessel_0<T>(x);
      case 1:
        return bessel_1<T>(x);
      case 2:
        return bessel_2<T>(x);
      case 3:
        return bessel_3<T>(x);
      default:
        return bessel_slow<T>(n, x);
    }
  }


  template<typename T>
  T neumann_0(const T &x) { return -cos(x) / x; }

  template<typename T>
  T neumann_1(const T &x) { return (-cos(x) / x - sin(x)) / x; }

  template<typename T>
  T neumann_2(const T &x) { return ((-3.0 / x / x + 1.0) * cos(x) - 3.0 * sin(x) / x) / x; }

  template<typename T>
  T neumann_3(const T &x) { return ((-15.0 / x / x + 6.0) * cos(x) / x - (15.0 / x / x - 1.0) * sin(x)) / x; }
  // write some more if necessary

  template<typename T>
  T neumann_slow(unsigned int n, const T &x)
  {
    return slow(n, x, neumann_0<T>, neumann_1<T>);
  }

  template<typename T>
  T neumann(unsigned int n, const T &x)
  {
    switch (n) {
      case 0:
        return neumann_0<T>(x);
      case 1:
        return neumann_1<T>(x);
      case 2:
        return neumann_2<T>(x);
      case 3:
        return neumann_3<T>(x);
      default:
        return neumann_slow<T>(n, x);
    }
  }

  std::complex<double> hankel_1(unsigned int n, const std::complex<double> &x)
  {
    using namespace std::complex_literals;
    return bessel<std::complex<double>>(n, x) + 1i * neumann<std::complex<double>>(n, x);
  }

  std::complex<double> hankel_2(unsigned int n, const std::complex<double> &x)
  {
    using namespace std::complex_literals;
    return bessel<std::complex<double>>(n, x) - 1i * neumann<std::complex<double>>(n, x);
  }
}


#endif //LATTICEGEN_BESSEL_H
