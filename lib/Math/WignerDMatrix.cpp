#include <iostream>
#include "WignerDMatrix.h"
#include "CommonFunctions.h"

using namespace Math;

namespace
{

}

Eigen::MatrixXcd WignerDMatrix::calculate(double alpha, double beta, double gamma, unsigned int l)
{
  using namespace std::complex_literals;

  auto l_int = (int) l;
  Eigen::VectorXcd alpha_rot(2 * l_int + 1);
  Eigen::VectorXcd gamma_rot(2 * l_int + 1);
  for (auto m = -l_int; m <= l_int; ++m) {
    alpha_rot(l + m) = std::exp(-1.0i * alpha * ((double) m));
    gamma_rot(l + m) = std::exp(-1.0i * gamma * ((double) m));
  }
  Eigen::MatrixXd small_d = calculate_small(beta, l);
  Eigen::MatrixXcd res = alpha_rot.asDiagonal() * small_d.cast<std::complex<double>>() * gamma_rot.asDiagonal();
  return res;
}

Eigen::MatrixXd WignerDMatrix::calculate_small(double beta, unsigned int l)
{
  const int l_int = (int) l;

  Eigen::MatrixXd res(2 * l_int + 1, 2 * l_int + 1);
  for (auto m_p = -l_int; m_p <= l_int; ++m_p) {
    for (auto m = -l_int; m <= l_int; ++m) {
      auto part = 0.0;

      auto s_min = std::max(0, m - m_p);
      auto s_max = std::min(l_int + m, l_int - m_p);

      for (auto s = s_min; s <= s_max; ++s) {
        part +=
          sign(m_p - m + s)
          / (factorial(l_int + m - s) * factorial(s) * factorial(m_p - m + s) * factorial(l_int - m_p - s))
          * Math::pow(cos(beta / 2.0), 2 * l_int + m - m_p - 2 * s)
          * Math::pow(sin(beta / 2.0), m_p - m + 2 * s);
      }

      res(l_int + m_p, l_int + m) =
        part * sqrt(factorial(l_int + m_p) * factorial(l_int - m_p) * factorial(l_int + m) * factorial(l_int - m));
    }
  }
  return res;
}
