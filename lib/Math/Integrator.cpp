#include <Math/Integrator.h>

#include <stdexcept>

using namespace Math;

double IntegratorEquidistant::simpson(const std::vector<double> &f, double dx)
{
  auto f_size = f.size();

  if (f_size < 2)
    throw std::invalid_argument("in Integrator::simpson, vector size is < 2");
  else if (f_size == 2)
    return dx * (f[0] + f[1]) / 2.0;
  else {
    double result = 0.0;

    const double dx_p_3 = dx / 3.0;
    const double dx_4_p_3 = 4.0 * dx / 3.0;
    const double dx_2_p_3 = 2.0 * dx / 3.0;

    result += dx_p_3 * f[0];
    auto i = 1u;
    while (i < f_size - 3) {
      result += dx_4_p_3 * f[i];
      result += dx_2_p_3 * f[i + 1];
      i += 2;
    }
    result += dx_4_p_3 * f[i++];
    result += dx_p_3 * f[i++];

    if (i == f_size - 1) // size is even
    {
      result += dx * (f[i - 1] + f[i]) / 2.0;
    }

    return result;
  }
}
