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
      // warning! if we get here the error may potentially increase
      result += dx * (f[i - 1] + f[i]) / 2.0;
    }

    return result;
  }
}

double IntegratorEquidistant::simpson_3_8(const std::vector<double> &f, double dx)
{
  int f_size = static_cast<int>(f.size());

  if (f_size < 4)
    return simpson(f, dx);
  else {
    double result = 0.0;

    const double dx_3_p_8 = 3.0 * dx / 8.0;
    const double dx_9_p_8 = 9.0 * dx / 8.0;
    const double dx_6_p_8 = 6.0 * dx / 8.0;

    result += dx_3_p_8 * f[0];
    int i = 1;
    while (i < f_size - 5) {
      result += dx_9_p_8 * f[i];
      result += dx_9_p_8 * f[i + 1];
      result += dx_6_p_8 * f[i + 2];
      i += 3;
    }
    result += dx_9_p_8 * f[i++];
    result += dx_9_p_8 * f[i++];
    result += dx_3_p_8 * f[i++];

    if (i == f_size - 2) {
      result += dx * (f[i - 1] + 4.0 * f[i] + f[i + 1]) / 3.0;
    } else if (i == f_size - 1) {
      // warning! if we get here the error may potentially increase
      result += dx * (f[i - 1] + f[i]) / 2.0;
    }

    return result;
  }
}
