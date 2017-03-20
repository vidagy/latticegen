#include <Math/Integrator.h>

using namespace Math;

double IntegratorEquidistant::simpson(const std::vector<double> &f, double dx)
{
  auto f_size = f.size();

  if (f_size < 2)
    THROW_INVALID_ARGUMENT("in IntegratorEquidistant::simpson, vector size is < 2");
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

double IntegratorEquidistant::simpson_alt(const std::vector<double> &f, double dx)
{
  auto f_size = f.size();

  if (f_size < 8)
    return simpson(f, dx);
  else {
    const double dx_17_p_48 = 17.0 * dx / 48.0;
    const double dx_59_p_48 = 59.0 * dx / 48.0;
    const double dx_43_p_48 = 43.0 * dx / 48.0;
    const double dx_49_p_48 = 49.0 * dx / 48.0;

    double pre = 0.0;
    pre += dx_17_p_48 * f[0];
    pre += dx_59_p_48 * f[1];
    pre += dx_43_p_48 * f[2];
    pre += dx_49_p_48 * f[3];

    auto i = 4u;
    double result = 0.0;
    while (i < f_size - 4) {
      result += dx * f[i++];
    }
    double post = 0.0;
    post += dx_49_p_48 * f[i];
    post += dx_43_p_48 * f[i + 1];
    post += dx_59_p_48 * f[i + 2];
    post += dx_17_p_48 * f[i + 3];

    return pre + result + post;
  }
}

namespace
{
  static const double consts_1[1] = {0.5};
  static const double consts_3[3] = {9.0 / 24.0, 28.0 / 24.0, 23.0 / 24.0};
  static const double consts_5[5] = {475.0 / 1440.0, 1902.0 / 1440.0, 1104.0 / 1440.0, 1586.0 / 1440.0,
                                     1413.0 / 1440.0};
  static const double consts_7[7] = {36799.0 / 120960.0, 176648.0 / 120960.0, 54851.0 / 120960.0, 177984.0 / 120960.0,
                                     89437.0 / 120960.0, 130936.0 / 120960.0, 119585.0 / 120960.0};

  static const double *trapezoid_consts[4] = {consts_1, consts_3, consts_5, consts_7};
}

double IntegratorEquidistant::trapezoidal(const std::vector<double> &f, double dx, unsigned int quadrature)
{
  if (!(quadrature == 1 || quadrature == 3 || quadrature == 5 || quadrature == 7))
    THROW_INVALID_ARGUMENT("invalid quadrature = " + std::to_string(quadrature) + " in trapezoidal rule ");
  if (f.size() < 2)
    THROW_INVALID_ARGUMENT("f.size() = " + std::to_string(f.size()) + " in trapezoidal rule, must be at least 2");
  if (f.size() < 2 * quadrature) {
    auto new_quadrature = static_cast<unsigned int>(f.size()) / 2;
    new_quadrature = (new_quadrature & 1) == 0 ? new_quadrature - 1 : new_quadrature;
    return trapezoidal(f, dx, new_quadrature);
  }

  const auto &consts = trapezoid_consts[quadrature / 2];

  auto i = 0u;
  double res = 0.0;
  for (; i < quadrature; ++i) {
    res += consts[i] * f[i];
  }
  for (; i < f.size() - quadrature; ++i) {
    res += f[i];
  }
  for (; i < f.size(); ++i) {
    res += consts[f.size() - i - 1] * f[i];
  }
  return res * dx;
}

double IntegratorExponential::simpson(const std::vector<double> &f, const Core::ExponentialMesh &x)
{
  auto f_size = f.size();
  const auto &d_points = x.d_points;
  const auto &dx = x.dx;

  if (f_size != d_points.size())
    THROW_INVALID_ARGUMENT(
      "in IntegratorExponential::simpson, vector sizes are unequal f=" + std::to_string(f_size) +
      " x=" + std::to_string(d_points.size()));
  if (f_size < 2)
    THROW_INVALID_ARGUMENT("in IntegratorExponential::simpson, vector size is < 2");
  else if (f_size == 2)
    return dx * (f[0] * d_points[0] + f[1] * d_points[1]) / 2.0;
  else {
    double result = 0.0;

    const double dx_p_3 = dx / 3.0;
    const double dx_4_p_3 = 4.0 * dx / 3.0;
    const double dx_2_p_3 = 2.0 * dx / 3.0;

    result += dx_p_3 * f[0] * d_points[0];
    auto i = 1u;
    while (i < f_size - 3) {
      result += dx_4_p_3 * f[i] * d_points[i];
      result += dx_2_p_3 * f[i + 1] * d_points[i + 1];
      i += 2;
    }
    result += dx_4_p_3 * f[i] * d_points[i];
    ++i;
    result += dx_p_3 * f[i] * d_points[i];
    ++i;

    if (i == f_size - 1) // size is even
    {
      // warning! if we get here the error may potentially increase
      result += dx * (f[i - 1] * d_points[i - 1] + f[i] * d_points[i]) / 2.0;
    }

    return result;
  }
}

double IntegratorExponential::simpson_alt(const std::vector<double> &f, const Core::ExponentialMesh &x)
{
  auto f_size = f.size();
  const auto &d_points = x.d_points;
  const auto &dx = x.dx;

  if (f_size != d_points.size())
    THROW_INVALID_ARGUMENT(
      "in IntegratorExponential::simpson, vector sizes are unequal f=" + std::to_string(f_size) +
      " x=" + std::to_string(d_points.size()));

  if (f_size < 8)
    return IntegratorExponential::simpson(f, x);
  else {
    const double dx_17_p_48 = 17.0 * dx / 48.0;
    const double dx_59_p_48 = 59.0 * dx / 48.0;
    const double dx_43_p_48 = 43.0 * dx / 48.0;
    const double dx_49_p_48 = 49.0 * dx / 48.0;

    double pre = 0.0;
    pre += dx_17_p_48 * f[0] * d_points[0];
    pre += dx_59_p_48 * f[1] * d_points[1];
    pre += dx_43_p_48 * f[2] * d_points[2];
    pre += dx_49_p_48 * f[3] * d_points[3];

    auto i = 4u;
    double result = 0.0;
    while (i < f_size - 4) {
      result += dx * f[i] * d_points[i];
      ++i;
    }
    double post = 0.0;
    post += dx_49_p_48 * f[i] * d_points[i];
    post += dx_43_p_48 * f[i + 1] * d_points[i + 1];
    post += dx_59_p_48 * f[i + 2] * d_points[i + 2];
    post += dx_17_p_48 * f[i + 3] * d_points[i + 3];

    return pre + result + post;
  }
}

double IntegratorGeneric::integrate(const std::vector<double> &f, const std::vector<double> &x)
{
  auto f_size = f.size();

  if (f_size != x.size())
    THROW_INVALID_ARGUMENT(
      "in IntegratorGeneric::simpson, vector sizes are unequal f=" + std::to_string(f_size) +
      " x=" + std::to_string(x.size()));
  else if (f_size < 2)
    THROW_INVALID_ARGUMENT("in IntegratorGeneric::simpson, vector size is < 2");
  else {
    double result = 0.0;

    result += f[0] * (x[1] - x[0]);
    auto i = 1u;
    while (i < f_size - 1) {
      result += f[i] * (x[i + 1] - x[i - 1]);
      i += 1;
    }
    result += f[i] * (x[i] - x[i - 1]);

    return result / 2.0;
  }
}
