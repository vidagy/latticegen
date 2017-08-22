#ifndef LATTICEGEN_INTEGRATOR_H
#define LATTICEGEN_INTEGRATOR_H

#include <vector>
#include <Core/RadialMesh.h>

namespace Math
{
  class IntegratorEquidistant
  {
  public:
    static double simpson(const std::vector<double> &f, double dx);
    static double simpson_3_8(const std::vector<double> &f, double dx);
    static double simpson_alt(const std::vector<double> &f, double dx);

    static double trapezoidal(const std::vector<double> &f, double dx, unsigned int quadrature = 7);
  };

  class IntegratorExponential
  {
  public:
    template<typename T>
    static T simpson(const std::vector<T> &f, const Core::ExponentialMesh &x)
    {
      auto f_size = f.size();
      const auto &d_points = x.get_d_points();
      const auto &dx = x.get_dx();

      if (f_size != d_points.size())
        THROW_INVALID_ARGUMENT(
          "in IntegratorExponential::simpson, vector sizes are unequal f=" + std::to_string(f_size) +
          " x=" + std::to_string(d_points.size()));
      if (f_size < 2)
        THROW_INVALID_ARGUMENT("in IntegratorExponential::simpson, vector size is < 2");
      else if (f_size == 2)
        return dx * (f[0] * d_points[0] + f[1] * d_points[1]) / 2.0;
      else {
        T result = 0.0;

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


    template<typename T>
    static T simpson_alt(const std::vector<T> &f, const Core::ExponentialMesh &x)
    {
      auto f_size = f.size();
      const auto &d_points = x.get_d_points();
      const auto &dx = x.get_dx();

      if (f_size != d_points.size())
        THROW_INVALID_ARGUMENT(
          "in IntegratorExponential::simpson, vector sizes are unequal f=" + std::to_string(f_size) +
          " x=" + std::to_string(d_points.size()));

      if (f_size < 8)
        return IntegratorExponential::simpson<T>(f, x);
      else {
        const double dx_17_p_48 = 17.0 * dx / 48.0;
        const double dx_59_p_48 = 59.0 * dx / 48.0;
        const double dx_43_p_48 = 43.0 * dx / 48.0;
        const double dx_49_p_48 = 49.0 * dx / 48.0;

        T pre = 0.0;
        pre += dx_17_p_48 * f[0] * d_points[0];
        pre += dx_59_p_48 * f[1] * d_points[1];
        pre += dx_43_p_48 * f[2] * d_points[2];
        pre += dx_49_p_48 * f[3] * d_points[3];

        auto i = 4u;
        T result = 0.0;
        while (i < f_size - 4) {
          result += dx * f[i] * d_points[i];
          ++i;
        }
        T post = 0.0;
        post += dx_49_p_48 * f[i] * d_points[i];
        post += dx_43_p_48 * f[i + 1] * d_points[i + 1];
        post += dx_59_p_48 * f[i + 2] * d_points[i + 2];
        post += dx_17_p_48 * f[i + 3] * d_points[i + 3];

        return pre + result + post;
      }
    }
  };

  class IntegratorGeneric
  {
  public:
    static double integrate(const std::vector<double> &f, const std::vector<double> &x);
  };
}

#endif //LATTICEGEN_INTEGRATOR_H
