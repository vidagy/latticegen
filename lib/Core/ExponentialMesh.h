#ifndef LATTICEGEN_EXPONENTIALMESH_H
#define LATTICEGEN_EXPONENTIALMESH_H

#include <vector>

namespace Core
{
  struct ExponentialMesh
  {
    // Exponential mesh in the range of [r0 dx, b]
    // dx is the logarithmic increment: (p[i+1]-p[i])/p[i] = dx

    ExponentialMesh(double b_, double dx_, size_t num_of_points_)
      : r0(get_r0(b_, dx_, num_of_points_)), dx(dx_), points(generate(r0, dx, num_of_points_)) {}

    const double r0;
    const double dx;
    const std::vector<double> points;

  private:
    static double get_r0(double b_, double dx_, size_t num_of_points_)
    {
      return b_ / (exp(dx_ * num_of_points_) - 1.0);
    }

    static std::vector<double> generate(double r0_, double dx_, size_t num_of_points_)
    {
      std::vector<double> points;
      points.reserve(num_of_points_);
      for (auto i = 1u; i <= num_of_points_; ++i) {
        points.push_back(r0_ * (exp(dx_ * i) - 1.0));
      }
      return points;
    }
  };
}

#endif //LATTICEGEN_EXPONENTIALMESH_H
