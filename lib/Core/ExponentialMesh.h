#ifndef LATTICEGEN_EXPONENTIALMESH_H
#define LATTICEGEN_EXPONENTIALMESH_H

#include <vector>
#include <cmath>
#include <cstdlib>
#include <Core/ComparisonHelpers.h>
#include <Core/Exceptions.h>

namespace Core
{
  struct ExponentialMesh
  {
    // Exponential mesh in the range of [a, b] having N points (endpoints included)
    // x_i = a * exp(dx * i) where dx = ln(b/a) / (N-1)

    ExponentialMesh(double a_, double b_, size_t N_)
      : a(a_), dx(get_dx(a_, b_, N_)), points(generate(a, dx, N_)) {}

    const double a;
    const double dx;
    const std::vector<double> points;

  private:
    static double get_dx(double a_, double b_, size_t N_)
    {
      if (nearlyZero(a_))
        THROW_INVALID_ARGUMENT("in ExponentialMesh a is zero");
      if (N_ < 2u)
        THROW_INVALID_ARGUMENT("in ExponentialMesh N < 2");

      return log(b_ / a_) / (N_ - 1);
    }

    static std::vector<double> generate(double a_, double dx_, size_t N_)
    {
      std::vector<double> points;
      points.reserve(N_);
      for (auto i = 0u; i < N_; ++i) {
        points.push_back(a_ * exp(dx_ * i));
      }
      return points;
    }
  };
}

#endif //LATTICEGEN_EXPONENTIALMESH_H
