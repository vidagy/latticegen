#ifndef LATTICEGEN_RADIALMESH_H
#define LATTICEGEN_RADIALMESH_H

#include <utility>
#include <vector>
#include <cmath>
#include <cstdlib>
#include <Core/ComparisonHelpers.h>
#include <Core/Exceptions.h>

namespace Core
{

  struct RadialMesh
  {
    virtual const std::vector<double> &get_points() const = 0;
  };

  struct GenericMesh : public RadialMesh
  {
    explicit GenericMesh(std::vector<double> points) : points(std::move(points)) {}

    const std::vector<double> &get_points() const override { return points; }

  private:
    std::vector<double> points;
  };

  struct ExponentialMesh : public RadialMesh
  {
    /// @brief Exponential mesh in the range of [a, b] having N points (endpoints included)
    ///
    /// x_i = multiplier * exp(dx * i) + shift
    /// where dx = ln(b/a) / (N-1) if scale == 1.0
    /// if scale is different, then dx = dx[scale = 1.0] * scale and
    /// the mesh is (linearly) rescaled to match endpoints
    ExponentialMesh(double a_, double b_, size_t N_, double scale_ = 0.7)
      : a(a_), b(b_), scale(scale_),
        dx(get_dx(a_, b_, N_, scale_)),
        multiplier((b - a) / (exp(dx * (N_ - 1.0)) - 1.0)),
        shift(a_ - multiplier),
        points(generate_points(multiplier, shift, dx, N_)),
        d_points(generate_d_points(multiplier, dx, N_)) {}

    double get_a() const { return a; }

    double get_b() const { return b; }

    double get_scale() const { return scale; }

    double get_dx() const { return dx; }

    double get_multiplier() const { return multiplier; }

    double get_shift() const { return shift; }

    const std::vector<double> &get_points() const override { return points; }

    const std::vector<double> &get_d_points() const { return d_points; }

  private:
    double a;
    double b;
    double scale;
    double dx;
    double multiplier;
    double shift;
    std::vector<double> points; /// multiplier_ * exp(dx_ * i) + shift_
    std::vector<double> d_points; /// (d points[i] / di) * dx =  multiplier_ * exp(dx_ * i)

    static double get_dx(double a_, double b_, size_t N_, double scale)
    {
      if (nearlyZero(a_))
        THROW_INVALID_ARGUMENT("a is nearly zero");
      if (nearlyZero(b_))
        THROW_INVALID_ARGUMENT("b is nearly zero");
      if (N_ < 2u)
        THROW_INVALID_ARGUMENT("N < 2");
      if (!strictlyPositive(scale))
        THROW_INVALID_ARGUMENT("scale is not strictly positive");

      return log(b_ / a_) / (N_ - 1) * scale;
    }

    static std::vector<double> generate_d_points(double multiplier_, double dx_, size_t N_)
    {
      std::vector<double> d_points;
      d_points.reserve(N_);

      for (auto i = 0u; i < N_; ++i) {
        d_points.push_back(multiplier_ * exp(dx_ * i));
      }
      return d_points;
    }

    static std::vector<double> generate_points(double multiplier_, double shift_, double dx_, size_t N_)
    {
      std::vector<double> points;
      points.reserve(N_);

      for (auto i = 0u; i < N_; ++i) {
        points.push_back(multiplier_ * exp(dx_ * i) + shift_);
      }
      return points;
    }
  };
}

#endif //LATTICEGEN_RADIALMESH_H
