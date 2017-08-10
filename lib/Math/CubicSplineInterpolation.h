#ifndef LATTICEGEN_INTERPOLATION_H
#define LATTICEGEN_INTERPOLATION_H

#include <vector>
#include <Core/Exceptions.h>
#include <algorithm>
#include <Eigen/Dense>

namespace Math
{
  class CubicSplineInterpolation
  {
  public:
    CubicSplineInterpolation(const std::vector<double> &knot_points_x_, const std::vector<double> &knot_points_y_);

    double interpolate(double x) const;

    std::vector<double> knot_points_x;
    std::vector<double> knot_points_y;
    Eigen::VectorXd k; // derivatives at knot points
  };
}


#endif //LATTICEGEN_INTERPOLATION_H
