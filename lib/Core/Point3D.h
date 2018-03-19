#ifndef LATTICEGEN_POINT3D_H_
#define LATTICEGEN_POINT3D_H_

#include <limits>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <string>
#include <cmath>
#include "Mute.h"

LATTICEGEN_MUTE_BEGIN
LATTICEGEN_MUTE_EIGEN
#include <Eigen/Dense>

LATTICEGEN_MUTE_END

#include <Core/ComparisonHelpers.h>

namespace Core
{
  typedef Eigen::Vector3d Point3D;
  typedef Eigen::Ref<Eigen::Vector3d> Point3DRef;
  typedef Eigen::Ref<const Eigen::Vector3d> Point3DCRef;

  [[maybe_unused]]
  static Point3D create_polar(double r, double theta, double phi) {
    return Point3D{r * sin(theta) * cos(phi), r * sin(theta) * sin(phi), r * cos(theta)};
  }

  struct Point3DComparator
  {
    explicit Point3DComparator(double abs_tol, double rel_tol)
      : x_abs_tol(abs_tol), y_abs_tol(abs_tol), z_abs_tol(abs_tol), x_rel_tol(rel_tol), y_rel_tol(rel_tol),
        z_rel_tol(rel_tol) {}

    Point3DComparator(
      double x_abs_tol_, double y_abs_tol_, double z_abs_tol_,
      double x_rel_tol_, double y_rel_tol_, double z_rel_tol_)
      : x_abs_tol(x_abs_tol_), y_abs_tol(y_abs_tol_), z_abs_tol(z_abs_tol_), x_rel_tol(x_rel_tol_),
        y_rel_tol(y_rel_tol_), z_rel_tol(z_rel_tol_) {}

    bool isEqual(const Point3DCRef &lhs, const Point3DCRef &rhs) const
    {
      return equalsWithTolerance(lhs(0), rhs(0), x_abs_tol, x_rel_tol) &&
             equalsWithTolerance(lhs(1), rhs(1), y_abs_tol, y_rel_tol) &&
             equalsWithTolerance(lhs(2), rhs(2), z_abs_tol, z_rel_tol);
    }

    const double x_abs_tol;
    const double y_abs_tol;
    const double z_abs_tol;
    const double x_rel_tol;
    const double y_rel_tol;
    const double z_rel_tol;
  };

  static const auto default_point_comparator =
    Point3DComparator(Core::default_abs_tolerance, Core::default_rel_tolerance);
}

namespace std
{
  template<typename Derived>
  std::string to_string(const Eigen::MatrixBase<Derived> &m)
  {
    std::stringstream ss;
    ss << std::scientific << std::setprecision(std::numeric_limits<double>::max_digits10) << m;
    return ss.str();
  }
}
#endif // LATTICEGEN_POINT3D_H_
