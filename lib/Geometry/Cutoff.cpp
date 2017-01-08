#include "Cutoff.h"

#include <math.h>
#include <stdexcept>
#include <string>

#include <Core/ComparisonHelpers.h>
using namespace Core;

namespace Geometry
{
  static const Point3D x = Point3D(1.0, 0.0, 0.0);
  static const Point3D y = Point3D(1.0, 1.0, 0.0);
  static const Point3D z = Point3D(1.0, 0.0, 1.0);

  namespace
  {
    size_t get_max_steps(const Point3D& v, const Point3D& perp1, const Point3D& perp2, const double max_distance)
    {
      Point3D perp = cross_product(perp1, perp2);
      const double perp_length = perp.getLength();
      perp = 1.0 / perp_length * perp;
      
      double perp_component = v * perp; 
      return max_distance / perp_component + 1;
    }
  }

  CutoffCube::CutoffCube(const UnitCell3D& unit_cell_, const double a_)
    : Cutoff(unit_cell_), a(a_) 
  {
    if (! strictlyPositive(a))
      throw std::invalid_argument("In CutoffCube::ctor: a must be non-negative but a = " + std::to_string(a));
  }

  bool CutoffCube::is_included(const Point3D& point) const
  {
    return lessEqualsWithTolerance(fabs(point.x), a) && 
      lessEqualsWithTolerance(fabs(point.y), a) && 
      lessEqualsWithTolerance(fabs(point.z), a);
  }
  
  std::tuple<size_t, size_t, size_t> CutoffCube::steps_to_cover() const
  {
    return std::make_tuple(
        get_max_steps(unit_cell.a, unit_cell.b, unit_cell.c, a*sqrt(3)),
        get_max_steps(unit_cell.b, unit_cell.c, unit_cell.a, a*sqrt(3)),
        get_max_steps(unit_cell.c, unit_cell.a, unit_cell.b, a*sqrt(3))
      );
  }

  CutoffSphere::CutoffSphere(const UnitCell3D& unit_cell_, const double r_)
    : Cutoff(unit_cell_), r(r_) 
  {
    if (! strictlyPositive(r))
      throw std::invalid_argument("In CutoffSphere::ctor: a must be non-negative but r = " + std::to_string(r));
  }

  bool CutoffSphere::is_included(const Point3D& point) const
  {
    return lessEqualsWithTolerance(point.getLength(), r);
  }

  std::tuple<size_t, size_t, size_t> CutoffSphere::steps_to_cover() const
  {
    return std::make_tuple(
        get_max_steps(unit_cell.a, unit_cell.b, unit_cell.c, r),
        get_max_steps(unit_cell.b, unit_cell.c, unit_cell.a, r),
        get_max_steps(unit_cell.c, unit_cell.a, unit_cell.b, r)
      );
  }

  CutoffUnitVectors::CutoffUnitVectors(
    const UnitCell3D& unit_cell_,
    const size_t a_max_, const size_t b_max_, const size_t c_max_)
    : Cutoff(unit_cell_)
    , a_max(a_max_), b_max(b_max_), c_max(c_max_)
  {
  }

  bool CutoffUnitVectors::is_included(const Point3D& point) const
  {
    long nx, ny, nz;
    std::tie(nx, ny, nz) = unit_cell.get_offsets(point);
    return abs(nx) <= a_max &&
      abs(ny) <= b_max &&
      abs(nz) <= c_max;
  }
  
  std::tuple<size_t, size_t, size_t> CutoffUnitVectors::steps_to_cover() const
  {
    return std::make_tuple(a_max, b_max, c_max);
  }
}
