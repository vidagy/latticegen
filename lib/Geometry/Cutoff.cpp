#include "Cutoff.h"

namespace Geometry
{
  namespace
  {
    std::tuple< long, long, long > get_max_steps_and_offsets(
      const Point3D &v, const Point3D &perpendicular1, const Point3D &perpendicular2, double max_distance
    )
    {
      Point3D perpendicular = cross_product(perpendicular1, perpendicular2);
      const double perpendicular_length = perpendicular.length();
      perpendicular = 1.0 / perpendicular_length * perpendicular;
      double perpendicular_component = v * perpendicular;

      long max_steps = long(max_distance / perpendicular_component) + 1;
      long offset1 = long( (perpendicular1*v*max_steps) / perpendicular1.length() ) + 1;
      long offset2 = long( (perpendicular2*v*max_steps) / perpendicular2.length() ) + 1;

      return std::make_tuple(max_steps, offset1, offset2);
    }

    Cutoff::StepsToCover get_steps_to_cover(
      const UnitCell3D &unit_cell, double max_distance
    )
    {
      long a_max_steps, a_offset1, a_offset2, b_max_steps, b_offset1, b_offset2, c_max_steps, c_offset1, c_offset2;
      const Vector3D& a = unit_cell.a; const Vector3D& b = unit_cell.b; const Vector3D& c = unit_cell.c;

      std::tie(a_max_steps, b_offset1, c_offset1) = get_max_steps_and_offsets(a, b, c, max_distance);
      std::tie(b_max_steps, c_offset2, a_offset1) = get_max_steps_and_offsets(b, c, a, max_distance);
      std::tie(c_max_steps, a_offset2, b_offset2) = get_max_steps_and_offsets(c, a, b, max_distance);

      long res_a = a_max_steps + std::max(a_offset1, a_offset2);
      long res_b = b_max_steps + std::max(b_offset1, b_offset2);
      long res_c = c_max_steps + std::max(c_offset1, c_offset2);

      return {-res_a, -res_b, -res_c, res_a, res_b, res_c};
    }
  }

  CutoffCube::CutoffCube(double a_)
    : a(a_)
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

  Cutoff::StepsToCover CutoffCube::steps_to_cover(const UnitCell3D& unit_cell) const
  {
    return get_steps_to_cover(unit_cell, a * sqrt(3));
  }

  CutoffSphere::CutoffSphere(double r_)
    : r(r_)
  {
    if (! strictlyPositive(r))
      throw std::invalid_argument("In CutoffSphere::ctor: a must be non-negative but r = " + std::to_string(r));
  }

  bool CutoffSphere::is_included(const Point3D& point) const
  {
    return lessEqualsWithTolerance(point.length(), r);
  }

  Cutoff::StepsToCover CutoffSphere::steps_to_cover(const UnitCell3D& unit_cell) const
  {
    return get_steps_to_cover(unit_cell, r);
  }

  CutoffUnitVectors::CutoffUnitVectors(
    const UnitCell3D& unit_cell_,
    size_t a_max_, size_t b_max_, size_t c_max_, bool positive_only_)
    : unit_cell(unit_cell_)
    , a_max(a_max_), b_max(b_max_), c_max(c_max_)
    , positive_only(positive_only_)
  {
  }

  bool CutoffUnitVectors::is_included(const Point3D& point) const
  {
    long nx, ny, nz;
    std::tie(nx, ny, nz) = unit_cell.get_offsets(point);
    return (size_t)labs(nx) <= a_max &&
           (size_t)labs(ny) <= b_max &&
           (size_t)labs(nz) <= c_max;
  }

  Cutoff::StepsToCover CutoffUnitVectors::steps_to_cover(const UnitCell3D& unit_cell) const
  {
    if (positive_only)
      return {0l, 0l, 0l, (long)a_max, (long)b_max, (long)c_max};
    else
      return {-(long)a_max, -(long)b_max, -(long)c_max, (long)a_max, (long)b_max, (long)c_max};
  }

  CutoffWSCell::CutoffWSCell(const UnitCell3D& unit_cell_)
  : unit_cell(unit_cell_)
  {
    const Point3D& a = unit_cell_.a;
    const Point3D& b = unit_cell_.b;
    const Point3D& c = unit_cell_.c;

    neighbors = {
      a,-a,b,-b,c,-c,

      a+b,a-b,-a+b,-a-b,
      a+c,a-c,-a+c,-a-c,
      b+c,b-c,-b+c,-b-c,

      a+b+c,a-b+c,a+b-c,-a+b+c,
      a-b-c,-a-b+c,-a+b-c,-a-b-c,
    };
  }

  bool CutoffWSCell::is_included(const Point3D& point) const
  {
    auto length = point.length();
    for (auto neighbor: neighbors)
    {
      if (length > (point-neighbor).length())
        return false;
    }
    return true;
  }

  Cutoff::StepsToCover CutoffWSCell::steps_to_cover(const UnitCell3D& unit_cell_) const
  {
    double r = 0.0;
    for (auto neighbor: neighbors)
    {
      auto length = neighbor.length();
      if (length > r)
        r = length;
    }
    return get_steps_to_cover(unit_cell_, r);
  }
}