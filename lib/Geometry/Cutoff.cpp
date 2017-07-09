#include <algorithm>
#include "Cutoff.h"

using namespace Geometry;

namespace
{
  std::tuple<long, long, long> get_max_steps_and_offsets(
    const Point3D &v, const Point3D &perpendicular1, const Point3D &perpendicular2, double max_distance
  )
  {
    Point3D perpendicular = cross_product(perpendicular1, perpendicular2);
    const double perpendicular_length = perpendicular.length();
    perpendicular = 1.0 / perpendicular_length * perpendicular;
    double perpendicular_component = v * perpendicular;

    long max_steps = long(max_distance / perpendicular_component) + 1;
    long offset1 = long((perpendicular1 * v * max_steps) / perpendicular1.length()) + 1;
    long offset2 = long((perpendicular2 * v * max_steps) / perpendicular2.length()) + 1;

    return std::make_tuple(max_steps, offset1, offset2);
  }

  Cutoff::StepsToCover get_steps_to_cover(
    const Cell3D &cell, double max_distance
  )
  {
    long a_max_steps, a_offset1, a_offset2, b_max_steps, b_offset1, b_offset2, c_max_steps, c_offset1, c_offset2;
    const Vector3D &a = cell.v1;
    const Vector3D &b = cell.v2;
    const Vector3D &c = cell.v3;

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
  if (!strictlyPositive(a))
    THROW_INVALID_ARGUMENT("In CutoffCube::ctor: a must be non-negative but a = " + std::to_string(a));
}

bool CutoffCube::is_included(const Point3D &point) const
{
  return lessEqualsWithTolerance(fabs(point.x), a) &&
         lessEqualsWithTolerance(fabs(point.y), a) &&
         lessEqualsWithTolerance(fabs(point.z), a);
}

Cutoff::StepsToCover CutoffCube::steps_to_cover(const Cell3D &cell) const
{
  return get_steps_to_cover(cell, a * sqrt(3));
}

CutoffSphere::CutoffSphere(double r_)
  : r(r_)
{
  if (!strictlyPositive(r))
    THROW_INVALID_ARGUMENT("In CutoffSphere::ctor: a must be non-negative but r = " + std::to_string(r));
}

bool CutoffSphere::is_included(const Point3D &point) const
{
  return lessEqualsWithTolerance(point.length(), r);
}

Cutoff::StepsToCover CutoffSphere::steps_to_cover(const Cell3D &cell) const
{
  return get_steps_to_cover(cell, r);
}

CutoffUnitVectors::CutoffUnitVectors(
  const Cell3D &cell_,
  size_t a_max_, size_t b_max_, size_t c_max_, bool positive_only_)
  : cell(cell_), a_max(a_max_), b_max(b_max_), c_max(c_max_), positive_only(positive_only_)
{
}

bool CutoffUnitVectors::is_included(const Point3D &point) const
{
  long nx, ny, nz;
  std::tie(nx, ny, nz) = cell.get_offsets(point);
  return (size_t) labs(nx) <= a_max &&
         (size_t) labs(ny) <= b_max &&
         (size_t) labs(nz) <= c_max;
}

Cutoff::StepsToCover CutoffUnitVectors::steps_to_cover(const Cell3D &cell) const
{
  if (positive_only)
    return {0l, 0l, 0l, (long) a_max, (long) b_max, (long) c_max};
  else
    return {-(long) a_max, -(long) b_max, -(long) c_max, (long) a_max, (long) b_max, (long) c_max};
}

namespace
{
  std::vector<Point3D> get_neighbors(const Cell3D &cell)
  {
    const Point3D &a = cell.v1;
    const Point3D &b = cell.v2;
    const Point3D &c = cell.v3;

    auto neighbors = std::vector<Point3D>{
      a, -a, b, -b, c, -c,

      a + b, a - b, -a + b, -a - b,
      a + c, a - c, -a + c, -a - c,
      b + c, b - c, -b + c, -b - c,

      a + b + c, a - b + c, a + b - c, -a + b + c,
      a - b - c, -a - b + c, -a + b - c, -a - b - c
    };
    return neighbors;
  }

  // TODO this is incorrect! min_r can be smaller than the half distance to the closest neighbor.
  double get_min_r(const std::vector<Point3D> &points)
  {
    auto min = std::min_element(
      points.begin(), points.end(),
      [](const Point3D &lhs, const Point3D &rhs)
      {
        return lhs.length2() < rhs.length2();
      });
    return 0.5 * min->length();
  }

  // TODO this is incorrect! max_r can be larger than the half distance to the furthest neighbor.
  double get_max_r(const std::vector<Point3D> &points)
  {
    auto max = std::max_element(
      points.begin(), points.end(),
      [](const Point3D &lhs, const Point3D &rhs)
      {
        return lhs.length2() < rhs.length2();
      });
    return 0.5 * max->length();
  }
}

CutoffWSCell::CutoffWSCell(const Cell3D &cell_)
  : cell(cell_), neighbors(get_neighbors(cell_)), r_in_for_sure(get_min_r(neighbors)),
    r_out_for_sure(get_max_r(neighbors)) {}

bool CutoffWSCell::is_included(const Point3D &point) const
{
  auto length2 = point.length2();
  if (strictlyLess(length2, r_in_for_sure * r_in_for_sure))
    return true;
  if (strictlyGreater(length2, r_out_for_sure * r_out_for_sure))
    return false;
  // TODO this should be optimized! we should not iterate over all the neighbors, only those that give faces to the
  // Voronoi cell! (could speed up is_included ~35% in best case scenario)
  for (const auto &neighbor: neighbors) {
    if (length2 > (point - neighbor).length2())
      return false;
  }
  return true;
}

Cutoff::StepsToCover CutoffWSCell::steps_to_cover(const Cell3D &cell_) const
{
  double r = 0.0;
  for (const auto &neighbor: neighbors) {
    auto length = neighbor.length();
    if (length > r)
      r = length;
  }
  return get_steps_to_cover(cell_, r);
}
