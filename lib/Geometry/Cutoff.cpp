#include <algorithm>
#include <Core/Mute.h>

LATTICEGEN_MUTE_BEGIN
LATTICEGEN_MUTE_EIGEN
#include <Eigen/Dense>

LATTICEGEN_MUTE_END

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
  std::vector<std::pair<Point3D, double>> get_all_face_points(const Cell3D &cell)
  {
    const Point3D a = 0.5 * cell.v1;
    const Point3D b = 0.5 * cell.v2;
    const Point3D c = 0.5 * cell.v3;

    // we count on inversion symmetry so we can store one vector in each direction
    const auto face_points = std::vector<Point3D>{
      // note that the order of these points matter, check get_face_points() function
      a, b, c,

      a + b, a - b,
      a + c, a - c,
      b + c, b - c,

      a + b + c, a - b + c, a + b - c, -a + b + c
    };
    auto res = std::vector<std::pair<Point3D, double>>();
    res.reserve(face_points.size());
    std::transform(
      face_points.begin(), face_points.end(), std::back_inserter(res),
      [](const Point3D &p)
      {
        auto len = p.length();
        return std::make_pair(p / len, len);
      }
    );

    return res;
  }

  std::vector<std::pair<Point3D, double>> get_face_points(const Cell3D &cell)
  {
    auto face_points = get_all_face_points(cell);

    std::sort(
      face_points.begin(), face_points.end(),
      [](const std::pair<Point3D, double> &lhs, const std::pair<Point3D, double> &rhs) -> bool
      {
        return lhs.second < rhs.second;
      }
    );

    if (face_points.size() < 3u)
      THROW_LOGIC_ERROR("face_points.size() = " + std::to_string(face_points.size()) + " < 3");
    if (face_points.size() == 3u)
      return face_points;

    // we have more than 3 face points, we can potentially remove a few
    for (auto it1 = 0u; it1 < face_points.size() - 2; ++it1) {
      for (auto it2 = it1 + 1u; it2 < face_points.size() - 1; ++it2) {
        for (auto it3 = it2 + 1u; it3 < face_points.size(); ++it3) {
          const auto a = face_points[it1];
          const auto b = face_points[it2];
          const auto c = face_points[it3];
          if (strictlyPositive(fabs(cross_product(a.first, b.first) * c.first))) {
            face_points.erase(
              std::remove_if(
                face_points.begin(), face_points.end(),
                [&a, &b, &c](const std::pair<Point3D, double> &p)
                {
                  auto outside_of_parallelepiped =
                    greaterEqualsWithTolerance(p.second * fabs(p.first * a.first), a.second) ||
                    greaterEqualsWithTolerance(p.second * fabs(p.first * b.first), b.second) ||
                    greaterEqualsWithTolerance(p.second * fabs(p.first * c.first), c.second);
                  auto is_current_point =
                    (a.first == p.first && a.second == p.second) ||
                    (b.first == p.first && b.second == p.second) ||
                    (c.first == p.first && c.second == p.second);
                  return outside_of_parallelepiped && !is_current_point;
                }
              ), face_points.end()
            );
          }
        }
      }
    }

    face_points.shrink_to_fit();
    return face_points;
  }

  double get_min_r(const std::vector<std::pair<Point3D, double>> &points)
  {
    auto min = std::min_element(
      points.begin(), points.end(),
      [](const std::pair<Point3D, double> &lhs, const std::pair<Point3D, double> &rhs) -> bool
      {
        return lhs.second < rhs.second;
      });
    return min->second;
  }

  double get_max_r(const std::vector<std::pair<Point3D, double>> &face_points)
  {
    // we have to include the inversion of all face points, too
    auto all_face_points = std::vector<std::pair<Point3D, double>>();
    all_face_points.reserve(face_points.size() * 2);
    for (const auto &p : face_points) {
      all_face_points.push_back(p);
      all_face_points.emplace_back(-1.0 * p.first, p.second);
    }

    auto possible_distances = std::vector<double>();
    possible_distances.reserve(all_face_points.size());
    for (auto i = 0u; i < all_face_points.size(); ++i) {
      for (auto j = i + 1; j < all_face_points.size(); ++j) {
        for (auto k = j + 1; k < all_face_points.size(); ++k) {
          auto a1 = all_face_points[i];
          auto a2 = all_face_points[j];
          auto a3 = all_face_points[k];

          if (!nearlyZero(fabs(cross_product(a1.first, a2.first) * a3.first))) { // non co-planar
            Eigen::Matrix3d m(3, 3);
            m <<
              a1.first.x, a1.first.y, a1.first.z,
              a2.first.x, a2.first.y, a2.first.z,
              a3.first.x, a3.first.y, a3.first.z;

            Eigen::Vector3d v(a1.second, a2.second, a3.second);
            Eigen::Vector3d res = m.colPivHouseholderQr().solve(v);

            // we only add it to the intersections if the point is on the surface
            auto is_included = true;
            for (const auto &p : face_points) {
              const auto scalar_prod = res(0) * p.first.x + res(1) * p.first.y + res(2) * p.first.z;
              if (strictlyGreater(fabs(scalar_prod), p.second)) {
                is_included = false;
                break;
              }
            }
            if (is_included)
              possible_distances.push_back(res.norm());
          }
        }
      }
    }

    auto max = std::max_element(possible_distances.begin(), possible_distances.end());
    return *max;
  }
}

CutoffWSCell::CutoffWSCell(const Cell3D &cell_)
  : cell(cell_),
    face_points(get_face_points(cell_)),
    r_mt(get_min_r(face_points)),
    r_bs(get_max_r(face_points)) {}

bool CutoffWSCell::is_included(const Point3D &point) const
{
  auto length2 = point.length2();
  if (strictlyLess(length2, r_mt * r_mt))
    return true;
  if (strictlyGreater(length2, r_bs * r_bs))
    return false;
  for (const auto &face_point: face_points) {
    if (strictlyGreater(fabs(point * face_point.first), face_point.second))
      return false;
  }
  return true;
}

Cutoff::StepsToCover CutoffWSCell::steps_to_cover(const Cell3D &cell_) const
{
  return get_steps_to_cover(cell_, r_bs);
}
