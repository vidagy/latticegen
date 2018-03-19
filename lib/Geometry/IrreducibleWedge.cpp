#include <algorithm>
#include <Core/Exceptions.h>
#include "IrreducibleWedge.h"
#include "Mesh.h"

using namespace Geometry;

namespace
{
  std::vector<Point3D> get_irreduc_points(
    std::vector<Point3D>::iterator &begin,
    const std::vector<Point3D>::iterator &end,
    SymmetryTransformationFactory::Transformations &transformations,
    const double abs_tolerance
  )
  {
    if (begin == end)
      THROW_INVALID_ARGUMENT("empty range");
    const auto size = static_cast<unsigned long>(std::distance(begin, end));
    if (size == 1u)
      return std::vector<Point3D>(begin, end);

    std::vector<bool> already_covered(size, false);
    std::vector<Point3D> result;
    result.reserve(size);
    auto still_false = size;

    for (size_t i = 0; begin != end; ++i, ++begin) {
      // for all points in shell
      if (!already_covered[i]) {
        result.push_back(*begin);
        already_covered[i] = true;
        --still_false;

        // account all points that are replicated by the transformation
        for (auto &transformation: transformations) {
          if (transformation.type != Transformation::Type::Identity) {
          auto transformed = transformation * (*begin);
            for (size_t j = i + 1; j < size; ++j) {
              if ((!already_covered[j]) && (begin - i + j)->isApprox(transformed, abs_tolerance)) {
                already_covered[j] = true;
                --still_false;
              }
            }
          }
        }
      }
      if (still_false == 0)
        return result;
    }
    return result;
  }
}

std::vector<Point3D>
IrreducibleWedge::reduce_by_symmetries(
  std::vector<Point3D> points,
  SymmetryTransformationFactory::Transformations &transformations,
  double abs_tolerance
)
{
  if (transformations.empty())
    return points;

  std::vector<Point3D> result;
  result.reserve(points.size() / transformations.size() * 2);

  std::sort(points.begin(), points.end(),
            [](const Point3D &lhs, const Point3D &rhs)
            {
              return lhs.norm() < rhs.norm();
            });

  auto pre_length = points[0].norm();
  auto begin = points.begin();
  for (auto it = points.begin(); it != points.end(); ++it) {
    if (!equalsWithTolerance(pre_length, it->norm(), 1e-10)) {
      auto current_irreduc_points = ::get_irreduc_points(begin, it, transformations, abs_tolerance);
      std::copy(current_irreduc_points.begin(), current_irreduc_points.end(), std::back_inserter(result));
      pre_length = it->norm();
      begin = it;
    }
  }
  auto current_irreduc_points = ::get_irreduc_points(begin, points.end(), transformations, abs_tolerance);
  std::copy(current_irreduc_points.begin(), current_irreduc_points.end(), std::back_inserter(result));

  return result;
}

double IrreducibleWedge::get_tolerance(const Cell3D &cell)
{
  auto v1v2 = cell.v1.cross(cell.v2);
  auto v2v3 = cell.v2.cross(cell.v3);
  auto v1v3 = cell.v1.cross(cell.v3);


  auto v1_tol = 0.5 * cell.v1.dot(v2v3) / v2v3.norm();
  auto v2_tol = 0.5 * cell.v2.dot(v1v3) / v1v3.norm();
  auto v3_tol = 0.5 * cell.v3.dot(v1v2) / v1v2.norm();

  auto abs_tolerance = std::min(std::min(fabs(v1_tol), fabs(v2_tol)), fabs(v3_tol));
  return abs_tolerance;
}

std::vector<Point3D>
IrreducibleWedge::get_irreducible_wedge(const Cell3D &cell, size_t sample, bool centered)
{
  if (sample == 0)
    THROW_INVALID_ARGUMENT("sample is zero in IrreducibleWedge::get_irreducible_wedge");

  auto transformations = SymmetryTransformationFactory::generate(cell);

  std::unique_ptr<Mesh> mesh;
  CrystalSystem crystal_system = get_crystal_system(cell.get_point_group());
  if (crystal_system == CrystalSystem::Hexagonal)
  {
    double sample_width = cell.v1.norm() / sample;
    mesh = std::make_unique<TetrahedronMesh>(sample_width);
  }
  else if (crystal_system == CrystalSystem::Trigonal)
  {
    double sample_width_a = cell.v1.norm() / sample;
    double sample_width_c = cell.v3.norm() / sample;
    mesh = std::make_unique<TrigonalMesh>(sample_width_a, sample_width_c);
  }
  else
  {
    double sample_width = cell.v1.norm() / sample;
    mesh = std::make_unique<CubicMesh>(sample_width);
  }

  double abs_tol = get_tolerance(mesh->cell);
  return reduce_by_symmetries(mesh->generate(CutoffWSCell(cell), centered), transformations, abs_tol);
}
