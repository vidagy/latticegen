#include <algorithm>
#include <Core/Exceptions.h>
#include "IrreducibleWedge.h"
#include "Mesh.h"

using namespace Geometry;

namespace
{
  std::vector<Point3D> get_irreduc_points(
    std::vector<Point3D>::const_iterator begin,
    const std::vector<Point3D>::const_iterator &end,
    const SymmetryTransformationFactory::Transformations &transformations
  )
  {
    if (begin == end)
      THROW_INVALID_ARGUMENT("empty range");
    const auto size = static_cast<unsigned long>(std::distance(begin, end));
    std::vector<bool> already_covered(size, false);
    std::vector<Point3D> result;
    result.reserve(size);
    auto still_false = size;

    for (size_t i = 0; begin != end; ++i, ++begin) {
      // for all points in shell
      if (!already_covered[i]) {
        result.push_back(*begin);

        // account all points that are replicated by the transformation
        for (auto transformation: transformations) {
          auto transformed = transformation * (*begin);
          if (!(transformed == (*begin))) {
            for (size_t j = i + 1; j < size; ++j) {
              if ((!already_covered[j]) && (*(begin + j - i) == transformed)) {
                already_covered[j] = true;
                --still_false;
                if (still_false == 0)
                  return result;
              }
            }
          }
        }
      }
      already_covered[i] = true;
      --still_false;
      if (still_false == 0)
        return result;
    }
    return result;
  }
}

std::vector<Point3D>
IrreducibleWedge::reduce_by_symmetries(std::vector<Point3D> points,
                                       const SymmetryTransformationFactory::Transformations &transformations)
{
  std::vector<Point3D> result;
  result.reserve(points.size() / transformations.size() * 2);

  std::sort(points.begin(), points.end(),
            [](const Point3D &lhs, const Point3D &rhs)
            {
              return lhs.length() < rhs.length();
            });

  auto pre_length = points[0].length();
  auto begin = points.begin();
  for (auto it = points.begin(); it != points.end(); ++it) {
    if (!equalsWithTolerance(pre_length, it->length())) {
      auto current_irreduc_points = ::get_irreduc_points(begin, it, transformations);
      std::copy(current_irreduc_points.begin(), current_irreduc_points.end(), std::back_inserter(result));
      pre_length = it->length();
      begin = it;
    }
  }
  auto current_irreduc_points = ::get_irreduc_points(begin, points.end(), transformations);
  std::copy(current_irreduc_points.begin(), current_irreduc_points.end(), std::back_inserter(result));

  return result;
}

std::vector<Point3D>
IrreducibleWedge::get_irreducible_wedge(const Cell3D &cell, size_t sample)
{
  if (sample == 0)
    THROW_INVALID_ARGUMENT("sample is zero in IrreducibleWedge::get_irreducible_wedge");

  auto transformations = SymmetryTransformationFactory::generate(cell);

  std::unique_ptr<Mesh> mesh;
  CrystalSystem crystal_system = get_crystal_system(cell.get_point_group());
  if (crystal_system == CrystalSystem::Hexagonal)
  {
    double sample_width = cell.v1.length() / sample;
    mesh = std::make_unique<TetrahedronMesh>(sample_width);
  }
  else if (crystal_system == CrystalSystem::Trigonal)
  {
    double sample_width_a = cell.v1.length() / sample;
    double sample_width_c = cell.v3.length() / sample;
    mesh = std::make_unique<TrigonalMesh>(sample_width_a, sample_width_c);
  }
  else
  {
    double sample_width = cell.v1.length() / sample;
    mesh = std::make_unique<CubicMesh>(sample_width);
  }

  return reduce_by_symmetries(mesh->generate(CutoffWSCell(cell)), transformations);
}
