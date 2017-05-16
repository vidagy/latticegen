#include <algorithm>
#include <Core/Exceptions.h>
#include "IrreducibleWedge.h"
#include "Mesh.h"

using namespace Geometry;

std::vector<Point3D>
IrreducibleWedge::reduce_by_symmetries(const std::vector<Point3D> &points,
                                       const SymmetryTransformationFactory::Transformations &transformations)
{
  std::vector<Point3D> result;
  result.reserve(points.size());
  std::vector<bool> already_covered(points.size(), false);
  for (size_t i = 0; i < points.size(); ++i)
  {
    if (! already_covered[i])
    {
      result.push_back(points[i]);

      for (size_t j = i+1; j < points.size(); ++j)
      {
        if (! already_covered[j])
        {
          for (auto transformation: transformations)
          {
            if (points[j] == transformation*points[i])
            {
              already_covered[j] = true;
              break;
            }
          }
        }
      }
    }
    already_covered[i] = true;
  }
  return result;
}

std::vector<Point3D>
IrreducibleWedge::get_irreducible_wedge(const UnitCell3D &unit_cell, size_t sample)
{
  if (sample == 0)
    THROW_INVALID_ARGUMENT("sample is zero in IrreducibleWedge::get_irreducible_wedge");

  std::unique_ptr<CrystallographicPointGroup> group = CrystallographicPointGroup::create(unit_cell.get_point_group());
  SymmetryTransformationFactory::Transformations transformations
    = SymmetryTransformationFactory::get(group->get_elements());

  std::unique_ptr<Mesh> mesh;
  CrystalSystem crystal_system = get_crystal_system(unit_cell.get_point_group());
  if (crystal_system == CrystalSystem::Hexagonal)
  {
    double sample_width = unit_cell.v1.length() / sample;
    mesh = std::make_unique<TetrahedronMesh>(sample_width);
  }
  else if (crystal_system == CrystalSystem::Trigonal)
  {
    double sample_width_a = unit_cell.v1.length() / sample;
    double sample_width_c = unit_cell.v3.length() / sample;
    mesh = std::make_unique<TrigonalMesh>(sample_width_a, sample_width_c);
  }
  else
  {
    double sample_width = unit_cell.v1.length() / sample;
    mesh = std::make_unique<CubicMesh>(sample_width);
  }

  return reduce_by_symmetries(mesh->generate(CutoffWSCell(unit_cell)), transformations);
}
