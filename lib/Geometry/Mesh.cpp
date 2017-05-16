#include <algorithm>
#include "Mesh.h"

using namespace Geometry;

namespace
{
  static const double pi = 3.14159265358979323846;

  std::vector<Point3D> generate_mesh(
    const Cell3D &cell3D, const Cutoff &cutoff, const Point3D &offset)
  {
    Cutoff::StepsToCover steps_to_cover = cutoff.steps_to_cover(cell3D);

    std::vector<Point3D> lattice;
    lattice.reserve((unsigned long) (
      (steps_to_cover.x_max - steps_to_cover.x_min) *
      (steps_to_cover.y_max - steps_to_cover.y_min) *
      (steps_to_cover.z_max - steps_to_cover.z_min)
    ));

    for (long n_z = steps_to_cover.z_min; n_z <= steps_to_cover.z_max; ++n_z) {
      auto c = n_z * cell3D.v3 + offset;
      for (long n_y = steps_to_cover.y_min; n_y <= steps_to_cover.y_max; ++n_y) {
        auto b = n_y * cell3D.v2;
        for (long n_x = steps_to_cover.x_min; n_x <= steps_to_cover.x_max; ++n_x) {
          auto a = n_x * cell3D.v1;
          auto v = a + b + c;
          if (cutoff.is_included(v))
            lattice.push_back(v);
        }
      }
    }
    lattice.shrink_to_fit();
    return lattice;
  }

  Point3D get_offset(const Cell3D &unit_cell, bool shift)
  {
    return shift ? (unit_cell.v1 + unit_cell.v2 + unit_cell.v3) / 2.0 : Point3D();
  }
}

std::vector<Point3D> LatticeMesh::generate(const Cutoff &cutoff, bool shift) const
{
  return generate_mesh(cell, cutoff, get_offset(cell, shift));
}

std::vector<Point3D> TrigonalMesh::generate(const Cutoff &cutoff, bool shift) const
{
  const UnitCell3D unit_cell = UnitCell3D::create_hexagonal_primitive(a, c);
  return generate_mesh(unit_cell, cutoff, get_offset(unit_cell, shift));
}

std::vector<Point3D> TetrahedronMesh::generate(const Cutoff &cutoff, bool shift) const
{
  const UnitCell3D unit_cell = UnitCell3D::create_rhombohedral_centered(a, pi / 3.0);
  return generate_mesh(unit_cell, cutoff, get_offset(unit_cell, shift));
}

std::vector<Point3D> CubicMesh::generate(const Cutoff &cutoff, bool shift) const
{
  const UnitCell3D unit_cell = UnitCell3D::create_cubic_primitive(a);
  return generate_mesh(unit_cell, cutoff, get_offset(unit_cell, shift));
}
