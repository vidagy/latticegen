#include <algorithm>
#include "Mesh.h"

namespace Geometry
{
  std::vector<Point3D> LatticeMesh::generate(const Cutoff& cutoff) const
  {
    Cutoff::StepsToCover steps_to_cover = cutoff.steps_to_cover(unit_cell);

    std::vector<Point3D> lattice;
    lattice.reserve((unsigned long) (
      (steps_to_cover.x_max - steps_to_cover.x_min) *
      (steps_to_cover.y_max - steps_to_cover.y_min) *
      (steps_to_cover.z_max - steps_to_cover.z_min)
    ));

    for (long n_z = steps_to_cover.z_min; n_z <= steps_to_cover.z_max; ++n_z)
    {
      auto c = n_z * unit_cell.c;
      for (long n_y = steps_to_cover.y_min; n_y <= steps_to_cover.y_max; ++n_y)
      {
        auto b = n_y * unit_cell.b;
        for (long n_x = steps_to_cover.x_min; n_x <= steps_to_cover.x_max; ++n_x)
        {
          auto a = n_x * unit_cell.a;
          auto v = a + b + c;
          if (cutoff.is_included(v))
            lattice.push_back(v);
        }
      }
    }
    lattice.shrink_to_fit();
    return lattice;
  }
}
