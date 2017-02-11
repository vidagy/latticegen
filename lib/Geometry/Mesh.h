#ifndef LATTICEGEN_MESH_H
#define LATTICEGEN_MESH_H

#include <Geometry/Point3D.h>
#include <Geometry/UnitCell3D.h>
#include <Geometry/Cutoff.h>
#include <vector>

namespace Geometry
{
  class Mesh
  {
  public:
    virtual std::vector<Point3D> generate(const Cutoff& cutoff) const = 0;
  };

  class LatticeMesh : public Mesh
  {
  public:
    LatticeMesh(const UnitCell3D& unit_cell_, const bool positive_only_ = false)
      : unit_cell(unit_cell_)
      , positive_only(positive_only_)
    {}

    std::vector<Point3D> generate(const Cutoff& cutoff) const final override;

    UnitCell3D unit_cell;
    bool positive_only;
  };

  // class TetrahedronMesh : public Mesh

  // class CubicMesh : public Mesh
}
#endif //LATTICEGEN_MESH_H
