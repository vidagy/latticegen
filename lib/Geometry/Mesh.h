#ifndef LATTICEGEN_MESH_H
#define LATTICEGEN_MESH_H

#include <Core/Point3D.h>
#include <Math/CommonFunctions.h>
#include <Geometry/UnitCell3D.h>
#include <Geometry/Cutoff.h>
#include <vector>
#include <Math/SphericalHarmonics.h>

using namespace Core;

namespace Geometry
{
  class Mesh
  {
  public:
    std::vector<Point3D> generate(const Cutoff &cutoff, bool shift = false) const;

    const Cell3D cell;
  protected:
    explicit Mesh(const Cell3D &cell_) : cell(cell_) {}
  };

  class LatticeMesh : public Mesh
  {
  public:
    explicit LatticeMesh(const Cell3D &cell_) : Mesh(cell_)
    {}
  };

  class TrigonalMesh : public Mesh
  {
  public:
    explicit TrigonalMesh(double a_, double c_)
      : Mesh(UnitCell3D::create_hexagonal_primitive(a_, c_))
    {}
  };

  class TetrahedronMesh : public Mesh
  {
  public:
    explicit TetrahedronMesh(double a_)
      : Mesh(UnitCell3D::create_rhombohedral_centered(a_, pi / 3.0))
    {}
  };

  class CubicMesh : public Mesh
  {
  public:
    explicit CubicMesh(double a_)
      : Mesh(UnitCell3D::create_cubic_primitive(a_))
    {}
  };
}
#endif //LATTICEGEN_MESH_H
