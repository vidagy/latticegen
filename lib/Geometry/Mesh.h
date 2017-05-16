#ifndef LATTICEGEN_MESH_H
#define LATTICEGEN_MESH_H

#include <Core/Point3D.h>
#include <Geometry/UnitCell3D.h>
#include <Geometry/Cutoff.h>
#include <vector>

using namespace Core;

namespace Geometry
{
  class Mesh
  {
  public:
    virtual std::vector<Point3D> generate(const Cutoff &cutoff, bool shift = false) const = 0;
  };

  class LatticeMesh : public Mesh
  {
  public:
    LatticeMesh(const Cell3D &cell_) : cell(cell_)
    {}

    std::vector<Point3D> generate(const Cutoff &cutoff, bool shift = false) const final override;

    Cell3D cell;
  };

  class TrigonalMesh : public Mesh
  {
  public:
    TrigonalMesh(double a_, double c_)
      : a(a_), c(c_)
    {}

    std::vector<Point3D> generate(const Cutoff &cutoff, bool shift = false) const final override;

    double a;
    double c;
  };

  class TetrahedronMesh : public Mesh
  {
  public:
    TetrahedronMesh(double a_)
      : a(a_)
    {}

    std::vector<Point3D> generate(const Cutoff &cutoff, bool shift = false) const final override;

    double a;
  };

  class CubicMesh : public Mesh
  {
  public:
    CubicMesh(double a_)
      : a(a_)
    {}

    std::vector<Point3D> generate(const Cutoff &cutoff, bool shift = false) const final override;

    double a;
  };
}
#endif //LATTICEGEN_MESH_H
