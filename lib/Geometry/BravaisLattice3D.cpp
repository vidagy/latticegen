#include "BravaisLattice3D.h" 

#include "Mesh.h"

using namespace Geometry;
using namespace Core;

BravaisLattice3D::BravaisLattice3D(const UnitCell3D& unit_cell_, const size_t x_width_, const size_t y_width_, const size_t z_width_)
  : unit_cell(unit_cell_)
  , x_width(x_width_)
  , y_width(y_width_)
  , z_width(z_width_)
  , lattice(LatticeMesh(unit_cell, true).generate(CutoffUnitVectors(unit_cell, x_width, y_width, z_width)))
{
}
