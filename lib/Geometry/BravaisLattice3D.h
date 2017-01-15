#ifndef LATTICEGEN_BRAVAISLATTICE3D_H_
#define LATTICEGEN_BRAVAISLATTICE3D_H_

#include "Point3D.h"
#include "UnitCell3D.h"
#include <vector>

namespace Geometry
{
  class BravaisLattice3D
  {
  public:
    typedef std::vector<Point3D> Point3DVec;

    BravaisLattice3D(const UnitCell3D& unit_cell_, const size_t x_width, const size_t y_width, const size_t z_width);

    UnitCell3D get_unit_cell() const { return unit_cell; };
    size_t     get_x_width()   const { return x_width; }
    size_t     get_y_width()   const { return y_width; }
    size_t     get_z_width()   const { return z_width; }
    Point3DVec get_lattice()   const { return lattice; }

  private:
    UnitCell3D unit_cell;
    size_t x_width;
    size_t y_width;
    size_t z_width;
    Point3DVec lattice;
  };
}

#endif // LATTICEGEN_BRAVAISLATTICE3D_H_