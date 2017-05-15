#ifndef LATTICEGEN_STRUCTURECONSTANTS_H
#define LATTICEGEN_STRUCTURECONSTANTS_H

#include <complex>
#include <Geometry/UnitCell3D.h>

using namespace Geometry;

namespace Physics
{
  namespace NonRelativistic
  {
    class RealSpaceStructureConstants
    {
      RealSpaceStructureConstants(const UnitCell3D &unit_cell_) : unit_cell(unit_cell_) {}

      std::complex<double> calculate(
        unsigned int l, unsigned int m, unsigned int lprime, unsigned int mprime,
        Coordinates3D n, Coordinates3D nprime, const std::complex<double> &z
      ) const;

    private:
      const UnitCell3D unit_cell;
    };
  }
}

#endif //LATTICEGEN_STRUCTURECONSTANTS_H
