#ifndef LATTICEGEN_IRREDUCIBLEWEDGE_H
#define LATTICEGEN_IRREDUCIBLEWEDGE_H

#include <vector>
#include "Point3D.h"
#include "SymmetryTransformationFactory.h"
#include "UnitCell3D.h"

namespace Geometry
{
  class IrreducibleWedge
  {
  public:
    static std::vector<Point3D> reduce_by_symmetries(
      const std::vector<Point3D>& points,
      const SymmetryTransformationFactory::Transformations& transformations);

    static std::vector<Point3D> get_irreducible_wedge(const UnitCell3D& unit_cell, const size_t sample);
  };
}

#endif //LATTICEGEN_IRREDUCIBLEWEDGE_H
