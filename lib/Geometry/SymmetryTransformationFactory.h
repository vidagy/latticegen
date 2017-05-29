#ifndef LATTICEGEN_SYMMETRYTARNSFORMATIONFACTORY_H
#define LATTICEGEN_SYMMETRYTARNSFORMATIONFACTORY_H

#include <Geometry/UnitCell3D.h>
#include <Geometry/Matrix3D.h>
#include <Geometry/CrystallographicPointGroups.h>
#include <Geometry/Transformations.h>

namespace Geometry
{
  class SymmetryTransformationFactory
  {
  public:
    typedef std::vector<Matrix3D> Transformations;

    static Transformation get(CrystallographicPointGroup::SymmetryElement symmetry_element);
    static Transformations get(const CrystallographicPointGroup::Elements& elements);

    static Transformations generate(const CrystallographicPointGroup::Elements& generators);

    static Transformations generate(const Cell3D &cell);
  };
}

#endif //LATTICEGEN_SYMMETRYTARNSFORMATIONFACTORY_H
