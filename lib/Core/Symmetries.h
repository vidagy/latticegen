#ifndef LATTICEGEN_SYMMETRIES_H_
#define LATTICEGEN_SYMMETRIES_H_

#include "Point2D.h"
#include <utility>

namespace Core
{
  namespace Geometry
  {
    class BravaisLattice2D
    {
    public:
      enum BravaisLattice2DType {
        Oblique,
        Rectangular,
        CenteredRectangular,
        Hexagonal,
        Square
      };

      // canonical unit cell: first cell vector is (1.0, 0.0), second is (x,y) where x > 0.0 and y >= 0.0 and (x^2 + y^2) < 1.0
      // since the first unit vector is trivial, we omit it.
      static std::pair<Point2D,double> get_canonical_unit_cell_and_scale(Point2D x, Point2D y);

      static BravaisLattice2DType find_lattice_type(const Point2D& a );

    };
  }
}
#endif // LATTICEGEN_SYMMETRIES_H_
