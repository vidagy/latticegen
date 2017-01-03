#ifndef LATTICEGEN_BRAVAISLATTICE2D_H_
#define LATTICEGEN_BRAVAISLATTICE2D_H_

#include "Point2D.h"
#include <utility>
#include <vector>

namespace Core
{
  namespace Geometry
  {
    class BravaisLattice2D
    {
    public:
      typedef std::vector<Point2D> Point2DVec;

      enum BravaisLattice2DType 
      {
        Oblique,
        Rectangular,
        CenteredRectangular,
        Hexagonal,
        Square
      };
      
      // canonical unit cell: first cell vector is (a, 0.0), second is (x,y) where x > 0.0 and y >= 0.0 and (x^2 + y^2) < a^2
      struct UnitCell 
      {
        Point2D a; 
        Point2D b;

      private: 
        UnitCell(const Point2D& a_, const Point2D& b_);

        friend BravaisLattice2D;
      };

      static UnitCell get_unit_cell(Point2D x, Point2D y);
      static BravaisLattice2DType find_lattice_type(const UnitCell& unit_cell);
      static Point2DVec get_irreducible_wedge(const UnitCell& unit_cell, const unsigned int xsample, const unsigned int ysample);
      
      static std::pair<Point2D,Point2D> get_inverse_unit_cell(const UnitCell& unitCell);
      
      BravaisLattice2D(const Point2D& a, const Point2D& b, const size_t x_width, const size_t y_width);

      UnitCell   get_unit_cell() const { return unit_cell; };
      size_t     get_x_width()   const { return x_width; }
      size_t     get_y_width()   const { return y_width; }
      Point2DVec get_lattice()   const { return lattice; }

    private:
      UnitCell unit_cell;
      size_t x_width;
      size_t y_width;
      Point2DVec lattice;
    };
  }
}
#endif // LATTICEGEN_BRAVAISLATTICE_H_
