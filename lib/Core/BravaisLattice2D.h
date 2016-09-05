#ifndef LATTICEGEN_BRAVAISLATTICE_H_
#define LATTICEGEN_BRAVAISLATTICE_H_

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

      enum BravaisLattice2DType {
        Oblique,
        Rectangular,
        CenteredRectangular,
        Hexagonal,
        Square
      };

      // canonical unit cell: first cell vector is (a, 0.0), second is (x,y) where x > 0.0 and y >= 0.0 and (x^2 + y^2) < a^2
      static std::pair<Point2D,Point2D> get_canonical_unit_cell(Point2D x, Point2D y);
      static BravaisLattice2DType find_lattice_type(const Point2D& a, const Point2D& b );
      static Point2DVec get_irreducible_wedge(const Point2D& a, const Point2D& b, const unsigned int xsample, const unsigned int ysample);
      
      static std::pair<Point2D,Point2D> get_inverse_unit_cell(const Point2D& a, const Point2D& b);
      
      BravaisLattice2D(const Point2D& a, const Point2D& b, const size_t width_x, const size_t width_y);

      Point2D    get_a()       const { return a_; };
      Point2D    get_b()       const { return b_; };
      size_t     get_width_x() const { return width_x_; }
      size_t     get_width_y() const { return width_y_; }
      Point2DVec get_lattice() const { return lattice_; }

    private:
      Point2D a_;
      Point2D b_;
      double scale_;
      size_t width_x_;
      size_t width_y_;
      Point2DVec lattice_;
    };
  }
}
#endif // LATTICEGEN_BRAVAISLATTICE_H_
