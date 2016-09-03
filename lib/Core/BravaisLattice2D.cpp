#include "BravaisLattice2D.h" 
#include "ComparisonHelpers.h"

#include <math.h>
#include <stdexcept>
#include <algorithm>
#include <cassert>

namespace Core
{
  namespace Geometry
  {

    std::pair<Point2D,double> BravaisLattice2D::get_canonical_unit_cell_and_scale(Point2D x, Point2D y)
    {
      // validate input
      if ( nearlyZero(x.getLength()) || nearlyZero(y.getLength()))
      {
        throw std::invalid_argument("Null vector on the input of get_canonical_unit_cell");
      }
      if ( equalsWithTolerance( x.x_ * y.y_, x.y_ * y.x_) )
      {
        throw std::invalid_argument("Non-independent vectors on the input of get_canonical_unit_cell");
      }

      // x will be the longer
      const double length_x = x.getLength();
      const double length_y = y.getLength();
      if (length_y > length_x)
        swap(x,y);

      const double length = std::max(length_x,length_y);
      // rotation and scaling to (1.0, 0.0)
      const double sin_theta = x.y_ / length / length;
      const double cos_theta = x.x_ / length / length;

      Point2D a = Point2D(
        cos_theta * y.x_ + sin_theta * y.y_,
        cos_theta * y.y_ - sin_theta * y.x_);
      
      // invert and shift to the (++) quarter plane
      if (a.y_ < 0.0)
        a = a * -1.0;
      while (a.x_ < 0.0)
      {
        a.x_ += 1.0;
      }
      return std::make_pair(a,length);
    }

    BravaisLattice2D::BravaisLattice2DType BravaisLattice2D::find_lattice_type(const Point2D& a)
    {
      double length_a = a.getLength();

      assert(length_a > 0.0);
      assert(a.x_ >= 0.0 || a.y_ >= 0.0);

      if ( nearlyZero(length_a - 1.0) )
      {
        if ( nearlyZero(a.x_) )
        {
          return Square;
        }
        else if ( equalsWithTolerance( a.x_, 0.5) )
        {
          return Hexagonal;
        }
      }
      else
      {
        if ( nearlyZero(a.x_) )
        {
          return Rectangular;
        }
        else if ( equalsWithTolerance( a.x_, 0.5) )
        {
          return CenteredRectangular;
        }
      }
      return Oblique;
    }

    BravaisLattice2D::BravaisLattice2D(const Point2D& unit_vector, const double scale, const size_t width_x, const size_t width_y)
    : unit_vector_(unit_vector), scale_(scale), width_x_(width_x), width_y_(width_y)
    {
      lattice_.reserve(width_x_ * width_y_);
      Point2D x0 = Point2D(1.0,0.0);
      
      for (unsigned int j = 0; j < width_y_; ++j)
      {
        for (unsigned int i = 0; i < width_x_; ++i)
        {
          lattice_.push_back( i*x0 + j*unit_vector_ );
        }
      }
    }
  }
}
