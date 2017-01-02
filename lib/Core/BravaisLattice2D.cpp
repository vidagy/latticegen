#include "BravaisLattice2D.h" 
#include "ComparisonHelpers.h"

#include <math.h>
#include <stdexcept>
#include <algorithm>
#include <cassert>

using namespace Core::Geometry;

namespace
{
  static const double pi = 3.14159265358979323846;
}

std::pair< Point2D, Point2D > BravaisLattice2D::get_canonical_unit_cell(Point2D x, Point2D y)
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
  // rotation
  const double sin_theta = x.y_ / length;
  const double cos_theta = x.x_ / length;

  Point2D b = Point2D(
    cos_theta * y.x_ + sin_theta * y.y_,
    cos_theta * y.y_ - sin_theta * y.x_);
  
  // invert and shift to the (++) quarter plane
  if (b.y_ < 0.0)
    b = b * -1.0;
  while (b.x_ < 0.0)
  {
    b.x_ += length;
  }
  return std::make_pair( Point2D(length, 0.0), b );
}

BravaisLattice2D::BravaisLattice2DType BravaisLattice2D::find_lattice_type(const Point2D& a, const Point2D& b)
{
  double length_a = a.getLength();
  double length_b = b.getLength();

  assert(length_a > 0.0);
  assert(length_b > 0.0);
  assert( nearlyZero(a.y_) );
  assert(b.x_ >= 0.0 || b.y_ >= 0.0);

  if ( nearlyZero(length_b - length_a) )
  {
    if ( nearlyZero(b.x_) )
    {
      return Square;
    }
    else if ( equalsWithTolerance( b.x_, 0.5 * a.x_) )
    {
      return Hexagonal;
    }
  }
  else
  {
    if ( nearlyZero(b.x_) )
    {
      return Rectangular;
    }
    else if ( equalsWithTolerance( b.x_, 0.5 * a.x_) )
    {
      return CenteredRectangular;
    }
  }
  return Oblique;
}

namespace
{
  BravaisLattice2D::Point2DVec sample_between(const Point2D& a, const Point2D& b, const unsigned int xsample, const unsigned int ysample)
  {
    assert( Core::nearlyZero(a.y_) );
    assert( b.x_ > 0.0 );
    assert( b.y_ > 0.0 );

    assert( b.x_ < a.x_ );
    
    // x is in [0.0, a.x_], y > 0.0 and < line (0,b) and y < line (a,b)

    BravaisLattice2D::Point2DVec result;
    result.reserve( xsample*ysample / 2 + 1);

    const double xstep = a.x_ / xsample;
    const double ystep = b.y_ / ysample;

    const double xshift = xstep / 2.0;
    const double yshift = ystep / 2.0;

    for (unsigned int i = 0; i < xsample; ++i)
    {
      const double x = xshift + i * xstep;
      for (unsigned int j = 0; j < ysample; ++j)
      {
        const double y = yshift + j * xstep;
        if ( 
          x < b.x_ 
          ? x*b.y_/b.x_ < y 
          : (a.x_-x)*b.y_/(a.x_-b.x_) < y )
          break;
        result.push_back(Point2D(x,y));
      }
    }

    return result;
  }
}

BravaisLattice2D::Point2DVec BravaisLattice2D::get_irreducible_wedge(const Point2D& a, const Point2D& b, const unsigned int xsample, const unsigned int ysample)
{
  // we will sample the convex combination of (0,0), x and y, where x and y is determined by the lattice type and a.

  assert( nearlyZero(a.y_) );
  assert( b.x_ > 0.0 );
  assert( b.y_ > 0.0 );

  assert( b.x_ < a.x_ );
  
  BravaisLattice2DType type = find_lattice_type(a,b);

  switch (type)
  {
    case Oblique:
    {
      return sample_between(a,b,xsample,ysample);
    }
    case Rectangular:
    {
      return sample_between(a,b,xsample,ysample);
    }
    case CenteredRectangular:
    {
      Point2D x = a * 0.5;
      Point2D y = b;
      return sample_between(x,y,xsample,ysample);
    }
    case Hexagonal:
    {
      Point2D x = a;
      Point2D y = (a + b) * 0.5;
      return sample_between(x,y,xsample,ysample);
    }
    case Square:
    {
      Point2D x = a * 0.5;
      Point2D y = (a + b) * 0.5;
      return sample_between(x,y,xsample,ysample);
    }
    default:
      throw std::logic_error("Unhandled BravaisLattice2DType.");
  }
}

std::pair<Point2D,Point2D> BravaisLattice2D::get_inverse_unit_cell(const Point2D& a, const Point2D& b)
{
  assert( nearlyZero(a.y_) );
  assert( ! nearlyZero(a.x_) );
  assert( ! nearlyZero(b.y_) );

  const Point2D k1 = Point2D(2.0*pi, -2.0*pi*b.x_/b.y_/a.x_);
  const Point2D k2 = Point2D(0.0,2.0*pi/b.y_/a.x_);

  return std::make_pair(k1,k2);
}

BravaisLattice2D::BravaisLattice2D(const Point2D& a, const Point2D& b, const size_t width_x, const size_t width_y)
: a_(a), b_(b), width_x_(width_x), width_y_(width_y)
{
  lattice_.reserve(width_x_ * width_y_);
  
  for (unsigned int j = 0; j < width_y_; ++j)
  {
    for (unsigned int i = 0; i < width_x_; ++i)
    {
      lattice_.push_back( i*a_ + j*b_ );
    }
  }
}
