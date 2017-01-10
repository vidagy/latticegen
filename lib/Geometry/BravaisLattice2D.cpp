#include "BravaisLattice2D.h" 

#include <algorithm>

using namespace Core;
using namespace Geometry;

namespace
{
  static const double pi = 3.14159265358979323846;
}

BravaisLattice2D::UnitCell::UnitCell(const Point2D& a_, const Point2D& b_)
  : a(a_), b(b_)
{
  if (a_.length() <= 0.0)
    throw std::invalid_argument("In BravaisLattice2D::UnitCell::ctor: a must be non null vector but a = " + std::to_string(a_));
  if (b_.length() <= 0.0)
    throw std::invalid_argument("In BravaisLattice2D::UnitCell::ctor: b must be non null vector but b = " + std::to_string(b_));
  
  if (a_.x <= 0.0)
    throw std::invalid_argument("In BravaisLattice2D::UnitCell::ctor: a.x must be positive but a = " + std::to_string(a_));
  if (! nearlyZero(a_.y) )
    throw std::invalid_argument("In BravaisLattice2D::UnitCell::ctor: a.y must be zero but a = " + std::to_string(a_));

  if (! positiveWithTolerance(b_.x))
    throw std::invalid_argument("In BravaisLattice2D::UnitCell::ctor: b.x must be positive with tolerance but b = " + std::to_string(b_));
  if (b_.y <= 0.0)
    throw std::invalid_argument("In BravaisLattice2D::UnitCell::ctor: b.y must be positive but b = " + std::to_string(b_));
  
  if (! greaterEqualsWithTolerance(a_.length(), b_.length()))
    throw std::invalid_argument(
      "In BravaisLattice2D::UnitCell::ctor: a must be longer than b but a = " + std::to_string(a_) + " and b = " + std::to_string(b)
      );
}

BravaisLattice2D::UnitCell BravaisLattice2D::get_unit_cell(Point2D x, Point2D y)
{
  // validate input
  if ( nearlyZero(x.length()) || nearlyZero(y.length()))
  {
    throw std::invalid_argument(
      "Null vector on the input of BravaisLattice2D::get_unit_cell: x = " + std::to_string(x) + " y = " + std::to_string(y)
      );
  }
  if ( nearlyZero( x.x * y.y - x.y * y.x) )
  {
    throw std::invalid_argument(
      "Non-independent vectors on the input of BravaisLattice2D::get_unit_cell: x = " + std::to_string(x) + " y = " + std::to_string(y)
      );
  }

  // x will be the longer
  const double length_x = x.length();
  const double length_y = y.length();
  if (length_y > length_x)
    swap(x,y);

  const double length = std::max(length_x,length_y);
  // rotation
  const double sin_theta = x.y / length;
  const double cos_theta = x.x / length;

  Point2D b = Point2D(
    cos_theta * y.x + sin_theta * y.y,
    cos_theta * y.y - sin_theta * y.x);
  
  // invert and shift to the (++) quarter plane
  if (b.y < 0.0)
    b = b * -1.0;
  while (b.x < 0.0)
  {
    b.x += length;
  }
  return UnitCell( Point2D(length, 0.0), b );
}

BravaisLattice2D::BravaisLattice2DType BravaisLattice2D::find_lattice_type(const UnitCell& unit_cell_)
{
  const Point2D& a = unit_cell_.a;
  const Point2D& b = unit_cell_.b;

  double length_a = a.length();
  double length_b = b.length();

  if ( nearlyZero(length_b - length_a) )
  {
    if ( nearlyZero(b.x) )
    {
      return Square;
    }
    else if ( equalsWithTolerance( b.x, 0.5 * a.x) )
    {
      return Hexagonal;
    }
  }
  else
  {
    if ( nearlyZero(b.x) )
    {
      return Rectangular;
    }
    else if ( equalsWithTolerance( b.x, 0.5 * a.x) )
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
    // x is in [0.0, a.x], y > 0.0 and < line (0,b) and y < line(a,b)

    BravaisLattice2D::Point2DVec result;
    result.reserve( xsample*ysample / 2 + 1);

    const double xstep = a.x / xsample;
    const double ystep = b.y / ysample;

    const double xshift = xstep / 2.0;
    const double yshift = ystep / 2.0;

    for (unsigned int i = 0; i < xsample; ++i)
    {
      const double x = xshift + i * xstep;
      for (unsigned int j = 0; j < ysample; ++j)
      {
        const double y = yshift + j * xstep;
        if ( 
          x < b.x 
          ? x*b.y/b.x < y 
          : (a.x-x)*b.y/(a.x-b.x) < y )
          break;
        result.push_back(Point2D(x,y));
      }
    }

    return result;
  }
}

BravaisLattice2D::Point2DVec BravaisLattice2D::get_irreducible_wedge(const UnitCell& unit_cell_, const unsigned int xsample, const unsigned int ysample)
{
  // we will sample the convex combination of (0,0), x and y, where x and y is determined by the lattice type and unit_cell_.a.
  BravaisLattice2DType type = find_lattice_type(unit_cell_);

  switch (type)
  {
    case Oblique:
    {
      return sample_between(unit_cell_.a, unit_cell_.b, xsample, ysample);
    }
    case Rectangular:
    {
      return sample_between(unit_cell_.a, unit_cell_.b, xsample, ysample);
    }
    case CenteredRectangular:
    {
      Point2D x = unit_cell_.a * 0.5;
      Point2D y = unit_cell_.b;
      return sample_between(x, y, xsample, ysample);
    }
    case Hexagonal:
    {
      Point2D x = unit_cell_.a;
      Point2D y = (unit_cell_.a + unit_cell_.b) * 0.5;
      return sample_between(x, y, xsample, ysample);
    }
    case Square:
    {
      Point2D x = unit_cell_.a * 0.5;
      Point2D y = (unit_cell_.a + unit_cell_.b) * 0.5;
      return sample_between(x, y, xsample, ysample);
    }
    default:
      throw std::logic_error("Unhandled BravaisLattice2DType.");
  }
}

std::pair<Point2D,Point2D> BravaisLattice2D::get_inverse_unit_cell(const UnitCell& unit_cell_)
{
  const Point2D& a = unit_cell_.a;
  const Point2D& b = unit_cell_.b;

  Point2D k1 = Point2D(2.0*pi, -2.0*pi*b.x/b.y/a.x);
  Point2D k2 = Point2D(0.0,2.0*pi/b.y/a.x);

  return std::make_pair(k1,k2);
}

BravaisLattice2D::BravaisLattice2D(const Point2D& a, const Point2D& b, const size_t x_width_, const size_t y_width_)
: unit_cell(get_unit_cell(a, b)), x_width(x_width_), y_width(y_width_)
{
  lattice.reserve(x_width * y_width);
  
  for (unsigned int j = 0; j < y_width; ++j)
  {
    for (unsigned int i = 0; i < x_width; ++i)
    {
      lattice.push_back( i*unit_cell.a + j*unit_cell.b );
    }
  }
}
