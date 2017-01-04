#include "BravaisLattice3D.h" 
#include "ComparisonHelpers.h"
#include "SymmetryOperations.h"

#include <cmath>
#include <stdexcept>
#include <string>

using namespace Core::Geometry;
using namespace Core;

namespace
{
  static const double pi = 3.14159265358979323846;
}

BravaisLattice3D::UnitCell::UnitCell(const Point3D& a_, const Point3D& b_, const Point3D& c_)
  : a(a_), b(b_), c(c_)
{
  if (! a_.getLength() > 0.0)
    throw std::invalid_argument("In BravaisLattice3D::UnitCell::ctor: a must be non null vector but a = " + a_.toString());
  if (! b_.getLength() > 0.0)
    throw std::invalid_argument("In BravaisLattice3D::UnitCell::ctor: b must be non null vector but b = " + b_.toString());
  if (! c_.getLength() > 0.0)
    throw std::invalid_argument("In BravaisLattice3D::UnitCell::ctor: c must be non null vector but c = " + c_.toString());
  
  if (! a_.x > 0.0)
    throw std::invalid_argument("In BravaisLattice3D::UnitCell::ctor: a.x must be positive but a = " + a_.toString());
  if (! nearlyZero(a_.y) )
    throw std::invalid_argument("In BravaisLattice3D::UnitCell::ctor: a.y must be zero but a = " + a_.toString());
  if (! nearlyZero(a_.z) )
    throw std::invalid_argument("In BravaisLattice3D::UnitCell::ctor: a.z must be zero but a = " + a_.toString());
  
  if (! positiveWithTolerance(b_.x))
    throw std::invalid_argument("In BravaisLattice3D::UnitCell::ctor: b.x must be positive with tolerance but b = " + b_.toString());
  if (! b_.y > 0.0)
    throw std::invalid_argument("In BravaisLattice3D::UnitCell::ctor: b.y must be positive but b = " + b_.toString());
  if (! nearlyZero(b_.z) )
    throw std::invalid_argument("In BravaisLattice3D::UnitCell::ctor: b.z must be zero but a = " + b_.toString());
  
  if (! positiveWithTolerance(c_.x))
    throw std::invalid_argument("In BravaisLattice3D::UnitCell::ctor: c.x must be positive with tolerance but c = " + c_.toString());
  if (! positiveWithTolerance(c_.y))
    throw std::invalid_argument("In BravaisLattice3D::UnitCell::ctor: c.y must be positive with tolerance but c = " + c_.toString());
  if (! c_.z > 0.0 )
    throw std::invalid_argument("In BravaisLattice3D::UnitCell::ctor: c.z must be positive but c = " + c_.toString());
  
  if (! greaterEqualsWithTolerance(a_.getLength(), b_.getLength()))
    throw std::invalid_argument(
      "In BravaisLattice3D::UnitCell::ctor: a must be longer than b but a = " + a_.toString() + " and b = " + b.toString()
      );
  if (! greaterEqualsWithTolerance(b_.getLength(), c_.getLength()))
    throw std::invalid_argument(
      "In BravaisLattice3D::UnitCell::ctor: b must be longer than c but b = " + a_.toString() + " and c = " + b.toString()
      );
}

namespace
{
  inline double get_volume(const Point3D& x, const Point3D& y, const Point3D& z)
  {
    return fabs(x * cross_product(y,z));
  }

  std::tuple<Point3D,Point3D,Point3D> to_same_octad(Point3D x, Point3D y, Point3D z)
  {
    const double xy = x*y;
    const double yz = y*z;
    const double xz = x*z;

    const double volume = get_volume(x,y,z);

    // all acute
    if (positiveWithTolerance(xy) && positiveWithTolerance(yz) && positiveWithTolerance(xz))
      return std::make_tuple(x,y,z);

    // all obtuse
    if ( xy < 0.0 && yz < 0.0 && xz < 0.0 )
    {
      Vector3D temp_x = (-1.0 * x) + y + z;
      Vector3D temp_y = x + (-1.0 * y) + z;
      Vector3D temp_z = x + y + (-1.0 * z);
      if (equalsWithTolerance(volume, get_volume(temp_x, temp_y, temp_z)))
        return std::make_tuple(x,y,z);    

      throw std::logic_error("could not get all vectors to have acute angles: x = " + x.toString()
        + " y = " + y.toString()
        + " z = " + z.toString()
        );      
    }
    
    // one is obtuse angle
    if ( negativeWithTolerance(xy * yz * xz) )
    {
      if ( negativeWithTolerance(yz) )
        swap(x,z);
      if ( negativeWithTolerance(xz) )
        swap(y,z);

      Vector3D temp = x+y;
      if ( positiveWithTolerance(temp*y) && positiveWithTolerance(temp*z) && equalsWithTolerance(volume, get_volume(temp,y,z)) )
        return std::make_tuple(temp,y,z);
      if ( positiveWithTolerance(temp*x) && positiveWithTolerance(temp*z) && equalsWithTolerance(volume, get_volume(x,temp,z)) )
        return std::make_tuple(x,temp,z);

      throw std::logic_error("could not get all vectors to have acute angles: x = " + x.toString()
        + " y = " + y.toString()
        + " z = " + z.toString()
        );
    }

    // two are obtuse angles
    if ( positiveWithTolerance(xy * yz * xz) )
    {
      if ( positiveWithTolerance(yz) )
        swap(x,z);
      if ( positiveWithTolerance(xz) )
        swap(y,z);

      Vector3D temp = x+z;
      if ( positiveWithTolerance(temp*x) && positiveWithTolerance(temp*y) && equalsWithTolerance(volume, get_volume(x,y,temp)) )
        return std::make_tuple(x,y,temp);
      temp = y+z;
      if ( positiveWithTolerance(temp*x) && positiveWithTolerance(temp*y) && equalsWithTolerance(volume, get_volume(x,y,temp)) )
        return std::make_tuple(x,y,temp);
      temp = x+y+z;
      if ( positiveWithTolerance(temp*x) && positiveWithTolerance(temp*y) && equalsWithTolerance(volume, get_volume(x,y,temp)) )
        return std::make_tuple(x,y,temp);

      throw std::logic_error("could not get all vectors to have acute angles: x = " + x.toString()
        + " y = " + y.toString()
        + " z = " + z.toString()
        );
    }

    throw std::logic_error("unhandled angles in to_same_octad:  x = " + x.toString()
      + " y = " + y.toString()
      + " z = " + z.toString()
      );
  }

  std::tuple<Point3D,Point3D,Point3D> order_vectors(Point3D x, Point3D y, Point3D z)
  {
    if (y.getLength() > x.getLength())
      swap(x,y);
    if (z.getLength() > x.getLength())
      swap(x,z);
    if (z.getLength() > y.getLength())
      swap(y,z);

    return std::make_tuple(x,y,z);
  }

  std::tuple<Point3D,Point3D,Point3D> get_primitive_cell(Point3D a, Point3D b, Point3D c)
  {
    // they are on the same octad
    // so only rotations are necessary.

    // rotate so a becomes parallel to x
    const double length_a = a.getLength();
    const double angle_x = acos(a.x / length_a);
    const Vector3D rotation_vec_x = -1.0 * angle_x / length_a * cross_product(Vector3D(1.0,0.0,0.0), a);
    const Rotation rotation_x = Rotation(rotation_vec_x);

    a = rotation_x(a);
    b = rotation_x(b);
    c = rotation_x(c);

    if ( (!nearlyZero(a.y)) || (!nearlyZero(a.z)) )
      throw std::logic_error("a must be parallel to x after rotation but a = " + a.toString());

    // rotate so b.z = 0
    const double length_b = b.getLength();
    const double angle_y = acos(b.y / length_b);
    const Vector3D rotation_vec_y = -1.0 * angle_y * Vector3D(1.0,0.0,0.0);
    const Rotation rotation_y = Rotation(rotation_vec_y);

    a = rotation_y(a);
    b = rotation_y(b);
    c = rotation_y(c);

    if ( (!nearlyZero(a.y)) || (!nearlyZero(a.z)) )
      throw std::logic_error("a must be parallel to x after second rotation but a = " + a.toString());
    if ( !nearlyZero(b.z) )
      throw std::logic_error("b must be perpendicular to z after second rotation but b = " + b.toString());

    return std::make_tuple(a,b,c);
  }
}

BravaisLattice3D::UnitCell BravaisLattice3D::get_unit_cell(Point3D x, Point3D y, Point3D z)
{
  // validate input
  if ( nearlyZero(x.getLength()) || nearlyZero(y.getLength()) || nearlyZero(z.getLength()))
  {
    throw std::invalid_argument(
      "Null vector on the input of BravaisLattice3D::get_unit_cell: x = " + 
      x.toString() + " y = " + y.toString() + " z = " + z.toString()
      );
  }
  const double determinant = 
      x.x * (y.y * z.z - z.y * y.z)
        - y.x * (x.y * z.z - z.y * x.z)
        + z.x * (x.y * y.z - y.y * x.z);

  if (nearlyZero(determinant))
  {
    throw std::invalid_argument(
      "Non-independent vectors on the input of BravaisLattice3D::get_unit_cell: x = " + 
      x.toString() + " y = " + y.toString() + " z = " + z.toString()
      );
  }

  std::tie(x, y, z) = to_same_octad(x, y, z);

  // ordering is length(x) >= length(y) >= length(z) 
  std::tie(x, y, z) = order_vectors(x, y, z);

  Point3D a, b, c;
  std::tie(a, b, c) = get_primitive_cell(x, y, z);

  return UnitCell(a, b, c);
}

/*
BravaisLattice3D::BravaisLattice3DType BravaisLattice3D::find_lattice_type(const UnitCell& unit_cell_)
{
  const Point3D& a = unit_cell_.a;
  const Point3D& b = unit_cell_.b;

  double length_a = a.getLength();
  double length_b = b.getLength();

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
  BravaisLattice3D::Point3DVec sample_between(const Point3D& a, const Point3D& b, const unsigned int xsample, const unsigned int ysample)
  {
    // x is in [0.0, a.x], y > 0.0 and < line (0,b) and y < line(a,b)

    BravaisLattice3D::Point3DVec result;
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
        result.push_back(Point3D(x,y));
      }
    }

    return result;
  }
}

BravaisLattice3D::Point3DVec BravaisLattice3D::get_irreducible_wedge(const UnitCell& unit_cell_, const unsigned int xsample, const unsigned int ysample)
{
  // we will sample the convex combination of (0,0), x and y, where x and y is determined by the lattice type and unit_cell_.a.
  BravaisLattice3DType type = find_lattice_type(unit_cell_);

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
      Point3D x = unit_cell_.a * 0.5;
      Point3D y = unit_cell_.b;
      return sample_between(x, y, xsample, ysample);
    }
    case Hexagonal:
    {
      Point3D x = unit_cell_.a;
      Point3D y = (unit_cell_.a + unit_cell_.b) * 0.5;
      return sample_between(x, y, xsample, ysample);
    }
    case Square:
    {
      Point3D x = unit_cell_.a * 0.5;
      Point3D y = (unit_cell_.a + unit_cell_.b) * 0.5;
      return sample_between(x, y, xsample, ysample);
    }
    default:
      throw std::logic_error("Unhandled BravaisLattice3DType.");
  }
}

std::pair<Point3D,Point3D> BravaisLattice3D::get_inverse_unit_cell(const UnitCell& unit_cell_)
{
  const Point3D& a = unit_cell_.a;
  const Point3D& b = unit_cell_.b;

  Point3D k1 = Point3D(2.0*pi, -2.0*pi*b.x/b.y/a.x);
  Point3D k2 = Point3D(0.0,2.0*pi/b.y/a.x);

  return std::make_pair(k1,k2);
}
*/

BravaisLattice3D::BravaisLattice3D( 
  const Point3D& a_, const Point3D& b_, const Point3D& c_, 
  const size_t x_width_, const size_t y_width_, const size_t z_width_
  )
  : unit_cell(get_unit_cell(a_, b_, c_)), x_width(x_width_), y_width(y_width_), z_width(z_width_)
{
  lattice.reserve(x_width * y_width * z_width);
  
  for (unsigned int k = 0; k < z_width; ++k)
  {
    auto z = k*unit_cell.c;
    for (unsigned int j = 0; j < y_width; ++j)
    {
      auto y = j*unit_cell.b;
      for (unsigned int i = 0; i < x_width; ++i)
      {
        lattice.push_back( i*unit_cell.a + y + z );
      }
    }
  }
}
