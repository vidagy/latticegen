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
  static const Vector3D x = Vector3D(1.0,0.0,0.0);
  static const Vector3D y = Vector3D(0.0,1.0,0.0);
  static const Vector3D z = Vector3D(0.0,0.0,1.0);
}

BravaisLattice3D::UnitCell::UnitCell(const BravaisLattice3DType& type_, const Point3D& a_, const Point3D& b_, const Point3D& c_)
  : type(type_), a(a_), b(b_), c(c_)
{
  if (! a_.getLength() > 0.0)
    throw std::invalid_argument("In BravaisLattice3D::UnitCell::ctor: a must be non null vector but a = " + a_.toString());
  if (! b_.getLength() > 0.0)
    throw std::invalid_argument("In BravaisLattice3D::UnitCell::ctor: b must be non null vector but b = " + b_.toString());
  if (! c_.getLength() > 0.0)
    throw std::invalid_argument("In BravaisLattice3D::UnitCell::ctor: c must be non null vector but c = " + c_.toString());
  
  if (nearlyZero(a*b))
    throw std::invalid_argument("In BravaisLattice3D::UnitCell::ctor: a*b must not be zero a = " + a_.toString() + " b = " + b_.toString());
  if (nearlyZero(a*c))
    throw std::invalid_argument("In BravaisLattice3D::UnitCell::ctor: a*c must not be zero a = " + a_.toString() + " c = " + c_.toString());
  if (nearlyZero(b*c))
    throw std::invalid_argument("In BravaisLattice3D::UnitCell::ctor: b*c must not be zero b = " + b_.toString() + " c = " + c_.toString());
}

BravaisLattice3D::UnitCell BravaisLattice3D::UnitCell::create_triclinic(
  const double a_, const double b_, const double c_, const double alpha_, const double beta_, const double gamma_)
{
  const Vector3D v1 = x;

  const Vector3D v2 = Rotation(gamma_ * z) * x;

  const double v3x = cos(beta_);
  const double v3y = ( cos(alpha_)-v1.x*cos(beta_) )/v2.y;
  const double v3z = sqrt(1.0 - v3x*v3x - v3y*v3y);
  const Vector3D v3 = Vector3D( v3x , v3y , v3z );

  return UnitCell(BravaisLattice3D::Triclinic, a_ * v1, b_ * v2, c_ * v3);
}
BravaisLattice3D::UnitCell BravaisLattice3D::UnitCell::create_monoclinic_primitive(
  const double a_, const double b_, const double c_, const double beta_)
{
  const Vector3D v1 = a_*x;
  const Vector3D v2 = b_*y;
  const Vector3D v3 = c_* (Rotation(beta_ * y)*z);

  return UnitCell(BravaisLattice3D::Monoclinic_Primitive, v1, v2, v3);
}
BravaisLattice3D::UnitCell BravaisLattice3D::UnitCell::create_monoclinic_base(
  const double a_, const double b_, const double c_, const double beta_)
{
  const Vector3D v1 = a_*x;
  const Vector3D v2 = 1.0 / 2.0 * (a_*x + b_*y);
  const Vector3D v3 = c_*(Rotation(beta_ * y)*z);

  return UnitCell(BravaisLattice3D::Monoclinic_Base, v1, v2, v3);
}
BravaisLattice3D::UnitCell BravaisLattice3D::UnitCell::create_orthorhombic_primitive(const double a_, const double b_, const double c_)
{
  const Vector3D v1 = a_*x;
  const Vector3D v2 = b_*y;
  const Vector3D v3 = c_*z;

  return UnitCell(BravaisLattice3D::Orthorhombic_Primitive, v1, v2, v3);
}
BravaisLattice3D::UnitCell BravaisLattice3D::UnitCell::create_orthorhombic_base(const double a_, const double b_, const double c_)
{
  const Vector3D v1 = a_*x;
  const Vector3D v2 = 1.0 / 2.0 * (a_*x + b_*y);
  const Vector3D v3 = c_*z;

  return UnitCell(BravaisLattice3D::Orthorhombic_Base, v1, v2, v3);
}
BravaisLattice3D::UnitCell BravaisLattice3D::UnitCell::create_orthorhombic_body(const double a_, const double b_, const double c_)
{
  const Vector3D v1 = a_*x;
  const Vector3D v2 = b_*y;
  const Vector3D v3 = (a_*x) + (b_*y) + 1.0 / 2.0 * (c_*z);

  return UnitCell(BravaisLattice3D::Orthorhombic_Body, v1, v2, v3);
}
BravaisLattice3D::UnitCell BravaisLattice3D::UnitCell::create_orthorhombic_face(const double a_, const double b_, const double c_)
{
  const Vector3D v1 = 1.0 / 2.0 * (a_*x + b_*y);
  const Vector3D v2 = 1.0 / 2.0 * (b_*y + c_*z);
  const Vector3D v3 = 1.0 / 2.0 * (a_*x + c_*z);

  return UnitCell(BravaisLattice3D::Orthorhombic_Face, v1, v2, v3);
}
BravaisLattice3D::UnitCell BravaisLattice3D::UnitCell::create_tetragonal_primitive(const double a_, const double c_)
{
  const Vector3D v1 = a_*x;
  const Vector3D v2 = a_*y;
  const Vector3D v3 = c_*z;

  return UnitCell(BravaisLattice3D::Tetragonal_Primitive, v1, v2, v3);
}
BravaisLattice3D::UnitCell BravaisLattice3D::UnitCell::create_tetragonal_body(const double a_, const double c_)
{
  const Vector3D v1 = a_*x;
  const Vector3D v2 = a_*y;
  const Vector3D v3 = a_*(x+y) + 1.0 / 2.0 * (c_*z);

  return UnitCell(BravaisLattice3D::Tetragonal_Body, v1, v2, v3);
}
BravaisLattice3D::UnitCell BravaisLattice3D::UnitCell::create_rhombohedral(const double a_, const double alpha_)
{
  const Vector3D v1 = x;
  const Vector3D v2 = Vector3D(cos(alpha_), sin(alpha_), 0.0);
  const Vector3D v3 = Vector3D(cos(alpha_)*cos(alpha_/2.0), cos(alpha_)*sin(alpha_/2.0) , sin(alpha_) );

  return UnitCell(BravaisLattice3D::Rhombohedral, a_*v1, a_*v2, a_*v3);
}
BravaisLattice3D::UnitCell BravaisLattice3D::UnitCell::create_hexagonal(const double a_, const double c_)
{
  const Vector3D v1 = x;
  const Vector3D v2 = Vector3D(cos(pi/3.0), sin(pi/3.0), 0.0);
  const Vector3D v3 = z;

  return UnitCell(BravaisLattice3D::Hexagonal, a_*v1, a_*v2, c_*v3);
}
BravaisLattice3D::UnitCell BravaisLattice3D::UnitCell::create_cubic_primitive(const double a_)
{
  const Vector3D v1 = a_*x;
  const Vector3D v2 = a_*y;
  const Vector3D v3 = a_*z;

  return UnitCell(BravaisLattice3D::Cubic_Primitive, v1, v2, v3);
}
BravaisLattice3D::UnitCell BravaisLattice3D::UnitCell::create_cubic_body(const double a_)
{
  const Vector3D v1 = a_*x;
  const Vector3D v2 = a_*y;
  const Vector3D v3 = a_*(x + y + 1.0/2.0*z);

  return UnitCell(BravaisLattice3D::Cubic_Body, v1, v2, v3);
}
BravaisLattice3D::UnitCell BravaisLattice3D::UnitCell::create_cubic_face(const double a_)
{
  const Vector3D v1 = a_ / 2.0 * (x+y);
  const Vector3D v2 = a_ / 2.0 * (y+z);
  const Vector3D v3 = a_ / 2.0 * (y+z);

  return UnitCell(BravaisLattice3D::Cubic_Face, v1, v2, v3);
}

BravaisLattice3D::BravaisLattice3D(const UnitCell& unit_cell_, const size_t x_width_, const size_t y_width_, const size_t z_width_)
  : unit_cell(unit_cell_), x_width(x_width_), y_width(y_width_), z_width(z_width_)
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
