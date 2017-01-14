#include "UnitCell3D.h" 

#include "Transformations.h"

#include <cmath>

using namespace Geometry;
using namespace Core;

#define CHECK_LENGTH_POSITIVE(x) \
    if (! strictlyPositive(x))         \
      { throw std::invalid_argument(std::string("Lattice vector length ") + #x  + " must be positive"); } ; 

#define CHECK_LENGTH_DIFFERENT(x, y)   \
    if (equalsWithTolerance(x, y))     \
      { throw std::invalid_argument(std::string("Lattice vector length ") + #x  + " and " + #y + " must differ"); } ; 

#define CHECK_ANGLE_ACUTE(x)            \
    if (! strictlyPositive(x) || ! strictlyLess(x, pi/2.0) )    \
      throw std::invalid_argument(std::string("Lattice angle ") + #x  + " must be acute");

#define CHECK_ANGLE_DIFFERENT(x, y)   \
    if (equalsWithTolerance(x, y))     \
      { throw std::invalid_argument(std::string("Lattice angle ") + #x  + " and " + #y + " must differ"); } ; 

#define CHECK_LENGTH2(x, y) \
    CHECK_LENGTH_POSITIVE(x); \
    CHECK_LENGTH_POSITIVE(y); \
    CHECK_LENGTH_DIFFERENT(x,y); 

#define CHECK_LENGTH(x, y, z) \
    CHECK_LENGTH_POSITIVE(x); \
    CHECK_LENGTH_POSITIVE(y); \
    CHECK_LENGTH_POSITIVE(z); \
    CHECK_LENGTH_DIFFERENT(x,y); \
    CHECK_LENGTH_DIFFERENT(x,z); \
    CHECK_LENGTH_DIFFERENT(y,z);

#define CHECK_ANGLE(x, y, z) \
    CHECK_ANGLE_ACUTE(x); \
    CHECK_ANGLE_ACUTE(y); \
    CHECK_ANGLE_ACUTE(z); \
    CHECK_ANGLE_DIFFERENT(x,y); \
    CHECK_ANGLE_DIFFERENT(x,z); \
    CHECK_ANGLE_DIFFERENT(y,z);

namespace
{
  static const double pi = 3.14159265358979323846;
  static const Vector3D x = Vector3D(1.0,0.0,0.0);
  static const Vector3D y = Vector3D(0.0,1.0,0.0);
  static const Vector3D z = Vector3D(0.0,0.0,1.0);
}

UnitCell3D::UnitCell3D(const BravaisLattice3DType& type_, const Point3D& a_, const Point3D& b_, const Point3D& c_)
  : type(type_), a(a_), b(b_), c(c_)
{
  if (a_.length() <= 0.0)
    throw std::invalid_argument("In UnitCell3D::ctor: a must be non null vector but a = " + std::to_string(a_));
  if (b_.length() <= 0.0)
    throw std::invalid_argument("In UnitCell3D::ctor: b must be non null vector but b = " + std::to_string(b_));
  if (c_.length() <= 0.0)
    throw std::invalid_argument("In UnitCell3D::ctor: c must be non null vector but c = " + std::to_string(c_));
}

UnitCell3D UnitCell3D::create_triclinic_primitive(
  const double a_, const double b_, const double c_, const double alpha_, const double beta_, const double gamma_)
{
  CHECK_LENGTH(a_, b_, c_);
  CHECK_ANGLE(alpha_, beta_, gamma_);

  const Vector3D v1 = x;

  const Vector3D v2 = Rotation(gamma_ * z) * x;

  const double v3x = cos(beta_);
  const double v3y = ( cos(alpha_)-cos(beta_)*cos(gamma_) )/sin(gamma_);
  const double v3z = sqrt(1.0 - v3x*v3x - v3y*v3y);
  const Vector3D v3 = Vector3D( v3x , v3y , v3z );

  return UnitCell3D(Triclinic_Primitive, a_ * v1, b_ * v2, c_ * v3);
}
UnitCell3D UnitCell3D::create_monoclinic_primitive(
  const double a_, const double b_, const double c_, const double beta_)
{
  CHECK_LENGTH(a_, b_, c_);
  CHECK_ANGLE_ACUTE(beta_);

  const Vector3D v1 = a_*x;
  const Vector3D v2 = b_*y;
  const Vector3D v3 = c_* (Rotation(beta_ * y)*z);

  return UnitCell3D(Monoclinic_Primitive, v1, v2, v3);
}
UnitCell3D UnitCell3D::create_monoclinic_base(
  const double a_, const double b_, const double c_, const double beta_)
{
  CHECK_LENGTH(a_, b_, c_);
  CHECK_ANGLE_ACUTE(beta_);

  const Vector3D v1 = a_*x;
  const Vector3D v2 = 1.0 / 2.0 * (a_*x + b_*y);
  const Vector3D v3 = c_*(Rotation(beta_ * y)*z);

  return UnitCell3D(Monoclinic_Base, v1, v2, v3);
}
UnitCell3D UnitCell3D::create_orthorhombic_primitive(const double a_, const double b_, const double c_)
{
  CHECK_LENGTH(a_, b_, c_);

  const Vector3D v1 = a_*x;
  const Vector3D v2 = b_*y;
  const Vector3D v3 = c_*z;

  return UnitCell3D(Orthorhombic_Primitive, v1, v2, v3);
}
UnitCell3D UnitCell3D::create_orthorhombic_base(const double a_, const double b_, const double c_)
{
  CHECK_LENGTH(a_, b_, c_);

  const Vector3D v1 = a_*x;
  const Vector3D v2 = 1.0 / 2.0 * (a_*x + b_*y);
  const Vector3D v3 = c_*z;

  return UnitCell3D(Orthorhombic_Base, v1, v2, v3);
}
UnitCell3D UnitCell3D::create_orthorhombic_body(const double a_, const double b_, const double c_)
{
  CHECK_LENGTH(a_, b_, c_);

  const Vector3D v1 = a_*x;
  const Vector3D v2 = b_*y;
  const Vector3D v3 = (a_*x) + (b_*y) + 1.0 / 2.0 * (c_*z);

  return UnitCell3D(Orthorhombic_Body, v1, v2, v3);
}
UnitCell3D UnitCell3D::create_orthorhombic_face(const double a_, const double b_, const double c_)
{
  CHECK_LENGTH(a_, b_, c_);

  const Vector3D v1 = 1.0 / 2.0 * (a_*x + b_*y);
  const Vector3D v2 = 1.0 / 2.0 * (b_*y + c_*z);
  const Vector3D v3 = 1.0 / 2.0 * (a_*x + c_*z);

  return UnitCell3D(Orthorhombic_Face, v1, v2, v3);
}
UnitCell3D UnitCell3D::create_tetragonal_primitive(const double a_, const double c_)
{
  CHECK_LENGTH2(a_, c_);

  const Vector3D v1 = a_*x;
  const Vector3D v2 = a_*y;
  const Vector3D v3 = c_*z;

  return UnitCell3D(Tetragonal_Primitive, v1, v2, v3);
}
UnitCell3D UnitCell3D::create_tetragonal_body(const double a_, const double c_)
{
  CHECK_LENGTH2(a_, c_);

  const Vector3D v1 = a_*x;
  const Vector3D v2 = a_*y;
  const Vector3D v3 = a_*(x+y) + 1.0 / 2.0 * (c_*z);

  return UnitCell3D(Tetragonal_Body, v1, v2, v3);
}
UnitCell3D UnitCell3D::create_rhombohedral_centered(const double a_, const double alpha_)
{
  CHECK_LENGTH_POSITIVE(a_);
  CHECK_ANGLE_ACUTE(alpha_);

  const Vector3D v1 = x;
  const Vector3D v2 = Vector3D(cos(alpha_), sin(alpha_), 0.0);
  const Vector3D v3 = Vector3D(cos(alpha_)*cos(alpha_/2.0), cos(alpha_)*sin(alpha_/2.0) , sin(alpha_) );

  return UnitCell3D(Rhombohedral_Centered, a_*v1, a_*v2, a_*v3);
}
UnitCell3D UnitCell3D::create_hexagonal_primitive(const double a_, const double c_)
{
  CHECK_LENGTH2(a_, c_);

  const Vector3D v1 = x;
  const Vector3D v2 = Vector3D(cos(pi/3.0), sin(pi/3.0), 0.0);
  const Vector3D v3 = z;

  return UnitCell3D(Hexagonal_Primitive, a_*v1, a_*v2, c_*v3);
}
UnitCell3D UnitCell3D::create_cubic_primitive(const double a_)
{
  CHECK_LENGTH_POSITIVE(a_);

  const Vector3D v1 = a_*x;
  const Vector3D v2 = a_*y;
  const Vector3D v3 = a_*z;

  return UnitCell3D(Cubic_Primitive, v1, v2, v3);
}
UnitCell3D UnitCell3D::create_cubic_body(const double a_)
{
  CHECK_LENGTH_POSITIVE(a_);

  const Vector3D v1 = a_*x;
  const Vector3D v2 = a_*y;
  const Vector3D v3 = a_*(x + y + 1.0/2.0*z);

  return UnitCell3D(Cubic_Body, v1, v2, v3);
}
UnitCell3D UnitCell3D::create_cubic_face(const double a_)
{
  CHECK_LENGTH_POSITIVE(a_);

  const Vector3D v1 = a_ / 2.0 * (x+y);
  const Vector3D v2 = a_ / 2.0 * (y+z);
  const Vector3D v3 = a_ / 2.0 * (x+z);

  return UnitCell3D(Cubic_Face, v1, v2, v3);
}

namespace
{
  long get_steps(const Point3D& unit_vector, const Point3D& perp1, const Point3D& perp2, const Point3D to_digest)
  {
    Point3D parallel = cross_product(perp1, perp2);
    const double parallel_length = parallel.length();
    parallel = 1.0 / parallel_length * parallel;
    
    double parallel_to_digest = to_digest * parallel; 
    double parallel_unit_ = unit_vector * parallel; 

    double ratio = parallel_to_digest / parallel_unit_;
    long l_ratio = std::lround(ratio);

    if (! (equalsWithTolerance(ratio, l_ratio) || nearlyZero( ratio - (double)l_ratio)))
      throw std::invalid_argument("Point " + std::to_string(to_digest) + " is not a lattice vector. Unit vectors are unit_vector = "
        + std::to_string(unit_vector) + " perp1 = "
        + std::to_string(perp1) + " perp2 = "
        + std::to_string(perp2) + " diff for n = " + std::to_string(ratio - (double)l_ratio)
        ); 

    return l_ratio;
  }
}

std::tuple<long, long, long> UnitCell3D::get_offsets(const Point3D& point) const
{
  long lna, lnb, lnc;
  lna = get_steps(a, b, c, point);
  lnb = get_steps(b, c, a, point);
  lnc = get_steps(c, a, b, point);

  return std::make_tuple(lna, lnb, lnc);
}
