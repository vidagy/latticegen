#include "UnitCell3D.h" 

#include "Transformations.h"

#include <cmath>

using namespace Geometry;
using namespace Core;

#define CHECK_LENGTH_POSITIVE(x) \
    if (! strictlyPositive(x))         \
      { THROW_INVALID_ARGUMENT(std::string("Lattice vector length ") + #x  + " must be positive"); } ;

#define CHECK_LENGTH_DIFFERENT(x, y)   \
    if (equalsWithTolerance(x, y))     \
      { THROW_INVALID_ARGUMENT(std::string("Lattice vector length ") + #x  + " and " + #y + " must differ"); } ;

#define CHECK_ANGLE_ACUTE(x)            \
    if (! strictlyPositive(x) || ! strictlyLess(x, pi/2.0) )    \
      THROW_INVALID_ARGUMENT(std::string("Lattice angle ") + #x  + " must be acute");

#define CHECK_ANGLE_DIFFERENT(x, y)   \
    if (equalsWithTolerance(x, y))     \
      { THROW_INVALID_ARGUMENT(std::string("Lattice angle ") + #x  + " and " + #y + " must differ"); } ;

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

UnitCell3D UnitCell3D::create_triclinic_primitive(
  double a_, double b_, double c_, double alpha_, double beta_, double gamma_)
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
  double a_, double b_, double c_, double beta_)
{
  CHECK_LENGTH(a_, b_, c_);
  CHECK_ANGLE_ACUTE(beta_);

  const Vector3D v1 = a_*x;
  const Vector3D v2 = b_*y;
  const Vector3D v3 = c_* (Rotation(beta_ * y)*z);

  return UnitCell3D(Monoclinic_Primitive, v1, v2, v3);
}
UnitCell3D UnitCell3D::create_monoclinic_base(
  double a_, double b_, double c_, double beta_)
{
  CHECK_LENGTH(a_, b_, c_);
  CHECK_ANGLE_ACUTE(beta_);

  const Vector3D v1 = a_*x;
  const Vector3D v2 = 1.0 / 2.0 * (a_*x + b_*y);
  const Vector3D v3 = c_*(Rotation(beta_ * y)*z);

  return UnitCell3D(Monoclinic_Base, v1, v2, v3);
}

UnitCell3D UnitCell3D::create_orthorhombic_primitive(double a_, double b_, double c_)
{
  CHECK_LENGTH(a_, b_, c_);

  const Vector3D v1 = a_*x;
  const Vector3D v2 = b_*y;
  const Vector3D v3 = c_*z;

  return UnitCell3D(Orthorhombic_Primitive, v1, v2, v3);
}

UnitCell3D UnitCell3D::create_orthorhombic_base(double a_, double b_, double c_)
{
  CHECK_LENGTH(a_, b_, c_);

  const Vector3D v1 = a_*x;
  const Vector3D v2 = 1.0 / 2.0 * (a_*x + b_*y);
  const Vector3D v3 = c_*z;

  return UnitCell3D(Orthorhombic_Base, v1, v2, v3);
}

UnitCell3D UnitCell3D::create_orthorhombic_body(double a_, double b_, double c_)
{
  CHECK_LENGTH(a_, b_, c_);

  const Vector3D v1 = a_*x;
  const Vector3D v2 = b_*y;
  const Vector3D v3 = (a_*x) + (b_*y) + 1.0 / 2.0 * (c_*z);

  return UnitCell3D(Orthorhombic_Body, v1, v2, v3);
}

UnitCell3D UnitCell3D::create_orthorhombic_face(double a_, double b_, double c_)
{
  CHECK_LENGTH(a_, b_, c_);

  const Vector3D v1 = 1.0 / 2.0 * (a_*x + b_*y);
  const Vector3D v2 = 1.0 / 2.0 * (b_*y + c_*z);
  const Vector3D v3 = 1.0 / 2.0 * (a_*x + c_*z);

  return UnitCell3D(Orthorhombic_Face, v1, v2, v3);
}

UnitCell3D UnitCell3D::create_tetragonal_primitive(double a_, double c_)
{
  CHECK_LENGTH2(a_, c_);

  const Vector3D v1 = a_*x;
  const Vector3D v2 = a_*y;
  const Vector3D v3 = c_*z;

  return UnitCell3D(Tetragonal_Primitive, v1, v2, v3);
}

UnitCell3D UnitCell3D::create_tetragonal_body(double a_, double c_)
{
  CHECK_LENGTH2(a_, c_);

  const Vector3D v1 = a_*x;
  const Vector3D v2 = a_*y;
  const Vector3D v3 = a_*(x+y) + 1.0 / 2.0 * (c_*z);

  return UnitCell3D(Tetragonal_Body, v1, v2, v3);
}

UnitCell3D UnitCell3D::create_rhombohedral_centered(double a_, double alpha_)
{
  CHECK_LENGTH_POSITIVE(a_);
  CHECK_ANGLE_ACUTE(alpha_);

  const Vector3D v1 = x;
  const Vector3D v2 = Vector3D(cos(alpha_), sin(alpha_), 0.0);
  const double cosBeta = cos(alpha_) / cos(alpha_/2.0);
  const Vector3D v3 = Vector3D(cosBeta*cos(alpha_/2.0), cosBeta*sin(alpha_/2.0) , sqrt(1.0 - pow(cosBeta,2)) );

  return UnitCell3D(Rhombohedral_Centered, a_*v1, a_*v2, a_*v3);
}

UnitCell3D UnitCell3D::create_hexagonal_primitive(double a_, double c_)
{
  CHECK_LENGTH2(a_, c_);

  const Vector3D v1 = x;
  const Vector3D v2 = Vector3D(cos(pi/3.0), sin(pi/3.0), 0.0);
  const Vector3D v3 = z;

  return UnitCell3D(Hexagonal_Primitive, a_*v1, a_*v2, c_*v3);
}

UnitCell3D UnitCell3D::create_cubic_primitive(double a_)
{
  CHECK_LENGTH_POSITIVE(a_);

  const Vector3D v1 = a_*x;
  const Vector3D v2 = a_*y;
  const Vector3D v3 = a_*z;

  return UnitCell3D(Cubic_Primitive, v1, v2, v3);
}

UnitCell3D UnitCell3D::create_cubic_body(double a_)
{
  CHECK_LENGTH_POSITIVE(a_);

  const Vector3D v1 = a_*x;
  const Vector3D v2 = a_*y;
  const Vector3D v3 = a_*(x + y + 1.0/2.0*z);

  return UnitCell3D(Cubic_Body, v1, v2, v3);
}

UnitCell3D UnitCell3D::create_cubic_face(double a_)
{
  CHECK_LENGTH_POSITIVE(a_);

  const Vector3D v1 = a_ / 2.0 * (x+y);
  const Vector3D v2 = a_ / 2.0 * (y+z);
  const Vector3D v3 = a_ / 2.0 * (x+z);

  return UnitCell3D(Cubic_Face, v1, v2, v3);
}

namespace
{
  long get_steps(const Point3D& unit_vector, const Point3D& perp1, const Point3D& perp2, const Point3D& to_digest)
  {
    Point3D parallel = cross_product(perp1, perp2);
    const double parallel_length = parallel.length();
    parallel = 1.0 / parallel_length * parallel;
    
    double parallel_to_digest = to_digest * parallel; 
    double parallel_unit_ = unit_vector * parallel; 

    double ratio = parallel_to_digest / parallel_unit_;
    long l_ratio = std::lround(ratio);

    if (! (equalsWithTolerance(ratio, l_ratio) || nearlyZero( ratio - (double)l_ratio)))
      THROW_INVALID_ARGUMENT(
        "Point " + std::to_string(to_digest) + " is not a lattice vector. Unit vectors are unit_vector = "
        + std::to_string(unit_vector) + " perp1 = "
        + std::to_string(perp1) + " perp2 = "
        + std::to_string(perp2) + " diff for n = " + std::to_string(ratio - (double)l_ratio)
        ); 

    return l_ratio;
  }
}

std::tuple<long, long, long> Cell3D::get_offsets(const Point3D &point) const
{
  long lna, lnb, lnc;
  lna = get_steps(v1, v2, v3, point);
  lnb = get_steps(v2, v3, v1, point);
  lnc = get_steps(v3, v1, v2, point);

  return std::make_tuple(lna, lnb, lnc);
}

CrystalClass Cell3D::get_point_group() const
{
  switch (type)
  {
    case Triclinic_Primitive:
      return CrystalClass::Triclinic_Pinacoid;
    case Monoclinic_Primitive:
    case Monoclinic_Base:
      return CrystalClass::Monoclinic_Prism;
    case Orthorhombic_Primitive:
    case Orthorhombic_Base:
    case Orthorhombic_Body:
    case Orthorhombic_Face:
      return CrystalClass::Orthorhombic_Dipyramid;
    case Tetragonal_Primitive:
    case Tetragonal_Body:
      return CrystalClass::Ditetragonal_Dipyramid;
    case Rhombohedral_Centered:
      return CrystalClass::Ditrigonal_Scalenohedron;
    case Hexagonal_Primitive:
      return CrystalClass::Dihexagonal_Dipyramid;
    case Cubic_Primitive:
    case Cubic_Body:
    case Cubic_Face:
      return CrystalClass::Hexaoctahedron;
    default:
      THROW_INVALID_ARGUMENT("Unhandled BravaisLattice3DType in UnitCell3D::get_crystal_class");
  }
}

ReciprocalUnitCell3D::ReciprocalUnitCell3D(const UnitCell3D &unit_cell)
  : Cell3D(
  unit_cell.type,
  2 * pi * (cross_product(unit_cell.v2, unit_cell.v3)) / (unit_cell.v1 * cross_product(unit_cell.v2, unit_cell.v3)),
  2 * pi * (cross_product(unit_cell.v3, unit_cell.v1)) / (unit_cell.v2 * cross_product(unit_cell.v3, unit_cell.v1)),
  2 * pi * (cross_product(unit_cell.v1, unit_cell.v2)) / (unit_cell.v3 * cross_product(unit_cell.v1, unit_cell.v2))
)
{
}
