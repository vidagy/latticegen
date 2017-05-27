#ifndef LATTICEGEN_UNITCELL3D_H_
#define LATTICEGEN_UNITCELL3D_H_

#include <Core/Point3D.h>
#include <Core/Exceptions.h>
#include "CrystallographicPointGroups.h"
#include <vector>

using namespace Core;

namespace Geometry
{
  enum BravaisLattice3DType
  {
    Triclinic_Primitive,
    Monoclinic_Primitive,
    Monoclinic_Base,
    Orthorhombic_Primitive,
    Orthorhombic_Base,
    Orthorhombic_Body,
    Orthorhombic_Face,
    Tetragonal_Primitive,
    Tetragonal_Body,
    Rhombohedral_Centered,
    Hexagonal_Primitive,
    Cubic_Primitive,
    Cubic_Body,
    Cubic_Face
  };

  struct Coordinates3D
  {
    Coordinates3D(long a_, long b_, long c_)
      : a(a_), b(b_), c(c_) {}

    bool operator==(const Coordinates3D &other) const
    {
      return (this->a == other.a) && (this->b == other.b) && (this->c == other.c);
    }

    Coordinates3D operator-(const Coordinates3D &other) const
    {
      return Coordinates3D(a - other.a, b - other.b, c - other.c);
    }

    long a;
    long b;
    long c;
  };

  class Cell3D
  {
  public:
    std::tuple<long, long, long> get_offsets(const Point3D &point) const;

    CrystalClass get_point_group() const;

    const BravaisLattice3DType type;

    const Point3D v1;
    const Point3D v2;
    const Point3D v3;

    Point3D at(const Coordinates3D &coordinates) const
    {
      return coordinates.a * v1 + coordinates.b * v2 + coordinates.c * v3;
    }
  protected:
    Cell3D(BravaisLattice3DType type_, const Point3D &v1_, const Point3D &v2_, const Point3D &v3_)
      : type(type_), v1(v1_), v2(v2_), v3(v3_)
    {
      if (!strictlyPositive(v1_.length()))
        THROW_INVALID_ARGUMENT("In UnitCell3D::ctor: a must be non null vector but a = " + std::to_string(v1_));
      if (!strictlyPositive(v2_.length()))
        THROW_INVALID_ARGUMENT("In UnitCell3D::ctor: b must be non null vector but b = " + std::to_string(v2_));
      if (!strictlyPositive(v3_.length()))
        THROW_INVALID_ARGUMENT("In UnitCell3D::ctor: c must be non null vector but c = " + std::to_string(v3_));
    }
  };

  class UnitCell3D : public Cell3D
  {
  public: 
    static UnitCell3D create_triclinic_primitive(
      double a_, double b_, double c_,
      double alpha_, double beta_, double gamma_);

    static UnitCell3D create_monoclinic_primitive(double a_, double b_, double c_, double beta_);

    static UnitCell3D create_monoclinic_base(double a_, double b_, double c_, double beta_);

    static UnitCell3D create_orthorhombic_primitive(double a_, double b_, double c_);

    static UnitCell3D create_orthorhombic_base(double a_, double b_, double c_);

    static UnitCell3D create_orthorhombic_body(double a_, double b_, double c_);

    static UnitCell3D create_orthorhombic_face(double a_, double b_, double c_);

    static UnitCell3D create_tetragonal_primitive(double a_, double c_);

    static UnitCell3D create_tetragonal_body(double a_, double c_);

    static UnitCell3D create_rhombohedral_centered(double a_, double alpha_);

    static UnitCell3D create_hexagonal_primitive(double a_, double c_);

    static UnitCell3D create_cubic_primitive(double a_);

    static UnitCell3D create_cubic_body(double a_);

    static UnitCell3D create_cubic_face(double a_);

  private:
    UnitCell3D(BravaisLattice3DType type_, const Point3D &a_, const Point3D &b_, const Point3D &c_)
      : Cell3D(type_, a_, b_, c_) {}
  };

  class ReciprocalUnitCell3D : public Cell3D
  {
  public:
    ReciprocalUnitCell3D(const UnitCell3D &unit_cell);
  };
}

namespace std
{
  inline std::string to_string(const Geometry::Coordinates3D &coords)
  {
    return "( " + std::to_string(coords.a) + " , " + std::to_string(coords.b) + " , " + std::to_string(coords.c) + " )";
  }
}

#endif // LATTICEGEN_UNITCELL3D_H_
