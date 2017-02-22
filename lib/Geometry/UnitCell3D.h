#ifndef LATTICEGEN_UNITCELL3D_H_
#define LATTICEGEN_UNITCELL3D_H_

#include <Core/Point3D.h>
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

  class UnitCell3D
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

    CrystalClass get_point_group() const;

    std::tuple<long, long, long> get_offsets(const Point3D& point) const;

    const BravaisLattice3DType type;

    const Point3D a; 
    const Point3D b;
    const Point3D c;

  private:
    UnitCell3D(BravaisLattice3DType type_, const Point3D &a_, const Point3D &b_, const Point3D &c_);
  };
}

#endif // LATTICEGEN_UNITCELL3D_H_
