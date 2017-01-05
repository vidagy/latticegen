#ifndef LATTICEGEN_BRAVAISLATTICE3D_H_
#define LATTICEGEN_BRAVAISLATTICE3D_H_

#include "Point3D.h"
#include <vector>

namespace Core
{
  namespace Geometry
  {
    class BravaisLattice3D
    {
    public:
      typedef std::vector<Point3D> Point3DVec;

      enum BravaisLattice3DType 
      {
        Triclinic,
        Monoclinic_Primitive,
        Monoclinic_Base,
        Orthorhombic_Primitive,
        Orthorhombic_Base,
        Orthorhombic_Body,
        Orthorhombic_Face,
        Tetragonal_Primitive,
        Tetragonal_Body,
        Rhombohedral,
        Hexagonal,
        Cubic_Primitive,
        Cubic_Body,
        Cubic_Face
      };
      
      struct UnitCell 
      {
        static UnitCell create_triclinic(const double a_, const double b_, const double c_, const double alpha_, const double beta_, const double gamma_);
        static UnitCell create_monoclinic_primitive(const double a_, const double b_, const double c_, const double beta_);
        static UnitCell create_monoclinic_base(const double a_, const double b_, const double c_, const double beta_);
        static UnitCell create_orthorhombic_primitive(const double a_, const double b_, const double c_);
        static UnitCell create_orthorhombic_base(const double a_, const double b_, const double c_);
        static UnitCell create_orthorhombic_body(const double a_, const double b_, const double c_);
        static UnitCell create_orthorhombic_face(const double a_, const double b_, const double c_);
        static UnitCell create_tetragonal_primitive(const double a_, const double c_);
        static UnitCell create_tetragonal_body(const double a_, const double c_);
        static UnitCell create_rhombohedral(const double a_, const double alpha_);
        static UnitCell create_hexagonal(const double a_, const double c_);
        static UnitCell create_cubic_primitive(const double a_);
        static UnitCell create_cubic_body(const double a_);
        static UnitCell create_cubic_face(const double a_);

        const BravaisLattice3DType type;

        const Point3D a; 
        const Point3D b;
        const Point3D c;

      private: 
        UnitCell(const BravaisLattice3DType& type_, const Point3D& a_, const Point3D& b_, const Point3D& c_);

        friend BravaisLattice3D;
      };

      BravaisLattice3D(const UnitCell& unit_cell_, const size_t x_width, const size_t y_width, const size_t z_width);

      UnitCell   get_unit_cell() const { return unit_cell; };
      size_t     get_x_width()   const { return x_width; }
      size_t     get_y_width()   const { return y_width; }
      size_t     get_z_width()   const { return z_width; }
      Point3DVec get_lattice()   const { return lattice; }

    private:
      UnitCell unit_cell;
      size_t x_width;
      size_t y_width;
      size_t z_width;
      Point3DVec lattice;
    };
  }
}

#endif // LATTICEGEN_BRAVAISLATTICE3D_H_
