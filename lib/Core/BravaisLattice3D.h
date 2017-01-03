#ifndef LATTICEGEN_BRAVAISLATTICE3D_H_
#define LATTICEGEN_BRAVAISLATTICE3D_H_

#include "Point3D.h"
#include <utility>
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
      
      // canonical unit cell: first cell vector is (a, 0.0), second is (x,y) where x > 0.0 and y >= 0.0 and (x^2 + y^2) < a^2
      struct UnitCell 
      {
        Point3D a; 
        Point3D b;
        Point3D c;

      private: 
        UnitCell(const Point3D& a_, const Point3D& b_, const Point3D& c_);

        friend BravaisLattice3D;
      };

      static UnitCell get_unit_cell(Point3D x, Point3D y, Point3D z);
      // static BravaisLattice3DType find_lattice_type(const UnitCell& unit_cell);
      // static Point3DVec get_irreducible_wedge(
      //   const UnitCell& unit_cell, const unsigned int xsample, const unsigned int ysample, const unsigned int zsample
      //   );
      
      // static std::tuple<Point3D,Point3D,Point3D> get_inverse_unit_cell(const UnitCell& unitCell);
      
      BravaisLattice3D(
        const Point3D& a, const Point3D& b, const Point3D& c, 
        const size_t x_width, const size_t y_width, const size_t z_width
        );

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
