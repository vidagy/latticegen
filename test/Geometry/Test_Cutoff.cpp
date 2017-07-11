#include <TestUtils/base.h>

#include <Geometry/Cutoff.h>

using namespace Core;
using namespace Geometry;

namespace
{
    static const double pi = 3.14159265358979323846;
}

TEST(TestCutoff,CubeCtorThrows)
{
  EXPECT_THROW(CutoffCube(0.0), std::invalid_argument);
  EXPECT_THROW(CutoffCube(-1.5), std::invalid_argument);
}

TEST(TestCutoff,CubeCutoffTrue)
{
  CutoffCube cutoff = CutoffCube(1.0);
  EXPECT_TRUE(cutoff.is_included( Point3D( 1.0, 1.0, 1.0) ) );
  EXPECT_TRUE(cutoff.is_included( Point3D( 0.0, 1.0, 1.0) ) );
  EXPECT_TRUE(cutoff.is_included( Point3D( 1.0, 0.0, 1.0) ) );
  EXPECT_TRUE(cutoff.is_included( Point3D( 1.0, 1.0, 0.0) ) );
  EXPECT_TRUE(cutoff.is_included( Point3D( 0.0, 0.0, 1.0) ) );
  EXPECT_TRUE(cutoff.is_included( Point3D( 0.0, 1.0, 0.0) ) );
  EXPECT_TRUE(cutoff.is_included( Point3D( 1.0, 0.0, 0.0) ) );
  EXPECT_TRUE(cutoff.is_included( Point3D( 0.0, 0.0, 0.0) ) );

  EXPECT_TRUE(cutoff.is_included( Point3D(-1.0,-1.0,-1.0) ) );
  EXPECT_TRUE(cutoff.is_included( Point3D( 0.0, 1.0,-1.0) ) );
  EXPECT_TRUE(cutoff.is_included( Point3D(-1.0, 0.0, 1.0) ) );
  EXPECT_TRUE(cutoff.is_included( Point3D( 1.0,-1.0, 0.0) ) );
  EXPECT_TRUE(cutoff.is_included( Point3D( 0.0, 0.0,-1.0) ) );
  EXPECT_TRUE(cutoff.is_included( Point3D( 0.0,-1.0, 0.0) ) );
  EXPECT_TRUE(cutoff.is_included( Point3D(-1.0, 0.0, 0.0) ) );
}

TEST(TestCutoff,CubeCutoffFalse)
{
  CutoffCube cutoff = CutoffCube(1.0);
  EXPECT_FALSE(cutoff.is_included( Point3D( 1.1, 1.1, 1.1) ) );
  EXPECT_FALSE(cutoff.is_included( Point3D( 1.0, 1.1, 1.1) ) );
  EXPECT_FALSE(cutoff.is_included( Point3D( 1.1, 1.0, 1.1) ) );
  EXPECT_FALSE(cutoff.is_included( Point3D( 1.1, 1.1, 1.0) ) );
  EXPECT_FALSE(cutoff.is_included( Point3D( 1.0, 1.0, 1.1) ) );
  EXPECT_FALSE(cutoff.is_included( Point3D( 1.0, 1.1, 1.0) ) );
  EXPECT_FALSE(cutoff.is_included( Point3D( 1.1, 1.0, 1.0) ) );

  EXPECT_FALSE(cutoff.is_included( Point3D(-1.1,-1.1,-1.1) ) );
  EXPECT_FALSE(cutoff.is_included( Point3D( 1.0, 1.1,-1.1) ) );
  EXPECT_FALSE(cutoff.is_included( Point3D(-1.1, 1.0, 1.1) ) );
  EXPECT_FALSE(cutoff.is_included( Point3D( 1.1,-1.1, 1.0) ) );
  EXPECT_FALSE(cutoff.is_included( Point3D( 1.0, 1.0,-1.1) ) );
  EXPECT_FALSE(cutoff.is_included( Point3D( 1.0,-1.1, 1.0) ) );
  EXPECT_FALSE(cutoff.is_included( Point3D(-1.1, 1.0, 1.0) ) );
}

TEST(TestCutoff,SphereCtorThrows)
{
  EXPECT_THROW(CutoffSphere(0.0), std::invalid_argument);
  EXPECT_THROW(CutoffSphere(-1.5), std::invalid_argument);
}

TEST(TestCutoff,SphereCutoffTrue)
{
  CutoffSphere cutoff = CutoffSphere(1.0);
  EXPECT_TRUE(cutoff.is_included( Point3D( 0.0, 0.0, 0.0) ) );
  EXPECT_TRUE(cutoff.is_included( Point3D( 1.0, 0.0, 0.0) ) );
  EXPECT_TRUE(cutoff.is_included( Point3D( 0.0, 1.0, 0.0) ) );
  EXPECT_TRUE(cutoff.is_included( Point3D( 0.0, 0.0, 1.0) ) );
  EXPECT_TRUE(cutoff.is_included( Point3D(-1.0, 0.0, 0.0) ) );
  EXPECT_TRUE(cutoff.is_included( Point3D( 0.0,-1.0, 0.0) ) );
  EXPECT_TRUE(cutoff.is_included( Point3D( 0.0, 0.0,-1.0) ) );

  EXPECT_TRUE(cutoff.is_included( Point3D( 0.5, 0.5, 0.0) ) );
  EXPECT_TRUE(cutoff.is_included( Point3D(-0.5, 0.5, 0.0) ) );
  EXPECT_TRUE(cutoff.is_included( Point3D( 0.5,-0.5, 0.0) ) );
  EXPECT_TRUE(cutoff.is_included( Point3D(-0.5,-0.5, 0.0) ) );

  EXPECT_TRUE(cutoff.is_included( Point3D( 0.5, 0.0, 0.5) ) );
  EXPECT_TRUE(cutoff.is_included( Point3D(-0.5, 0.0, 0.5) ) );
  EXPECT_TRUE(cutoff.is_included( Point3D( 0.5, 0.0,-0.5) ) );
  EXPECT_TRUE(cutoff.is_included( Point3D(-0.5, 0.0,-0.5) ) );

  EXPECT_TRUE(cutoff.is_included( Point3D( 0.0, 0.5, 0.5) ) );
  EXPECT_TRUE(cutoff.is_included( Point3D( 0.0,-0.5, 0.5) ) );
  EXPECT_TRUE(cutoff.is_included( Point3D( 0.0, 0.5,-0.5) ) );
  EXPECT_TRUE(cutoff.is_included( Point3D( 0.0,-0.5,-0.5) ) );

  EXPECT_TRUE(cutoff.is_included( Point3D( 0.5, 0.5, 0.5) ) );
  EXPECT_TRUE(cutoff.is_included( Point3D(-0.5, 0.5, 0.5) ) );
  EXPECT_TRUE(cutoff.is_included( Point3D( 0.5,-0.5, 0.5) ) );
  EXPECT_TRUE(cutoff.is_included( Point3D( 0.5, 0.5,-0.5) ) );
  EXPECT_TRUE(cutoff.is_included( Point3D(-0.5,-0.5, 0.5) ) );
  EXPECT_TRUE(cutoff.is_included( Point3D(-0.5, 0.5,-0.5) ) );
  EXPECT_TRUE(cutoff.is_included( Point3D( 0.5,-0.5,-0.5) ) );
  EXPECT_TRUE(cutoff.is_included( Point3D(-0.5,-0.5,-0.5) ) );
}

TEST(TestCutoff,SphereCutoffFalse)
{
  CutoffSphere cutoff = CutoffSphere(1.0);
  EXPECT_FALSE(cutoff.is_included( Point3D( 1.1, 0.0, 0.0) ) );
  EXPECT_FALSE(cutoff.is_included( Point3D( 0.0, 1.1, 0.0) ) );
  EXPECT_FALSE(cutoff.is_included( Point3D( 0.0, 0.0, 1.1) ) );
  EXPECT_FALSE(cutoff.is_included( Point3D(-1.1, 0.0, 0.0) ) );
  EXPECT_FALSE(cutoff.is_included( Point3D( 0.0,-1.1, 0.0) ) );
  EXPECT_FALSE(cutoff.is_included( Point3D( 0.0, 0.0,-1.1) ) );

  EXPECT_FALSE(cutoff.is_included( Point3D( 0.75, 0.75, 0.0) ) );
  EXPECT_FALSE(cutoff.is_included( Point3D(-0.75, 0.75, 0.0) ) );
  EXPECT_FALSE(cutoff.is_included( Point3D( 0.75,-0.75, 0.0) ) );
  EXPECT_FALSE(cutoff.is_included( Point3D(-0.75,-0.75, 0.0) ) );

  EXPECT_FALSE(cutoff.is_included( Point3D( 0.75, 0.0, 0.75) ) );
  EXPECT_FALSE(cutoff.is_included( Point3D(-0.75, 0.0, 0.75) ) );
  EXPECT_FALSE(cutoff.is_included( Point3D( 0.75, 0.0,-0.75) ) );
  EXPECT_FALSE(cutoff.is_included( Point3D(-0.75, 0.0,-0.75) ) );

  EXPECT_FALSE(cutoff.is_included( Point3D( 0.0, 0.75, 0.75) ) );
  EXPECT_FALSE(cutoff.is_included( Point3D( 0.0,-0.75, 0.75) ) );
  EXPECT_FALSE(cutoff.is_included( Point3D( 0.0, 0.75,-0.75) ) );
  EXPECT_FALSE(cutoff.is_included( Point3D( 0.0,-0.75,-0.75) ) );

  EXPECT_FALSE(cutoff.is_included( Point3D( 0.75, 0.75, 0.75) ) );
  EXPECT_FALSE(cutoff.is_included( Point3D(-0.75, 0.75, 0.75) ) );
  EXPECT_FALSE(cutoff.is_included( Point3D( 0.75,-0.75, 0.75) ) );
  EXPECT_FALSE(cutoff.is_included( Point3D( 0.75, 0.75,-0.75) ) );
  EXPECT_FALSE(cutoff.is_included( Point3D(-0.75,-0.75, 0.75) ) );
  EXPECT_FALSE(cutoff.is_included( Point3D(-0.75, 0.75,-0.75) ) );
  EXPECT_FALSE(cutoff.is_included( Point3D( 0.75,-0.75,-0.75) ) );
  EXPECT_FALSE(cutoff.is_included( Point3D(-0.75,-0.75,-0.75) ) );
}

TEST(TestCutoff,UnitVectorsCutoffTrue)
{
  UnitCell3D unit_cell = UnitCell3D::create_rhombohedral_centered(1.0, pi / 3.0);
  const Vector3D &a = unit_cell.v1;
  const Vector3D &b = unit_cell.v2;
  const Vector3D &c = unit_cell.v3;

  CutoffUnitVectors cutoff = CutoffUnitVectors(unit_cell, 1, 2, 3);

  EXPECT_TRUE(cutoff.is_included( Point3D( 0.0, 0.0, 0.0) ) );
  EXPECT_TRUE(cutoff.is_included( a ) );
  EXPECT_TRUE(cutoff.is_included( b ) );
  EXPECT_TRUE(cutoff.is_included( c ) );
  EXPECT_TRUE(cutoff.is_included( -1.0 * a ) );
  EXPECT_TRUE(cutoff.is_included( -1.0 * b ) );
  EXPECT_TRUE(cutoff.is_included( -1.0 * c ) );

  EXPECT_TRUE(cutoff.is_included( 1.0 * a + 1.0 * b ) );
  EXPECT_TRUE(cutoff.is_included(-1.0 * a + 1.0 * b ) );
  EXPECT_TRUE(cutoff.is_included( 1.0 * a - 1.0 * b ) );
  EXPECT_TRUE(cutoff.is_included(-1.0 * a - 1.0 * b ) );

  EXPECT_TRUE(cutoff.is_included( 1.0 * a + 1.0 * c ) );
  EXPECT_TRUE(cutoff.is_included(-1.0 * a + 1.0 * c ) );
  EXPECT_TRUE(cutoff.is_included( 1.0 * a - 1.0 * c ) );
  EXPECT_TRUE(cutoff.is_included(-1.0 * a - 1.0 * c ) );

  EXPECT_TRUE(cutoff.is_included( 1.0 * b + 1.0 * c ) );
  EXPECT_TRUE(cutoff.is_included(-1.0 * b + 1.0 * c ) );
  EXPECT_TRUE(cutoff.is_included( 1.0 * b - 1.0 * c ) );
  EXPECT_TRUE(cutoff.is_included(-1.0 * b - 1.0 * c ) );

  EXPECT_TRUE(cutoff.is_included( 1.0 * a + 2.0 * b + 3.0 * c ) );
  EXPECT_TRUE(cutoff.is_included(-1.0 * a + 2.0 * b + 3.0 * c ) );
  EXPECT_TRUE(cutoff.is_included( 1.0 * a - 2.0 * b + 3.0 * c ) );
  EXPECT_TRUE(cutoff.is_included( 1.0 * a + 2.0 * b - 3.0 * c ) );
  EXPECT_TRUE(cutoff.is_included(-1.0 * a - 2.0 * b + 3.0 * c ) );
  EXPECT_TRUE(cutoff.is_included(-1.0 * a + 2.0 * b - 3.0 * c ) );
  EXPECT_TRUE(cutoff.is_included( 1.0 * a - 2.0 * b - 3.0 * c ) );
  EXPECT_TRUE(cutoff.is_included(-1.0 * a - 2.0 * b - 3.0 * c ) );
}

TEST(TestCutoff,UnitVectorsCutoffFalse)
{
  UnitCell3D unit_cell = UnitCell3D::create_rhombohedral_centered(1.0, pi / 3.0);
  const Vector3D &a = unit_cell.v1;
  const Vector3D &b = unit_cell.v2;
  const Vector3D &c = unit_cell.v3;

  CutoffUnitVectors cutoff = CutoffUnitVectors(unit_cell, 1, 2, 3);

  EXPECT_FALSE(cutoff.is_included( 2.0 * a ) );
  EXPECT_FALSE(cutoff.is_included( 3.0 * b ) );
  EXPECT_FALSE(cutoff.is_included( 4.0 * c ) );
  EXPECT_FALSE(cutoff.is_included( -2.0 * a ) );
  EXPECT_FALSE(cutoff.is_included( -3.0 * b ) );
  EXPECT_FALSE(cutoff.is_included( -4.0 * c ) );

  EXPECT_FALSE(cutoff.is_included( 2.0 * a + 2.0 * b + 3.0 * c ) );
  EXPECT_FALSE(cutoff.is_included(-2.0 * a + 2.0 * b + 3.0 * c ) );
  EXPECT_FALSE(cutoff.is_included( 2.0 * a - 2.0 * b + 3.0 * c ) );
  EXPECT_FALSE(cutoff.is_included( 2.0 * a + 2.0 * b - 3.0 * c ) );

  EXPECT_FALSE(cutoff.is_included( 1.0 * a + 3.0 * b + 3.0 * c ) );
  EXPECT_FALSE(cutoff.is_included(-1.0 * a + 3.0 * b + 3.0 * c ) );
  EXPECT_FALSE(cutoff.is_included( 1.0 * a - 3.0 * b + 3.0 * c ) );
  EXPECT_FALSE(cutoff.is_included( 1.0 * a + 3.0 * b - 3.0 * c ) );

  EXPECT_FALSE(cutoff.is_included( 1.0 * a + 2.0 * b + 4.0 * c ) );
  EXPECT_FALSE(cutoff.is_included(-1.0 * a + 2.0 * b + 4.0 * c ) );
  EXPECT_FALSE(cutoff.is_included( 1.0 * a - 2.0 * b + 4.0 * c ) );
  EXPECT_FALSE(cutoff.is_included( 1.0 * a + 2.0 * b - 4.0 * c ) );
}

TEST(TestCutoff,WSCellCutoffTrue)
{
  UnitCell3D unit_cell = UnitCell3D::create_cubic_primitive(1.0);
  const Vector3D &a = unit_cell.v1;
  const Vector3D &b = unit_cell.v2;
  const Vector3D &c = unit_cell.v3;

  CutoffWSCell cutoff = CutoffWSCell(UnitCell3D::create_cubic_primitive(3.0));

  EXPECT_TRUE(cutoff.is_included( {0,0,0} ));
  EXPECT_TRUE(cutoff.is_included( a ));
  EXPECT_TRUE(cutoff.is_included(-a ));
  EXPECT_TRUE(cutoff.is_included( b ));
  EXPECT_TRUE(cutoff.is_included(-b ));
  EXPECT_TRUE(cutoff.is_included( c ));
  EXPECT_TRUE(cutoff.is_included(-c ));
  EXPECT_TRUE(cutoff.is_included( a +b ));
  EXPECT_TRUE(cutoff.is_included(-a +b ));
  EXPECT_TRUE(cutoff.is_included( a -b ));
  EXPECT_TRUE(cutoff.is_included(-a -b ));
  EXPECT_TRUE(cutoff.is_included( a +c ));
  EXPECT_TRUE(cutoff.is_included(-a +c ));
  EXPECT_TRUE(cutoff.is_included( a -c ));
  EXPECT_TRUE(cutoff.is_included(-a -c ));
  EXPECT_TRUE(cutoff.is_included( b +c ));
  EXPECT_TRUE(cutoff.is_included(-b +c ));
  EXPECT_TRUE(cutoff.is_included( b -c ));
  EXPECT_TRUE(cutoff.is_included(-b -c ));
  EXPECT_TRUE(cutoff.is_included( a +b +c ));
  EXPECT_TRUE(cutoff.is_included( a -b +c ));
  EXPECT_TRUE(cutoff.is_included( a +b -c ));
  EXPECT_TRUE(cutoff.is_included( a -b -c ));
  EXPECT_TRUE(cutoff.is_included(-a +b +c ));
  EXPECT_TRUE(cutoff.is_included(-a -b +c ));
  EXPECT_TRUE(cutoff.is_included(-a +b -c ));
  EXPECT_TRUE(cutoff.is_included(-a -b -c ));
}

TEST(TestCutoff,WSCellCutoffFalse)
{
  UnitCell3D unit_cell = UnitCell3D::create_cubic_primitive(1.0);
  const Vector3D &a = 2.0 * unit_cell.v1;
  const Vector3D &b = 2.0 * unit_cell.v2;
  const Vector3D &c = 2.0 * unit_cell.v3;

  CutoffWSCell cutoff = CutoffWSCell(UnitCell3D::create_cubic_primitive(3.0));

  EXPECT_FALSE(cutoff.is_included( a ));
  EXPECT_FALSE(cutoff.is_included(-a ));
  EXPECT_FALSE(cutoff.is_included( b ));
  EXPECT_FALSE(cutoff.is_included(-b ));
  EXPECT_FALSE(cutoff.is_included( c ));
  EXPECT_FALSE(cutoff.is_included(-c ));
  EXPECT_FALSE(cutoff.is_included( a +b ));
  EXPECT_FALSE(cutoff.is_included(-a +b ));
  EXPECT_FALSE(cutoff.is_included( a -b ));
  EXPECT_FALSE(cutoff.is_included(-a -b ));
  EXPECT_FALSE(cutoff.is_included( a +c ));
  EXPECT_FALSE(cutoff.is_included(-a +c ));
  EXPECT_FALSE(cutoff.is_included( a -c ));
  EXPECT_FALSE(cutoff.is_included(-a -c ));
  EXPECT_FALSE(cutoff.is_included( b +c ));
  EXPECT_FALSE(cutoff.is_included(-b +c ));
  EXPECT_FALSE(cutoff.is_included( b -c ));
  EXPECT_FALSE(cutoff.is_included(-b -c ));
  EXPECT_FALSE(cutoff.is_included( a +b +c ));
  EXPECT_FALSE(cutoff.is_included( a -b +c ));
  EXPECT_FALSE(cutoff.is_included( a +b -c ));
  EXPECT_FALSE(cutoff.is_included( a -b -c ));
  EXPECT_FALSE(cutoff.is_included(-a +b +c ));
  EXPECT_FALSE(cutoff.is_included(-a -b +c ));
  EXPECT_FALSE(cutoff.is_included(-a +b -c ));
  EXPECT_FALSE(cutoff.is_included(-a -b -c ));
}

TEST(TestCutoff, WSCellRadius)
{
  //CutoffWSCell cutoff_cube = CutoffWSCell(UnitCell3D::create_cubic_primitive(1.0));
  //EXPECT_DOUBLE_EQ(cutoff_cube.r_in_for_sure, 0.5);
  //EXPECT_DOUBLE_EQ(cutoff_cube.r_out_for_sure, 0.5*sqrt(3.0));

  CutoffWSCell cutoff_hexagonal_short = CutoffWSCell(UnitCell3D::create_hexagonal_primitive(1.0, 1e-10));
  EXPECT_DOUBLE_EQ(cutoff_hexagonal_short.r_in_for_sure, 5e-11);
  EXPECT_DOUBLE_EQ(cutoff_hexagonal_short.r_out_for_sure, 1.0 / sqrt(3.0));

}
