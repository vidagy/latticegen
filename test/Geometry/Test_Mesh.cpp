#include <gtest/gtest.h>
#include <gmock/gmock.h>

#include <Geometry/Mesh.h>

using namespace Core;
using namespace Geometry;

namespace
{
  static const double pi = 3.14159265358979323846;
}

TEST(TestMesh,LatticeMeshContainsOrigin)
{
  UnitCell3D unit_cell = UnitCell3D::create_cubic_primitive(1.0);
  LatticeMesh mesh = LatticeMesh(unit_cell);

  CutoffCube cutoff = CutoffCube(0.1);

  const static std::vector< Point3D > reference = { {0.0,0.0,0.0} };

  EXPECT_THAT( mesh.generate(cutoff), ::testing::ContainerEq(reference) );
}

TEST(TestMesh,LatticeMeshContainsCube)
{
  UnitCell3D unit_cell = UnitCell3D::create_cubic_primitive(1.0);
  CutoffCube cutoff = CutoffCube(1.0);
  LatticeMesh mesh = LatticeMesh(unit_cell);

  const static std::vector< Point3D > reference =
    {
      {-1.0,-1.0,-1.0},{ 0.0,-1.0,-1.0},{ 1.0,-1.0,-1.0},
      {-1.0, 0.0,-1.0},{ 0.0, 0.0,-1.0},{ 1.0, 0.0,-1.0},
      {-1.0, 1.0,-1.0},{ 0.0, 1.0,-1.0},{ 1.0, 1.0,-1.0},

      {-1.0,-1.0, 0.0},{ 0.0,-1.0, 0.0},{ 1.0,-1.0, 0.0},
      {-1.0, 0.0, 0.0},{ 0.0, 0.0, 0.0},{ 1.0, 0.0, 0.0},
      {-1.0, 1.0, 0.0},{ 0.0, 1.0, 0.0},{ 1.0, 1.0, 0.0},

      {-1.0,-1.0, 1.0},{ 0.0,-1.0, 1.0},{ 1.0,-1.0, 1.0},
      {-1.0, 0.0, 1.0},{ 0.0, 0.0, 1.0},{ 1.0, 0.0, 1.0},
      {-1.0, 1.0, 1.0},{ 0.0, 1.0, 1.0},{ 1.0, 1.0, 1.0}
    };

  EXPECT_THAT( mesh.generate(cutoff), ::testing::ContainerEq(reference) );
}

namespace
{
  const static double a = 0.5;
  const static double b = sqrt(3.0) / 2.0;

  const static double c = sqrt(3.0) / 6.0;
  const static double d = 1.0 / sqrt(3.0);
  const static double h = sqrt(2.0 / 3.0);
}

TEST(TestMesh,LatticeMeshContainsSphere)
{
  CutoffSphere cutoff = CutoffSphere(1.0);
  UnitCell3D unit_cell = UnitCell3D::create_rhombohedral_centered(1.0, pi / 3.0);
  LatticeMesh mesh = LatticeMesh(unit_cell);

  const static std::vector< Point3D > reference =
    {
      { 0.0, 0.0, 0.0},
      { 1.0, 0.0, 0.0}, { a, b, 0.0}, { a,-b, 0.0},
      {-1.0, 0.0, 0.0}, {-a, b, 0.0}, {-a,-b, 0.0},

      {0.0,-d, h}, { a, c, h}, {-a, c, h},
      {0.0, d,-h}, { a,-c,-h}, {-a,-c,-h}
    };

  EXPECT_THAT( mesh.generate(cutoff), ::testing::UnorderedElementsAreArray(reference) );
}

TEST(TestMesh,CubicMeshContainsSphere)
{
  CutoffSphere cutoff = CutoffSphere(sqrt(2.0));
  CubicMesh mesh = CubicMesh(1.0);

  const static std::vector< Point3D > reference =
    {
      { 0.0,-1.0,-1.0},
      {-1.0, 0.0,-1.0},{ 0.0, 0.0,-1.0},{ 1.0, 0.0,-1.0},
      { 0.0, 1.0,-1.0},

      {-1.0,-1.0, 0.0},{ 0.0,-1.0, 0.0},{ 1.0,-1.0, 0.0},
      {-1.0, 0.0, 0.0},{ 0.0, 0.0, 0.0},{ 1.0, 0.0, 0.0},
      {-1.0, 1.0, 0.0},{ 0.0, 1.0, 0.0},{ 1.0, 1.0, 0.0},

      { 0.0,-1.0, 1.0},
      {-1.0, 0.0, 1.0},{ 0.0, 0.0, 1.0},{ 1.0, 0.0, 1.0},
      { 0.0, 1.0, 1.0}
    };

  EXPECT_THAT( mesh.generate(cutoff), ::testing::UnorderedElementsAreArray(reference) );
}

TEST(TestMesh,TetrahedronMeshContainsSphere)
{
  CutoffSphere cutoff = CutoffSphere(1.0);
  TetrahedronMesh mesh = TetrahedronMesh(1.0);

  const static std::vector< Point3D > reference =
    {
      { 0.0, 0.0, 0.0},
      { 1.0, 0.0, 0.0}, { a, b, 0.0}, { a,-b, 0.0},
      {-1.0, 0.0, 0.0}, {-a, b, 0.0}, {-a,-b, 0.0},

      {0.0,-d, h}, { a, c, h}, {-a, c, h},
      {0.0, d,-h}, { a,-c,-h}, {-a,-c,-h}
    };

  EXPECT_THAT( mesh.generate(cutoff), ::testing::UnorderedElementsAreArray(reference) );
}
