#include <gtest/gtest.h>
#include <gmock/gmock.h>

#include <Geometry/Mesh.h>

using namespace Core;
using namespace Geometry;

namespace
{
  static const double pi = 3.14159265358979323846;
}

TEST(TestMesh,MeshContainsOrigin)
{
  CutoffCube cutoff = CutoffCube(0.1);
  UnitCell3D unit_cell = UnitCell3D::create_cubic_primitive(1.0);

  LatticeMesh mesh = LatticeMesh(unit_cell);
  LatticeMesh mesh_pos = LatticeMesh(unit_cell, true);

  const static std::vector< Point3D > reference = { {0.0,0.0,0.0} };

  EXPECT_THAT( mesh.generate(cutoff), ::testing::ContainerEq(reference) );
  EXPECT_THAT( mesh_pos.generate(cutoff), ::testing::ContainerEq(reference) );
}

TEST(TestMesh,MeshContainsCube)
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

TEST(TestMesh,MeshContainsCubePositive)
{
  UnitCell3D unit_cell = UnitCell3D::create_cubic_primitive(1.0);
  CutoffCube cutoff = CutoffCube(1.0);
  LatticeMesh mesh = LatticeMesh(unit_cell, true);

  const static std::vector< Point3D > reference =
    {
      { 0.0, 0.0, 0.0},{ 1.0, 0.0, 0.0},
      { 0.0, 1.0, 0.0},{ 1.0, 1.0, 0.0},
      { 0.0, 0.0, 1.0},{ 1.0, 0.0, 1.0},
      { 0.0, 1.0, 1.0},{ 1.0, 1.0, 1.0}
    };

  EXPECT_THAT( mesh.generate(cutoff), ::testing::ContainerEq(reference) );
}

