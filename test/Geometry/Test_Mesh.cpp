#include <TestUtils/base.h>

#include <Geometry/Mesh.h>

using namespace Core;
using namespace Geometry;
using namespace Testing;

TEST(TestMesh,LatticeMeshContainsOrigin)
{
  UnitCell3D unit_cell = UnitCell3D::create_cubic_primitive(1.0);
  auto mesh = LatticeMesh(unit_cell);

  auto cutoff = CutoffCube(0.1);

  const static std::vector< Point3D > reference = { {0.0,0.0,0.0} };

  EXPECT_THAT(wrap(mesh.generate(cutoff)), ::testing::ContainerEq(wrap(reference)));
}

TEST(TestMesh,LatticeMeshContainsCube)
{
  UnitCell3D unit_cell = UnitCell3D::create_cubic_primitive(1.0);
  auto cutoff = CutoffCube(1.0);
  auto mesh = LatticeMesh(unit_cell);

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

  EXPECT_THAT(wrap(mesh.generate(cutoff)), ::testing::ContainerEq(wrap(reference)));
}

namespace
{
  const double a = 0.5;
  const double b = sqrt(3.0) / 2.0;

  const double c = sqrt(3.0) / 6.0;
  const double d = 1.0 / sqrt(3.0);
  const double h = sqrt(2.0 / 3.0);
}

TEST(TestMesh,LatticeMeshContainsSphere)
{
  auto cutoff = CutoffSphere(1.0);
  UnitCell3D unit_cell = UnitCell3D::create_rhombohedral_centered(1.0, pi / 3.0);
  auto mesh = LatticeMesh(unit_cell);

  const static std::vector< Point3D > reference =
    {
      { 0.0, 0.0, 0.0},
      { 1.0, 0.0, 0.0}, { a, b, 0.0}, { a,-b, 0.0},
      {-1.0, 0.0, 0.0}, {-a, b, 0.0}, {-a,-b, 0.0},

      {0.0,-d, h}, { a, c, h}, {-a, c, h},
      {0.0, d,-h}, { a,-c,-h}, {-a,-c,-h}
    };

  EXPECT_THAT(wrap(mesh.generate(cutoff)), ::testing::UnorderedElementsAreArray(wrap(reference)));
}

TEST(TestMesh,CubicMeshContainsSphere)
{
  auto cutoff = CutoffSphere(sqrt(2.0));
  auto mesh = CubicMesh(1.0);

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

  EXPECT_THAT(wrap(mesh.generate(cutoff)), ::testing::UnorderedElementsAreArray(wrap(reference)));
}

TEST(TestMesh, CubicMeshWithOffsetContainsSphere)
{
  auto cutoff = CutoffSphere(sqrt(2.0));
  auto mesh = CubicMesh(1.0);

  const static std::vector<Point3D> reference =
    {
      {-0.5, -0.5, -0.5},
      {0.5,  -0.5, -0.5},
      {-0.5, 0.5,  -0.5},
      {0.5,  0.5,  -0.5},
      {-0.5, -0.5, 0.5},
      {0.5,  -0.5, 0.5},
      {-0.5, 0.5,  0.5},
      {0.5,  0.5,  0.5}
    };

  EXPECT_THAT(wrap(mesh.generate(cutoff, true)), ::testing::UnorderedElementsAreArray(wrap(reference)));
}

TEST(TestMesh,TetrahedronMeshContainsSphere)
{
  auto cutoff = CutoffSphere(1.0);
  auto mesh = TetrahedronMesh(1.0);

  const static std::vector< Point3D > reference =
    {
      { 0.0, 0.0, 0.0},
      { 1.0, 0.0, 0.0}, { a, b, 0.0}, { a,-b, 0.0},
      {-1.0, 0.0, 0.0}, {-a, b, 0.0}, {-a,-b, 0.0},

      {0.0,-d, h}, { a, c, h}, {-a, c, h},
      {0.0, d,-h}, { a,-c,-h}, {-a,-c,-h}
    };

  EXPECT_THAT(wrap(mesh.generate(cutoff)), ::testing::UnorderedElementsAreArray(wrap(reference)));
}
