#include <gtest/gtest.h>
#include <gmock/gmock.h>

#include <Geometry/Mesh.h>
#include <Geometry/SymmetryTransformationFactory.h>
#include <Geometry/IrreducibleWedge.h>

using namespace Geometry;

TEST(TestIrreducibleWedge,Cubic)
{
  CubicMesh mesh = CubicMesh(1.0);
  CutoffCube cutoff = CutoffCube(2.0);

  auto points = mesh.generate(cutoff);
  auto symmetries = SymmetryTransformationFactory::generate(Oh().get_generators());
  auto irreducible_points = IrreducibleWedge::reduce_by_symmetries(points, symmetries);

  const std::vector<Point3D> reference {
    {-2, -2, -2},
    {-1, -2, -2},
    { 0, -2, -2},
    {-1, -1, -2},
    { 0, -1, -2},
    { 0,  0, -2},
    {-1, -1, -1},
    { 0, -1, -1},
    { 0,  0, -1},
    { 0,  0,  0}
  };
  EXPECT_THAT( irreducible_points, ::testing::UnorderedElementsAreArray(reference) );
}