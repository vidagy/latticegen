#include <TestUtils/base.h>

#include <Geometry/Mesh.h>
#include <Geometry/SymmetryTransformationFactory.h>
#include <Geometry/IrreducibleWedge.h>

using namespace Geometry;

TEST(TestIrreducibleWedge,ReduceBySymmetriesCubic)
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

TEST(TestIrreducibleWedge,GetCubic)
{
  auto irreducible_points = IrreducibleWedge::get_irreducible_wedge(UnitCell3D::create_cubic_primitive(5.0), 5);

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