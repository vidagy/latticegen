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

namespace
{
  struct NoData
  {
    bool operator==(const NoData &other) const { return true; }
  };
}

TEST(TestIrreducibleWedge, Replicate)
{
  auto cell = UnitCell3D::create_cubic_primitive(2.0);
  auto group = CrystallographicPointGroup::create(cell.get_point_group());
  auto transformations = SymmetryTransformationFactory::get(group->get_elements());

  auto irreducible_points = IrreducibleWedge::get_irreducible_wedge(cell, 2);

  auto irreduc = std::vector<std::pair<Point3D, NoData>>();
  std::transform(irreducible_points.begin(), irreducible_points.end(), std::back_inserter(irreduc),
                 [](const Point3D &x)
                 {
                   return std::make_pair(x, NoData());
                 });

  auto replicated = IrreducibleWedge::replicate(irreduc, transformations);
  auto replicated_points = std::vector<Point3D>();
  std::transform(replicated.begin(), replicated.end(), std::back_inserter(replicated_points),
                 [](const std::pair<Point3D, NoData> &x)
                 {
                   return x.first;
                 });

  auto reference_mesh = CubicMesh(1.0).generate(CutoffCube(1.0));

  EXPECT_THAT(replicated_points, ::testing::UnorderedElementsAreArray(reference_mesh));
}
