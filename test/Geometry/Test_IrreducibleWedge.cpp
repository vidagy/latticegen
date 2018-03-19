#include <TestUtils/base.h>

#include <Geometry/Mesh.h>
#include <Geometry/SymmetryTransformationFactory.h>
#include <Geometry/IrreducibleWedge.h>

using namespace Geometry;
using namespace Testing;

TEST(TestIrreducibleWedge,ReduceBySymmetriesCubic)
{
  auto mesh = CubicMesh(1.0);
  auto cutoff = CutoffCube(2.0);

  auto points = mesh.generate(cutoff);
  auto symmetries = SymmetryTransformationFactory::generate(Oh().get_generators());
  auto irreducible_points = IrreducibleWedge::reduce_by_symmetries(points, symmetries);

  const std::vector<Point3D> reference {
    {0, 0, 0},
    {0, 0, -1},
    {-1, -1, 0},
    {-1, 1, -1},
    {2, 0, 0},
    {-2, -1, 0},
    {-1, -2, 1},
    {0, -2, 2},
    {-2, 1, 2},
    {2, 2, 2}
  };
  EXPECT_THAT(wrap(irreducible_points), ::testing::UnorderedElementsAreArray(wrap(reference)));
}

TEST(TestIrreducibleWedge,GetCubic)
{
  auto irreducible_points = IrreducibleWedge::get_irreducible_wedge(UnitCell3D::create_cubic_primitive(5.0), 5);

  const std::vector<Point3D> reference {
    {0,  0,  0},
    {0,  0,  -1},
    {-1, -1, 0},
    {-1, 1,  -1},
    {2,  0,  0},
    {-2, -1, 0},
    {-1, -2, 1},
    {0,  -2, 2},
    {-2, 1,  2},
    {2,  2,  2}
  };
  EXPECT_THAT(wrap(irreducible_points), ::testing::UnorderedElementsAreArray(wrap(reference)));
}

namespace
{
  struct NoData
  {
    LATTICEGEN_MUTE_BEGIN
    LATTICEGEN_MUTE_UNUSED_VAR

    bool operator==(const NoData &other) const {
      LATTICEGEN_MUTE_END
      return true;
    }
  };

  std::vector<std::pair<Point3D, NoData>> zip(const std::vector<Point3D> &points)
  {
    auto res = std::vector<std::pair<Point3D, NoData>>();
    res.reserve(points.size());
    std::transform(points.begin(), points.end(), std::back_inserter(res),
                   [](const Point3D &x)
                   {
                     return std::make_pair(x, NoData());
                   });
    return res;
  }

  std::vector<Point3D> unzip(const std::vector<std::pair<Point3D, NoData>> &points_with_nodata)
  {
    auto res = std::vector<Point3D>();
    std::transform(points_with_nodata.begin(), points_with_nodata.end(), std::back_inserter(res),
                   [](const std::pair<Point3D, NoData> &x)
                   {
                     return x.first;
                   });
    return res;
  }
}

TEST(TestIrreducibleWedge, Replicate)
{
  auto cell = UnitCell3D::create_cubic_primitive(2.0);
  auto transformations = SymmetryTransformationFactory::generate(cell);

  auto irreducible_points = IrreducibleWedge::get_irreducible_wedge(cell, 2);

  auto irreduc = zip(irreducible_points);
  auto replicated = IrreducibleWedge::replicate(irreduc, transformations);
  auto replicated_points = unzip(replicated);

  auto reference_mesh = CubicMesh(1.0).generate(CutoffCube(1.0));

  EXPECT_THAT(wrap(replicated_points), ::testing::UnorderedElementsAreArray(wrap(reference_mesh)));
}

TEST(TestIrreducibleWedge, Roundtrip)
{
  auto a = 1.0;
  auto n = 5;
  auto abs_tol = IrreducibleWedge::get_tolerance(UnitCell3D::create_cubic_primitive(a));
  auto mesh = CubicMesh(a);
  auto cutoff = CutoffCube(n * a);
  auto all_points = mesh.generate(cutoff);
  auto points = std::vector<Point3D>();
  std::copy_if(all_points.begin(), all_points.end(), std::back_inserter(points),
               [](const Point3D &p)
               {
                 return p.norm() > sqrt(3.0) * 4.42 && p.norm() < sqrt(3.0) * 5.1;
               });
  auto transformations = SymmetryTransformationFactory::generate(Oh().get_generators());
  auto irreducible_points = IrreducibleWedge::reduce_by_symmetries(points, transformations, abs_tol);

  auto irreduc = zip(irreducible_points);
  auto replicated = IrreducibleWedge::replicate(irreduc, transformations, abs_tol);
  auto replicated_points = unzip(replicated);

  EXPECT_THAT(wrap(replicated_points), ::testing::UnorderedElementsAreArray(wrap(points)));
}


TEST(DISABLED_TestIrreducibleWedge, Speed)
{
  auto mesh = CubicMesh(1.0);
  auto cutoff = CutoffWSCell(UnitCell3D::create_cubic_primitive(10.0));

  std::cout << "generating mesh..." << std::endl;
  auto points = mesh.generate(cutoff);
  std::cout << "generating mesh... Done" << std::endl;
  std::cout << "mesh size = " << points.size() << std::endl;

  std::cout << "generating symmetries..." << std::endl;
  auto symmetries = SymmetryTransformationFactory::generate(Oh().get_generators());
  std::cout << "generating symmetries... Done" << std::endl;

  std::cout << "reduce_by_symmetries..." << std::endl;
  auto irreducible_points = IrreducibleWedge::reduce_by_symmetries(points, symmetries);
  std::cout << "reduce_by_symmetries...Done" << std::endl;

  std::cout << irreducible_points.size() << std::endl;
}
