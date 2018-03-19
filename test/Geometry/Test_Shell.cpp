#include <TestUtils/base.h>

#include <Geometry/Shell.h>

using namespace Geometry;
using namespace Testing;

TEST(TestShell, SmallCube)
{
  auto unit_cell = UnitCell3D::create_cubic_primitive(1.0);
  auto mesh = LatticeMesh(unit_cell);
  auto cutoff = CutoffSphere(1.8);
  auto points = mesh.generate(cutoff);
  auto transformations = SymmetryTransformationFactory::generate(unit_cell);

  auto shells = Shell::get_shells(transformations, points);
  EXPECT_EQ(shells.size(), 4u);

  const auto reference_0 = std::vector<Point3D>{{0, 0, 0}};
  EXPECT_THAT(wrap(shells[0].points), ::testing::UnorderedElementsAreArray(wrap(reference_0)));

  const auto reference_1 = std::vector<Point3D>{
    {1,  0,  0},
    {0,  1,  0},
    {0,  0,  1},
    {-1, 0,  0},
    {0,  -1, 0},
    {0,  0,  -1}
  };
  EXPECT_THAT(wrap(shells[1].points), ::testing::UnorderedElementsAreArray(wrap(reference_1)));

  const auto reference_2 = std::vector<Point3D>{
    {1,  1,  0},
    {-1, 1,  0},
    {1,  -1, 0},
    {-1, -1, 0},
    {1,  0,  1},
    {-1, 0,  1},
    {1,  0,  -1},
    {-1, 0,  -1},
    {0,  1,  1},
    {0,  -1, 1},
    {0,  1,  -1},
    {0,  -1, -1}
  };
  EXPECT_THAT(wrap(shells[2].points), ::testing::UnorderedElementsAreArray(wrap(reference_2)));

  const auto reference_3 = std::vector<Point3D>{
    {1,  1,  1},
    {-1, 1,  1},
    {1,  -1, 1},
    {1,  1,  -1},
    {-1, -1, 1},
    {-1, 1,  -1},
    {1,  -1, -1},
    {-1, -1, -1}
  };
  EXPECT_THAT(wrap(shells[3].points), ::testing::UnorderedElementsAreArray(wrap(reference_3)));
}

TEST(TestShell, SameDistance)
{
  auto unit_cell = UnitCell3D::create_cubic_primitive(1.0);
  auto mesh = LatticeMesh(unit_cell);
  auto cutoff = CutoffSphere(3.05);
  auto points = mesh.generate(cutoff);
  auto transformations = SymmetryTransformationFactory::generate(unit_cell);

  auto skin = std::vector<Point3D>();
  std::copy_if(points.begin(), points.end(), std::back_inserter(skin),
               [](const Point3D &p)
               {
                 return p.norm() > 2.95;
               });

  auto shells = Shell::get_shells(transformations, skin);
  ASSERT_EQ(shells.size(), 2u);
  auto shell_0 = shells[0].points.size() == 24 ? shells[0] : shells[1];
  auto shell_1 = shells[0].points.size() == 24 ? shells[1] : shells[0];

  const auto reference_0 = std::vector<Point3D>{
    {1,  2,  2},
    {2,  1,  2},
    {2,  2,  1},
    {-1, 2,  2},
    {-2, 1,  2},
    {-2, 2,  1},
    {1,  -2, 2},
    {2,  -1, 2},
    {2,  -2, 1},
    {1,  2,  -2},
    {2,  1,  -2},
    {2,  2,  -1},
    {-1, -2, 2},
    {-2, -1, 2},
    {-2, -2, 1},
    {-1, 2,  -2},
    {-2, 1,  -2},
    {-2, 2,  -1},
    {1,  -2, -2},
    {2,  -1, -2},
    {2,  -2, -1},
    {-1, -2, -2},
    {-2, -1, -2},
    {-2, -2, -1}
  };
  EXPECT_THAT(wrap(shell_0.points), ::testing::UnorderedElementsAreArray(wrap(reference_0)));

  const auto reference_1 = std::vector<Point3D>{
    {3,  0,  0},
    {0,  3,  0},
    {0,  0,  3},
    {-3, 0,  0},
    {0,  -3, 0},
    {0,  0,  -3}
  };
  EXPECT_THAT(wrap(shell_1.points), ::testing::UnorderedElementsAreArray(wrap(reference_1)));

}
