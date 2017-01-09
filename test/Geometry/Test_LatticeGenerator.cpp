#include <gtest/gtest.h>
#include <gmock/gmock.h>

#include <Geometry/LatticeGenerator.h>

using namespace Core;
using namespace Geometry;

namespace
{
    static const double pi = 3.14159265358979323846;
}

TEST(TestLatticeGenerator,LatticeGeneratorContainsOrigin)
{
  UnitCell3D unit_cell = UnitCell3D::create_cubic_primitive(1.0);
  std::shared_ptr<Cutoff> cutoff = std::make_shared<CutoffCube>(unit_cell, 0.1);
  LatticeGenerator latticeGenerator(cutoff);

  const static std::vector< Point3D > reference = { {0.0,0.0,0.0} };

  EXPECT_THAT( latticeGenerator.generate(), ::testing::ContainerEq(reference) );
  EXPECT_THAT( latticeGenerator.generate(true), ::testing::ContainerEq(reference) );
}

TEST(TestLatticeGenerator,LatticeGeneratorContainsCube)
{
  UnitCell3D unit_cell = UnitCell3D::create_cubic_primitive(1.0);
  std::shared_ptr<Cutoff> cutoff = std::make_shared<CutoffCube>(unit_cell, 1.0);
  LatticeGenerator latticeGenerator(cutoff);

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

  EXPECT_THAT( latticeGenerator.generate(), ::testing::ContainerEq(reference) );
}

TEST(TestLatticeGenerator,LatticeGeneratorContainsCubePositive)
{
  UnitCell3D unit_cell = UnitCell3D::create_cubic_primitive(1.0);
  std::shared_ptr<Cutoff> cutoff = std::make_shared<CutoffCube>(unit_cell, 1.0);
  LatticeGenerator latticeGenerator(cutoff);

  const static std::vector< Point3D > reference =
  {
    { 0.0, 0.0, 0.0},{ 1.0, 0.0, 0.0},
    { 0.0, 1.0, 0.0},{ 1.0, 1.0, 0.0},
    { 0.0, 0.0, 1.0},{ 1.0, 0.0, 1.0},
    { 0.0, 1.0, 1.0},{ 1.0, 1.0, 1.0}
  };

  EXPECT_THAT( latticeGenerator.generate(true), ::testing::ContainerEq(reference) );
}
