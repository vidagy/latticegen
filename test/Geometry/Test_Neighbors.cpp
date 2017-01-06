#include <gtest/gtest.h>
#include <gmock/gmock.h>

#include <Geometry/Neighbors.h>

using namespace Core;
using namespace Geometry;

TEST(TestCutoffCube,CtorThrows)
{
  EXPECT_THROW(CutoffCube( 0.0), std::invalid_argument);
  EXPECT_THROW(CutoffCube(-1.5), std::invalid_argument);
}

TEST(TestCutoffCube,CutoffTrue)
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

TEST(TestCutoffCube,CutoffFalse)
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

TEST(TestCutoffSphere,CtorThrows)
{
  EXPECT_THROW(CutoffSphere( 0.0), std::invalid_argument);
  EXPECT_THROW(CutoffSphere(-1.5), std::invalid_argument);
}

TEST(TestCutoffSphere,CutoffTrue)
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

TEST(TestCutoffSphere,CutoffFalse)
{
  CutoffSphere cutoff = CutoffSphere(1.0);
  EXPECT_FALSE(cutoff.is_included( Point3D( 1.1, 0.0, 0.0) ) );
  EXPECT_FALSE(cutoff.is_included( Point3D( 0.0, 1.1, 0.0) ) );
  EXPECT_FALSE(cutoff.is_included( Point3D( 0.0, 0.0, 1.1) ) );
  EXPECT_FALSE(cutoff.is_included( Point3D(-1.1, 0.0, 0.0) ) );
  EXPECT_FALSE(cutoff.is_included( Point3D( 0.0,-1.1, 0.0) ) );
  EXPECT_FALSE(cutoff.is_included( Point3D( 0.0, 0.0,-1.1) ) );

  EXPECT_FALSE(cutoff.is_included( Point3D( 0.8, 0.8, 0.0) ) );
  EXPECT_FALSE(cutoff.is_included( Point3D(-0.8, 0.8, 0.0) ) );
  EXPECT_FALSE(cutoff.is_included( Point3D( 0.8,-0.8, 0.0) ) );
  EXPECT_FALSE(cutoff.is_included( Point3D(-0.8,-0.8, 0.0) ) );

  EXPECT_FALSE(cutoff.is_included( Point3D( 0.8, 0.0, 0.8) ) );
  EXPECT_FALSE(cutoff.is_included( Point3D(-0.8, 0.0, 0.8) ) );
  EXPECT_FALSE(cutoff.is_included( Point3D( 0.8, 0.0,-0.8) ) );
  EXPECT_FALSE(cutoff.is_included( Point3D(-0.8, 0.0,-0.8) ) );

  EXPECT_FALSE(cutoff.is_included( Point3D( 0.0, 0.8, 0.8) ) );
  EXPECT_FALSE(cutoff.is_included( Point3D( 0.0,-0.8, 0.8) ) );
  EXPECT_FALSE(cutoff.is_included( Point3D( 0.0, 0.8,-0.8) ) );
  EXPECT_FALSE(cutoff.is_included( Point3D( 0.0,-0.8,-0.8) ) );

  EXPECT_FALSE(cutoff.is_included( Point3D( 0.6, 0.6, 0.6) ) );
  EXPECT_FALSE(cutoff.is_included( Point3D(-0.6, 0.6, 0.6) ) );
  EXPECT_FALSE(cutoff.is_included( Point3D( 0.6,-0.6, 0.6) ) );
  EXPECT_FALSE(cutoff.is_included( Point3D( 0.6, 0.6,-0.6) ) );
  EXPECT_FALSE(cutoff.is_included( Point3D(-0.6,-0.6, 0.6) ) );
  EXPECT_FALSE(cutoff.is_included( Point3D(-0.6, 0.6,-0.6) ) );
  EXPECT_FALSE(cutoff.is_included( Point3D( 0.6,-0.6,-0.6) ) );
  EXPECT_FALSE(cutoff.is_included( Point3D(-0.6,-0.6,-0.6) ) );
}