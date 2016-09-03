#include <gtest/gtest.h>
#include <gmock/gmock.h>

#include <Core/BravaisLattice2D.h>

#include <stdexcept>
#include <tuple>
#include <math.h>

using namespace Core::Geometry;

TEST(BravaisLattice2D,GetCanonicalUnitCellAndScaleThrows)
{
  const Point2D p00 = Point2D(0.0,0.0);
  const Point2D p01 = Point2D(0.0,1.0);
  const Point2D p10 = Point2D(1.0,0.0);
  const Point2D p11 = Point2D(1.0,1.0);

  EXPECT_THROW(BravaisLattice2D::get_canonical_unit_cell_and_scale(p00,p11), std::invalid_argument);
  EXPECT_THROW(BravaisLattice2D::get_canonical_unit_cell_and_scale(p00,p01), std::invalid_argument);
  EXPECT_THROW(BravaisLattice2D::get_canonical_unit_cell_and_scale(p00,p10), std::invalid_argument);

  EXPECT_THROW(BravaisLattice2D::get_canonical_unit_cell_and_scale(p11,p00), std::invalid_argument);
  EXPECT_THROW(BravaisLattice2D::get_canonical_unit_cell_and_scale(p01,p00), std::invalid_argument);
  EXPECT_THROW(BravaisLattice2D::get_canonical_unit_cell_and_scale(p10,p00), std::invalid_argument);

  EXPECT_THROW(BravaisLattice2D::get_canonical_unit_cell_and_scale(p11,p11), std::invalid_argument);
  EXPECT_THROW(BravaisLattice2D::get_canonical_unit_cell_and_scale(p01,p01), std::invalid_argument);
  EXPECT_THROW(BravaisLattice2D::get_canonical_unit_cell_and_scale(p10,p10), std::invalid_argument);
  EXPECT_THROW(BravaisLattice2D::get_canonical_unit_cell_and_scale(p00,p00), std::invalid_argument);

  EXPECT_NO_THROW(BravaisLattice2D::get_canonical_unit_cell_and_scale(p01,p10));
  EXPECT_NO_THROW(BravaisLattice2D::get_canonical_unit_cell_and_scale(p01,p11));
  EXPECT_NO_THROW(BravaisLattice2D::get_canonical_unit_cell_and_scale(p10,p01));
  EXPECT_NO_THROW(BravaisLattice2D::get_canonical_unit_cell_and_scale(p10,p11));
  EXPECT_NO_THROW(BravaisLattice2D::get_canonical_unit_cell_and_scale(p11,p01));
  EXPECT_NO_THROW(BravaisLattice2D::get_canonical_unit_cell_and_scale(p11,p10));
}

TEST(BravaisLattice2D,GetCanonicalUnitCellAndScaleSimpleNoScale)
{
  const Point2D p01 = Point2D(0.0,1.0);
  const Point2D p10 = Point2D(1.0,0.0);
  Point2D res;
  double scale;

  std::tie(res, scale) = BravaisLattice2D::get_canonical_unit_cell_and_scale(p01,p10);

  EXPECT_EQ(res, Point2D(0.0, 1.0));
  EXPECT_DOUBLE_EQ(scale, 1.0);

  std::tie(res, scale) = BravaisLattice2D::get_canonical_unit_cell_and_scale(p10,p01);

  EXPECT_EQ(res, Point2D(0.0, 1.0));
  EXPECT_DOUBLE_EQ(scale, 1.0);
}

TEST(BravaisLattice2D,GetCanonicalUnitCellAndScaleSimpleScale)
{
  const Point2D p01 = Point2D(0.0,5.0);
  const Point2D p10 = Point2D(1.0,0.0);
  Point2D res;
  double scale;

  std::tie(res, scale) = BravaisLattice2D::get_canonical_unit_cell_and_scale(p01,p10);

  EXPECT_EQ(res, Point2D(0.0, 1.0 / 5.0));
  EXPECT_DOUBLE_EQ(scale, 5.0);

  std::tie(res, scale) = BravaisLattice2D::get_canonical_unit_cell_and_scale(p10,p01);

  EXPECT_EQ(res, Point2D(0.0, 1.0 / 5.0));
  EXPECT_DOUBLE_EQ(scale, 5.0);
}

TEST(BravaisLattice2D,GetCanonicalUnitCellAndScaleRotate)
{
  const Point2D p11  = Point2D(1.0,1.0);
  const Point2D pm11 = Point2D(-1.0,1.0);
  Point2D res;
  double scale;

  std::tie(res, scale) = BravaisLattice2D::get_canonical_unit_cell_and_scale(p11,pm11);

  EXPECT_EQ(res, Point2D(0.0, 1.0));
  EXPECT_DOUBLE_EQ(scale, sqrt(2.0));

  std::tie(res, scale) = BravaisLattice2D::get_canonical_unit_cell_and_scale(pm11,p11);

  EXPECT_EQ(res, Point2D(0.0, 1.0));
  EXPECT_DOUBLE_EQ(scale, sqrt(2.0));
}

TEST(BravaisLattice2D,GetCanonicalUnitCellAndScaleReflect)
{
  const Point2D p11  = Point2D(1.0,-2.0);
  const Point2D pm11 = Point2D(0.0,-1.0);
  Point2D res;
  double scale;

  std::tie(res, scale) = BravaisLattice2D::get_canonical_unit_cell_and_scale(p11,pm11);

  EXPECT_EQ(res, Point2D(3.0/5.0, 1.0/5.0));
  EXPECT_DOUBLE_EQ(scale, sqrt(5.0));

  std::tie(res, scale) = BravaisLattice2D::get_canonical_unit_cell_and_scale(pm11,p11);

  EXPECT_EQ(res, Point2D(3.0/5.0, 1.0/5.0));
  EXPECT_DOUBLE_EQ(scale, sqrt(5.0));
}

TEST(BravaisLattice2D,FindLatticeType)
{
  EXPECT_EQ(BravaisLattice2D::find_lattice_type(Point2D(0.2,0.3)), BravaisLattice2D::Oblique);
  EXPECT_EQ(BravaisLattice2D::find_lattice_type(Point2D(0.0,0.75)), BravaisLattice2D::Rectangular);
  EXPECT_EQ(BravaisLattice2D::find_lattice_type(Point2D(0.5,0.25)), BravaisLattice2D::CenteredRectangular);
  EXPECT_EQ(BravaisLattice2D::find_lattice_type(Point2D(0.5,sqrt(3)/2.0)), BravaisLattice2D::Hexagonal);
  EXPECT_EQ(BravaisLattice2D::find_lattice_type(Point2D(0.0,1.0)), BravaisLattice2D::Square);
}

TEST(BravaisLattice2D,CtorOblique)
{
  const BravaisLattice2D lattice = BravaisLattice2D( Point2D(0.2,0.3), 1.0, 3, 3 );
  
  const static std::vector< Point2D > reference = 
  {
    {0.0,0.0},{1.0,0.0},{2.0,0.0},
    {0.2,0.3},{1.2,0.3},{2.2,0.3},
    {0.4,0.6},{1.4,0.6},{2.4,0.6}
  };
  
  EXPECT_THAT( lattice.get_lattice(), ::testing::ContainerEq(reference) );
}

TEST(BravaisLattice2D,CtorRectangular)
{
  const BravaisLattice2D lattice = BravaisLattice2D( Point2D(0.0,0.5), 1.0, 3, 3 );
  
  const static std::vector< Point2D > reference = 
  {
    {0.0,0.0},{1.0,0.0},{2.0,0.0},
    {0.0,0.5},{1.0,0.5},{2.0,0.5},
    {0.0,1.0},{1.0,1.0},{2.0,1.0}
  };
  
  EXPECT_THAT( lattice.get_lattice(), ::testing::ContainerEq(reference) );
}

TEST(BravaisLattice2D,CtorCenteredRectangular)
{
  const BravaisLattice2D lattice = BravaisLattice2D( Point2D(0.5,0.25), 1.0, 3, 6 );
  
  const static std::vector< Point2D > reference = 
  {
    {0.0,0.00},{1.0,0.00},{2.0,0.00},
    {0.5,0.25},{1.5,0.25},{2.5,0.25},
    {1.0,0.50},{2.0,0.50},{3.0,0.50},
    {1.5,0.75},{2.5,0.75},{3.5,0.75},
    {2.0,1.00},{3.0,1.00},{4.0,1.00},
    {2.5,1.25},{3.5,1.25},{4.5,1.25}
  };
  
  EXPECT_THAT( lattice.get_lattice(), ::testing::ContainerEq(reference) );
}

TEST(BravaisLattice2D,CtorHexagonal)
{
  const BravaisLattice2D lattice = BravaisLattice2D( Point2D(0.5,sqrt(3)/2.0), 1.0, 3, 3 );
  
  const static std::vector< Point2D > reference = 
  {
    {0.0,0.0},{1.0,0.0},{2.0,0.0},
    {0.5,sqrt(3)/2.0},{1.5,sqrt(3)/2.0},{2.5,sqrt(3)/2.0},
    {1.0,sqrt(3)},{2.0,sqrt(3)},{3.0,sqrt(3)}
  };
  
  EXPECT_THAT( lattice.get_lattice(), ::testing::ContainerEq(reference) );
}

TEST(BravaisLattice2D,CtorSquare)
{
  const BravaisLattice2D lattice = BravaisLattice2D( Point2D(0.0,1.0), 1.0, 3, 3 );
  
  const static std::vector< Point2D > reference = 
  {
    {0.0,0.0},{1.0,0.0},{2.0,0.0},
    {0.0,1.0},{1.0,1.0},{2.0,1.0},
    {0.0,2.0},{1.0,2.0},{2.0,2.0}
  };
  
  EXPECT_THAT( lattice.get_lattice(), ::testing::ContainerEq(reference) );
}