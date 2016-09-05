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

  EXPECT_THROW(BravaisLattice2D::get_canonical_unit_cell(p00,p11), std::invalid_argument);
  EXPECT_THROW(BravaisLattice2D::get_canonical_unit_cell(p00,p01), std::invalid_argument);
  EXPECT_THROW(BravaisLattice2D::get_canonical_unit_cell(p00,p10), std::invalid_argument);

  EXPECT_THROW(BravaisLattice2D::get_canonical_unit_cell(p11,p00), std::invalid_argument);
  EXPECT_THROW(BravaisLattice2D::get_canonical_unit_cell(p01,p00), std::invalid_argument);
  EXPECT_THROW(BravaisLattice2D::get_canonical_unit_cell(p10,p00), std::invalid_argument);

  EXPECT_THROW(BravaisLattice2D::get_canonical_unit_cell(p11,p11), std::invalid_argument);
  EXPECT_THROW(BravaisLattice2D::get_canonical_unit_cell(p01,p01), std::invalid_argument);
  EXPECT_THROW(BravaisLattice2D::get_canonical_unit_cell(p10,p10), std::invalid_argument);
  EXPECT_THROW(BravaisLattice2D::get_canonical_unit_cell(p00,p00), std::invalid_argument);

  EXPECT_NO_THROW(BravaisLattice2D::get_canonical_unit_cell(p01,p10));
  EXPECT_NO_THROW(BravaisLattice2D::get_canonical_unit_cell(p01,p11));
  EXPECT_NO_THROW(BravaisLattice2D::get_canonical_unit_cell(p10,p01));
  EXPECT_NO_THROW(BravaisLattice2D::get_canonical_unit_cell(p10,p11));
  EXPECT_NO_THROW(BravaisLattice2D::get_canonical_unit_cell(p11,p01));
  EXPECT_NO_THROW(BravaisLattice2D::get_canonical_unit_cell(p11,p10));
}

TEST(BravaisLattice2D,GetCanonicalUnitCellAndScaleSimpleNoScale)
{
  const Point2D p01 = Point2D(0.0,1.0);
  const Point2D p10 = Point2D(1.0,0.0);
  Point2D resa, resb;

  std::tie(resa, resb) = BravaisLattice2D::get_canonical_unit_cell(p01,p10);

  EXPECT_EQ(resa, Point2D(1.0, 0.0));
  EXPECT_EQ(resb, Point2D(0.0, 1.0));

  std::tie(resa, resb) = BravaisLattice2D::get_canonical_unit_cell(p10,p01);

  EXPECT_EQ(resa, Point2D(1.0, 0.0));
  EXPECT_EQ(resb, Point2D(0.0, 1.0));
}

TEST(BravaisLattice2D,GetCanonicalUnitCellAndScaleSimpleScale)
{
  const Point2D p01 = Point2D(0.0,5.0);
  const Point2D p10 = Point2D(1.0,0.0);
  Point2D resa, resb;

  std::tie(resa, resb) = BravaisLattice2D::get_canonical_unit_cell(p01,p10);

  EXPECT_EQ(resa, Point2D(5.0,0.0));
  EXPECT_EQ(resb, Point2D(0.0,1.0));

  std::tie(resa, resb) = BravaisLattice2D::get_canonical_unit_cell(p10,p01);

  EXPECT_EQ(resa, Point2D(5.0,0.0));
  EXPECT_EQ(resb, Point2D(0.0,1.0));
}

TEST(BravaisLattice2D,GetCanonicalUnitCellAndScaleRotate)
{
  const Point2D p11  = Point2D(1.0,1.0);
  const Point2D pm11 = Point2D(-1.0,1.0);
  Point2D resa, resb;

  std::tie(resa, resb) = BravaisLattice2D::get_canonical_unit_cell(p11,pm11);

  EXPECT_EQ(resa, Point2D(sqrt(2.0), 0.0));
  EXPECT_EQ(resb, Point2D(0.0, sqrt(2.0)));

  std::tie(resa, resb) = BravaisLattice2D::get_canonical_unit_cell(pm11,p11);
  
  EXPECT_EQ(resa, Point2D(sqrt(2.0), 0.0));
  EXPECT_EQ(resb, Point2D(0.0, sqrt(2.0)));
}

TEST(BravaisLattice2D,GetCanonicalUnitCellAndScaleReflect)
{
  const Point2D p11  = Point2D(1.0,-2.0);
  const Point2D pm11 = Point2D(0.0,-1.0);
  Point2D resa, resb;

  std::tie(resa, resb) = BravaisLattice2D::get_canonical_unit_cell(p11,pm11);

  EXPECT_EQ(resa, Point2D(sqrt(5.0), 0.0));
  EXPECT_EQ(resb, Point2D(3.0/sqrt(5.0), 1.0/sqrt(5.0)));

  std::tie(resa, resb) = BravaisLattice2D::get_canonical_unit_cell(pm11,p11);

  EXPECT_EQ(resa, Point2D(sqrt(5.0), 0.0));
  EXPECT_EQ(resb, Point2D(3.0/sqrt(5.0), 1.0/sqrt(5.0)));
}

TEST(BravaisLattice2D,FindLatticeType)
{
  EXPECT_EQ(BravaisLattice2D::find_lattice_type(Point2D(1.0,0.0), Point2D(0.2,0.3)), BravaisLattice2D::Oblique);
  EXPECT_EQ(BravaisLattice2D::find_lattice_type(Point2D(1.0,0.0), Point2D(0.0,0.75)), BravaisLattice2D::Rectangular);
  EXPECT_EQ(BravaisLattice2D::find_lattice_type(Point2D(1.0,0.0), Point2D(0.5,0.25)), BravaisLattice2D::CenteredRectangular);
  EXPECT_EQ(BravaisLattice2D::find_lattice_type(Point2D(1.0,0.0), Point2D(0.5,sqrt(3)/2.0)), BravaisLattice2D::Hexagonal);
  EXPECT_EQ(BravaisLattice2D::find_lattice_type(Point2D(1.0,0.0), Point2D(0.0,1.0)), BravaisLattice2D::Square);
}

TEST(BravaisLattice2D,ObliqueLattice)
{
  const BravaisLattice2D lattice = BravaisLattice2D( Point2D(1.0,0.0), Point2D(0.2,0.3), 3, 3 );
  
  const static std::vector< Point2D > reference = 
  {
    {0.0,0.0},{1.0,0.0},{2.0,0.0},
    {0.2,0.3},{1.2,0.3},{2.2,0.3},
    {0.4,0.6},{1.4,0.6},{2.4,0.6}
  };
  
  EXPECT_THAT( lattice.get_lattice(), ::testing::ContainerEq(reference) );
}

TEST(BravaisLattice2D,RectangularLattice)
{
  const BravaisLattice2D lattice = BravaisLattice2D( Point2D(1.0,0.0), Point2D(0.0,0.5), 3, 3 );
  
  const static std::vector< Point2D > reference = 
  {
    {0.0,0.0},{1.0,0.0},{2.0,0.0},
    {0.0,0.5},{1.0,0.5},{2.0,0.5},
    {0.0,1.0},{1.0,1.0},{2.0,1.0}
  };
  
  EXPECT_THAT( lattice.get_lattice(), ::testing::ContainerEq(reference) );
}

TEST(BravaisLattice2D,CenteredRectangularLattice)
{
  const BravaisLattice2D lattice = BravaisLattice2D( Point2D(1.0,0.0), Point2D(0.5,0.25), 3, 6 );
  
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

TEST(BravaisLattice2D,HexagonalLattice)
{
  const BravaisLattice2D lattice = BravaisLattice2D( Point2D(1.0,0.0), Point2D(0.5,sqrt(3)/2.0), 3, 3 );
  
  const static std::vector< Point2D > reference = 
  {
    {0.0,0.0},{1.0,0.0},{2.0,0.0},
    {0.5,sqrt(3)/2.0},{1.5,sqrt(3)/2.0},{2.5,sqrt(3)/2.0},
    {1.0,sqrt(3)},{2.0,sqrt(3)},{3.0,sqrt(3)}
  };
  
  EXPECT_THAT( lattice.get_lattice(), ::testing::ContainerEq(reference) );
}

TEST(BravaisLattice2D,SquareLattice)
{
  const BravaisLattice2D lattice = BravaisLattice2D( Point2D(1.0,0.0), Point2D(0.0,1.0), 3, 3 );
  
  const static std::vector< Point2D > reference = 
  {
    {0.0,0.0},{1.0,0.0},{2.0,0.0},
    {0.0,1.0},{1.0,1.0},{2.0,1.0},
    {0.0,2.0},{1.0,2.0},{2.0,2.0}
  };
  
  EXPECT_THAT( lattice.get_lattice(), ::testing::ContainerEq(reference) );
}

TEST(BravaisLattice2D,WedgeOblique)
{
  Point2D a = Point2D(1.0,0.0);
  Point2D b = Point2D(0.4,0.6);
  const BravaisLattice2D::BravaisLattice2DType latticeType = BravaisLattice2D::find_lattice_type(a, b);
  EXPECT_THAT(latticeType, BravaisLattice2D::Oblique);
  
  const static std::vector< Point2D > reference = 
  {
    {0.1 , 0.06 },
    {0.3 , 0.06 }, {0.3 , 0.26 },
    {0.5 , 0.06 }, {0.5 , 0.26 }, {0.5 , 0.46 },
    {0.7 , 0.06 }, {0.7 , 0.26 },
    {0.9 , 0.06 }
  };

  //std::copy(res.begin(), res.end(), std::ostream_iterator<Point2D>(std::cout, "\n"));
  EXPECT_THAT( BravaisLattice2D::get_irreducible_wedge(a, b, 5, 5), ::testing::ContainerEq(reference) );
}

