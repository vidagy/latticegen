#include <gtest/gtest.h>
#include <gmock/gmock.h>

#include <Core/Point2D.h>

using namespace Core::Geometry;

TEST(TestPoint2D,Ctor)
{
  const Point2D p = Point2D(1.0,2.0);
  EXPECT_DOUBLE_EQ(1.0, p.x );
  EXPECT_DOUBLE_EQ(2.0, p.y );
}

TEST(TestPoint2D,Swap)
{
  Point2D p = Point2D(1.0,2.0);
  Point2D q = Point2D(10.0,20.0);
  swap(p,q);
  EXPECT_DOUBLE_EQ(10.0, p.x );
  EXPECT_DOUBLE_EQ(20.0, p.y );
  EXPECT_DOUBLE_EQ(1.0, q.x );
  EXPECT_DOUBLE_EQ(2.0, q.y );
}

TEST(TestPoint2D,Copy)
{
  const Point2D p = Point2D(1.0,2.0);
  const Point2D q = p;
  EXPECT_TRUE( &q != &p);
  EXPECT_DOUBLE_EQ(p.x, q.x );
  EXPECT_DOUBLE_EQ(p.y, q.y );

  const Point2D r(p);
  EXPECT_TRUE( &r != &p);
  EXPECT_DOUBLE_EQ(p.x, r.x );
  EXPECT_DOUBLE_EQ(p.y, r.y );
}

TEST(TestPoint2D,Plus)
{
  const Point2D p = Point2D(1.0,2.0);
  const Point2D q = Point2D(10.0,20.0);
  
  const Point2D pq = p + q;
  EXPECT_DOUBLE_EQ(pq.x, p.x + q.x );
  EXPECT_DOUBLE_EQ(pq.y, p.y + q.y );
  
  const Point2D qp = q + p;
  EXPECT_DOUBLE_EQ(qp.x, p.x + q.x );
  EXPECT_DOUBLE_EQ(qp.y, p.y + q.y );

  Point2D r = Point2D(100.0,200.0);
  r += p;
  EXPECT_DOUBLE_EQ(r.x, p.x + 100.0 );
  EXPECT_DOUBLE_EQ(r.y, p.y + 200.0 );
  r += q;
  EXPECT_DOUBLE_EQ(r.x, p.x + q.x + 100.0 );
  EXPECT_DOUBLE_EQ(r.y, p.y + q.y + 200.0 );
}

TEST(TestPoint2D,Minus)
{
  const Point2D p = Point2D(1.0,2.0);
  const Point2D q = Point2D(10.0,20.0);
  
  const Point2D pq = p - q;
  EXPECT_DOUBLE_EQ(pq.x, p.x - q.x );
  EXPECT_DOUBLE_EQ(pq.y, p.y - q.y );
  
  const Point2D qp = q - p;
  EXPECT_DOUBLE_EQ(qp.x, q.x - p.x );
  EXPECT_DOUBLE_EQ(qp.y, q.y - p.y );

  Point2D r = Point2D(100.0,200.0);
  r -= p;
  EXPECT_DOUBLE_EQ(r.x, 100.0 - p.x );
  EXPECT_DOUBLE_EQ(r.y, 200.0 - p.y );
  r -= q;
  EXPECT_DOUBLE_EQ(r.x, 100.0 - p.x - q.x);
  EXPECT_DOUBLE_EQ(r.y, 200.0 - p.y - q.y);
}

TEST(TestPoint2D,MultiplicationWithDouble)
{
  const Point2D p = Point2D(1.0,2.0);

  const Point2D q = p * 2.0;
  EXPECT_DOUBLE_EQ(2.0, q.x );
  EXPECT_DOUBLE_EQ(4.0, q.y );

  const Point2D r = 2.0 * p;
  EXPECT_DOUBLE_EQ(2.0, r.x );
  EXPECT_DOUBLE_EQ(4.0, r.y );

  Point2D s = Point2D(1.0,2.0);
  s *= 3.0;
  EXPECT_DOUBLE_EQ(3.0, s.x );
  EXPECT_DOUBLE_EQ(6.0, s.y );
}

TEST(TestPoint2D,ScalarMultiplication)
{
  const Point2D p = Point2D(1.0,2.0);
  const Point2D q = Point2D(10.0,20.0);

  const double pq = p * q;
  EXPECT_DOUBLE_EQ(50.0, pq );
  
  const double qp = q * p;
  EXPECT_DOUBLE_EQ(50.0, qp );

}

TEST(TestPoint2D,Compare)
{
  const Point2D p = Point2D(1.0,2.0);
  const Point2D q = Point2D(1.0,2.0);
  const Point2D r = Point2D(2.0,2.0);
  const Point2D s = Point2D(1.0,3.0);
  const Point2D t = Point2D(2.0,3.0);

  EXPECT_TRUE(p == q);
  EXPECT_FALSE(p == r);
  EXPECT_FALSE(p == s);
  EXPECT_FALSE(p == t);
}

TEST(TestPoint2D,IsRectangular)
{
  const Point2D p = Point2D(1.0,2.0);
  const Point2D q = Point2D(1.0,2.0);
  
  EXPECT_FALSE(isRectangular(p,q));

  const Point2D r = Point2D(2.0,-1.0);

  EXPECT_TRUE(isRectangular(p,r));
  EXPECT_TRUE(isRectangular(r,p));
}
