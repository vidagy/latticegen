#include <gtest/gtest.h>
#include <gmock/gmock.h>

#include <Core/Point2D.h>

using namespace Core::Geometry;

TEST(TestPoint2D,Ctor)
{
  const Point2D p = Point2D(1.0,2.0);
  EXPECT_DOUBLE_EQ(1.0, p.x_ );
  EXPECT_DOUBLE_EQ(2.0, p.y_ );
}

TEST(TestPoint2D,Copy)
{
  const Point2D p = Point2D(1.0,2.0);
  const Point2D q = p;
  EXPECT_TRUE( &q != &p);
  EXPECT_DOUBLE_EQ(p.x_, q.x_ );
  EXPECT_DOUBLE_EQ(p.y_, q.y_ );

  const Point2D r(p);
  EXPECT_TRUE( &r != &p);
  EXPECT_DOUBLE_EQ(p.x_, r.x_ );
  EXPECT_DOUBLE_EQ(p.y_, r.y_ );
}

TEST(TestPoint2D,Plus)
{
  const Point2D p = Point2D(1.0,2.0);
  const Point2D q = Point2D(10.0,20.0);
  
  const Point2D pq = p + q;
  EXPECT_DOUBLE_EQ(pq.x_, p.x_ + q.x_ );
  EXPECT_DOUBLE_EQ(pq.y_, p.y_ + q.y_ );
  
  const Point2D qp = q + p;
  EXPECT_DOUBLE_EQ(qp.x_, p.x_ + q.x_ );
  EXPECT_DOUBLE_EQ(qp.y_, p.y_ + q.y_ );

  Point2D r = Point2D(100.0,200.0);
  r += p;
  EXPECT_DOUBLE_EQ(r.x_, p.x_ + 100.0 );
  EXPECT_DOUBLE_EQ(r.y_, p.y_ + 200.0 );
  r += q;
  EXPECT_DOUBLE_EQ(r.x_, p.x_ + q.x_ + 100.0 );
  EXPECT_DOUBLE_EQ(r.y_, p.y_ + q.y_ + 200.0 );
}

TEST(TestPoint2D,Minus)
{
  const Point2D p = Point2D(1.0,2.0);
  const Point2D q = Point2D(10.0,20.0);
  
  const Point2D pq = p - q;
  EXPECT_DOUBLE_EQ(pq.x_, p.x_ - q.x_ );
  EXPECT_DOUBLE_EQ(pq.y_, p.y_ - q.y_ );
  
  const Point2D qp = q - p;
  EXPECT_DOUBLE_EQ(qp.x_, q.x_ - p.x_ );
  EXPECT_DOUBLE_EQ(qp.y_, q.y_ - p.y_ );

  Point2D r = Point2D(100.0,200.0);
  r -= p;
  EXPECT_DOUBLE_EQ(r.x_, 100.0 - p.x_ );
  EXPECT_DOUBLE_EQ(r.y_, 200.0 - p.y_ );
  r -= q;
  EXPECT_DOUBLE_EQ(r.x_, 100.0 - p.x_ - q.x_);
  EXPECT_DOUBLE_EQ(r.y_, 200.0 - p.y_ - q.y_);
}

TEST(TestPoint2D,MultiplicationWithDouble)
{
  const Point2D p = Point2D(1.0,2.0);

  const Point2D q = p * 2.0;
  ASSERT_DOUBLE_EQ(2.0, q.x_ );
  ASSERT_DOUBLE_EQ(4.0, q.y_ );

  const Point2D r = 2.0 * p;
  ASSERT_DOUBLE_EQ(2.0, r.x_ );
  ASSERT_DOUBLE_EQ(4.0, r.y_ );

  Point2D s = Point2D(1.0,2.0);
  s *= 3.0;
  ASSERT_DOUBLE_EQ(3.0, s.x_ );
  ASSERT_DOUBLE_EQ(6.0, s.y_ );
}

TEST(TestPoint2D,ScalarMultiplication)
{
  const Point2D p = Point2D(1.0,2.0);
  const Point2D q = Point2D(10.0,20.0);

  const Point2D r = p * q;
  ASSERT_DOUBLE_EQ(10.0, r.x_ );
  ASSERT_DOUBLE_EQ(40.0, r.y_ );
  
  const Point2D s = q * p;
  ASSERT_DOUBLE_EQ(10.0, s.x_ );
  ASSERT_DOUBLE_EQ(40.0, s.y_ );

}
