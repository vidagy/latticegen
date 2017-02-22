#include <gtest/gtest.h>
#include <gmock/gmock.h>

#include <Core/Point3D.h>

using namespace Core;

TEST(TestPoint3D,Ctor)
{
  const Point3D p = Point3D(1.0, 2.0, 3.0);
  EXPECT_DOUBLE_EQ(1.0, p.x );
  EXPECT_DOUBLE_EQ(2.0, p.y );
  EXPECT_DOUBLE_EQ(3.0, p.z );
}

TEST(TestPoint3D,Swap)
{
  Point3D p = Point3D(1.0,2.0,3.0);
  Point3D q = Point3D(10.0,20.0,30.0);
  swap(p,q);
  EXPECT_DOUBLE_EQ(10.0, p.x );
  EXPECT_DOUBLE_EQ(20.0, p.y );
  EXPECT_DOUBLE_EQ(30.0, p.z );
  EXPECT_DOUBLE_EQ(1.0, q.x );
  EXPECT_DOUBLE_EQ(2.0, q.y );
  EXPECT_DOUBLE_EQ(3.0, q.z );
}

TEST(TestPoint3D,Copy)
{
  const Point3D p = Point3D(1.0,2.0,3.0);
  const Point3D q = p;
  EXPECT_TRUE( &q != &p);
  EXPECT_DOUBLE_EQ(p.x, q.x );
  EXPECT_DOUBLE_EQ(p.y, q.y );
  EXPECT_DOUBLE_EQ(p.z, q.z );

  const Point3D r(p);
  EXPECT_TRUE( &r != &p);
  EXPECT_DOUBLE_EQ(p.x, r.x );
  EXPECT_DOUBLE_EQ(p.y, r.y );
  EXPECT_DOUBLE_EQ(p.z, r.z );
}

TEST(TestPoint3D,Plus)
{
  const Point3D p = Point3D(1.0,2.0,3.0);
  const Point3D q = Point3D(10.0,20.0,30.0);
  
  const Point3D pq = p + q;
  EXPECT_DOUBLE_EQ(pq.x, p.x + q.x );
  EXPECT_DOUBLE_EQ(pq.y, p.y + q.y );
  EXPECT_DOUBLE_EQ(pq.z, p.z + q.z );
  
  const Point3D qp = q + p;
  EXPECT_DOUBLE_EQ(qp.x, p.x + q.x );
  EXPECT_DOUBLE_EQ(qp.y, p.y + q.y );
  EXPECT_DOUBLE_EQ(qp.z, p.z + q.z );

  Point3D r = Point3D(100.0,200.0,300.0);
  r += p;
  EXPECT_DOUBLE_EQ(r.x, p.x + 100.0 );
  EXPECT_DOUBLE_EQ(r.y, p.y + 200.0 );
  EXPECT_DOUBLE_EQ(r.z, p.z + 300.0 );
  r += q;
  EXPECT_DOUBLE_EQ(r.x, p.x + q.x + 100.0 );
  EXPECT_DOUBLE_EQ(r.y, p.y + q.y + 200.0 );
  EXPECT_DOUBLE_EQ(r.z, p.z + q.z + 300.0 );
}

TEST(TestPoint3D,Minus)
{
  const Point3D p = Point3D(1.0,2.0,3.0);
  const Point3D q = Point3D(10.0,20.0,30.0);
  
  const Point3D pq = p - q;
  EXPECT_DOUBLE_EQ(pq.x, p.x - q.x );
  EXPECT_DOUBLE_EQ(pq.y, p.y - q.y );
  EXPECT_DOUBLE_EQ(pq.z, p.z - q.z );
  
  const Point3D qp = q - p;
  EXPECT_DOUBLE_EQ(qp.x, q.x - p.x );
  EXPECT_DOUBLE_EQ(qp.y, q.y - p.y );
  EXPECT_DOUBLE_EQ(qp.z, q.z - p.z );

  Point3D r = Point3D(100.0,200.0,300.0);
  r -= p;
  EXPECT_DOUBLE_EQ(r.x, 100.0 - p.x );
  EXPECT_DOUBLE_EQ(r.y, 200.0 - p.y );
  EXPECT_DOUBLE_EQ(r.z, 300.0 - p.z );
  r -= q;
  EXPECT_DOUBLE_EQ(r.x, 100.0 - p.x - q.x);
  EXPECT_DOUBLE_EQ(r.y, 200.0 - p.y - q.y);
  EXPECT_DOUBLE_EQ(r.z, 300.0 - p.z - q.z);
}

TEST(TestPoint3D,MultiplicationWithDouble)
{
  const Point3D p = Point3D(1.0,2.0,3.0);

  const Point3D q = p * 2.0;
  EXPECT_DOUBLE_EQ(2.0, q.x );
  EXPECT_DOUBLE_EQ(4.0, q.y );
  EXPECT_DOUBLE_EQ(6.0, q.z );

  const Point3D r = 2.0 * p;
  EXPECT_DOUBLE_EQ(2.0, r.x );
  EXPECT_DOUBLE_EQ(4.0, r.y );
  EXPECT_DOUBLE_EQ(6.0, r.z );

  Point3D s = Point3D(1.0,2.0,3.0);
  s *= 3.0;
  EXPECT_DOUBLE_EQ(3.0, s.x );
  EXPECT_DOUBLE_EQ(6.0, s.y );
  EXPECT_DOUBLE_EQ(9.0, s.z );
}

TEST(TestPoint3D,ScalarProduct)
{
  const Point3D p = Point3D(1.0,2.0,3.0);
  const Point3D q = Point3D(10.0,20.0,30.0);

  const double pq = p * q;
  EXPECT_DOUBLE_EQ(140.0, pq );
  
  const double qp = q * p;
  EXPECT_DOUBLE_EQ(140.0, qp );

}

TEST(TestPoint3D,CrossProduct)
{
  const Point3D p = Point3D(1.0,0.0,0.0);
  const Point3D q = Point3D(0.0,1.0,0.0);

  const Point3D pxq = cross_product(p, q);
  EXPECT_EQ(Point3D(0.0, 0.0, 1.0), pxq );
  
  const Point3D qxp = cross_product(q, p);
  EXPECT_EQ(Point3D(0.0, 0.0, -1.0), qxp );

}

TEST(TestPoint3D,Compare)
{
  const Point3D p = Point3D( 1.0, 2.0, 3.0);
  const Point3D q = Point3D( 1.0, 2.0, 3.0);
  const Point3D r = Point3D(10.0, 2.0, 3.0);
  const Point3D s = Point3D( 1.0,20.0, 3.0);
  const Point3D t = Point3D( 1.0, 2.0,30.0);
  const Point3D u = Point3D(10.0,20.0, 3.0);
  const Point3D v = Point3D( 1.0,20.0,30.0);
  const Point3D w = Point3D(10.0, 2.0,30.0);
  const Point3D x = Point3D(10.0, 2.0,30.0);

  EXPECT_TRUE(p == q);
  EXPECT_FALSE(p == r);
  EXPECT_FALSE(p == s);
  EXPECT_FALSE(p == t);
  EXPECT_FALSE(p == u);
  EXPECT_FALSE(p == v);
  EXPECT_FALSE(p == w);
  EXPECT_FALSE(p == x);
}
