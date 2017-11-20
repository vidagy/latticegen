#include <TestUtils/base.h>

#include <Core/Point3D.h>
#include <Math/CommonFunctions.h>

using namespace Core;

TEST(TestPoint3D,Ctor)
{
  const auto p = Point3D(1.0, 2.0, 3.0);
  EXPECT_DOUBLE_EQ(1.0, p(0));
  EXPECT_DOUBLE_EQ(2.0, p(1));
  EXPECT_DOUBLE_EQ(3.0, p(2));
}

TEST(TestPoint3D, CreatPolar)
{
  auto xp = create_polar(1.0, pi / 2.0, 0.0);
  auto x = Point3D(1.0, 0.0, 0.0);
  EXPECT_TRUE(default_point_comparator.isEqual(xp, x));

  auto yp = create_polar(1.0, pi / 2.0, pi / 2.0);
  auto y = Point3D(0.0, 1.0, 0.0);
  EXPECT_TRUE(default_point_comparator.isEqual(yp, y));

  auto zp = create_polar(1.0, 0.0, 0.0);
  auto z = Point3D(0.0, 0.0, 1.0);
  EXPECT_TRUE(default_point_comparator.isEqual(zp, z));
}

TEST(TestPoint3D,Copy)
{
  const auto p = Point3D(1.0, 2.0, 3.0);
  const auto q = p;
  EXPECT_TRUE( &q != &p);
  EXPECT_DOUBLE_EQ(p(0), q(0));
  EXPECT_DOUBLE_EQ(p(1), q(1));
  EXPECT_DOUBLE_EQ(p(2), q(2));

  const Point3D r(p);
  EXPECT_TRUE( &r != &p);
  EXPECT_DOUBLE_EQ(p(0), r(0));
  EXPECT_DOUBLE_EQ(p(1), r(1));
  EXPECT_DOUBLE_EQ(p(2), r(2));
}

TEST(TestPoint3D,Plus)
{
  const Point3D p = Point3D(1.0,2.0,3.0);
  const Point3D q = Point3D(10.0,20.0,30.0);
  
  const Point3D pq = p + q;
  EXPECT_DOUBLE_EQ(pq(0), p(0) + q(0));
  EXPECT_DOUBLE_EQ(pq(1), p(1) + q(1));
  EXPECT_DOUBLE_EQ(pq(2), p(2) + q(2));
  
  const Point3D qp = q + p;
  EXPECT_DOUBLE_EQ(qp(0), p(0) + q(0));
  EXPECT_DOUBLE_EQ(qp(1), p(1) + q(1));
  EXPECT_DOUBLE_EQ(qp(2), p(2) + q(2));

  Point3D r = Point3D(100.0,200.0,300.0);
  r += p;
  EXPECT_DOUBLE_EQ(r(0), p(0) + 100.0);
  EXPECT_DOUBLE_EQ(r(1), p(1) + 200.0);
  EXPECT_DOUBLE_EQ(r(2), p(2) + 300.0);
  r += q;
  EXPECT_DOUBLE_EQ(r(0), p(0) + q(0) + 100.0);
  EXPECT_DOUBLE_EQ(r(1), p(1) + q(1) + 200.0);
  EXPECT_DOUBLE_EQ(r(2), p(2) + q(2) + 300.0);
}

TEST(TestPoint3D,Minus)
{
  const auto p = Point3D(1.0, 2.0, 3.0);
  const auto q = Point3D(10.0, 20.0, 30.0);

  const auto pq = p - q;
  EXPECT_DOUBLE_EQ(pq(0), p(0) - q(0));
  EXPECT_DOUBLE_EQ(pq(1), p(1) - q(1));
  EXPECT_DOUBLE_EQ(pq(2), p(2) - q(2));

  const auto qp = q - p;
  EXPECT_DOUBLE_EQ(qp(0), q(0) - p(0));
  EXPECT_DOUBLE_EQ(qp(1), q(1) - p(1));
  EXPECT_DOUBLE_EQ(qp(2), q(2) - p(2));

  Point3D r = Point3D(100.0,200.0,300.0);
  r -= p;
  EXPECT_DOUBLE_EQ(r(0), 100.0 - p(0));
  EXPECT_DOUBLE_EQ(r(1), 200.0 - p(1));
  EXPECT_DOUBLE_EQ(r(2), 300.0 - p(2));
  r -= q;
  EXPECT_DOUBLE_EQ(r(0), 100.0 - p(0) - q(0));
  EXPECT_DOUBLE_EQ(r(1), 200.0 - p(1) - q(1));
  EXPECT_DOUBLE_EQ(r(2), 300.0 - p(2) - q(2));
}

TEST(TestPoint3D,MultiplicationWithDouble)
{
  const auto p = Point3D(1.0, 2.0, 3.0);

  const auto q = p * 2.0;
  EXPECT_DOUBLE_EQ(2.0, q(0));
  EXPECT_DOUBLE_EQ(4.0, q(1));
  EXPECT_DOUBLE_EQ(6.0, q(2));

  const auto r = 2.0 * p;
  EXPECT_DOUBLE_EQ(2.0, r(0));
  EXPECT_DOUBLE_EQ(4.0, r(1));
  EXPECT_DOUBLE_EQ(6.0, r(2));

  auto s = Point3D(1.0, 2.0, 3.0);
  s *= 3.0;
  EXPECT_DOUBLE_EQ(3.0, s(0));
  EXPECT_DOUBLE_EQ(6.0, s(1));
  EXPECT_DOUBLE_EQ(9.0, s(2));
}

TEST(TestPoint3D,ScalarProduct)
{
  const auto p = Point3D(1.0, 2.0, 3.0);
  const auto q = Point3D(10.0, 20.0, 30.0);

  const double pq = p.dot(q);
  EXPECT_DOUBLE_EQ(140.0, pq );

  const double qp = q.dot(p);
  EXPECT_DOUBLE_EQ(140.0, qp );

}

TEST(TestPoint3D,CrossProduct)
{
  const auto p = Point3D(1.0, 0.0, 0.0);
  const auto q = Point3D(0.0, 1.0, 0.0);

  const auto pxq = p.cross(q);
  EXPECT_EQ(Point3D(0.0, 0.0, 1.0), pxq );

  const auto qxp = q.cross(p);
  EXPECT_EQ(Point3D(0.0, 0.0, -1.0), qxp );

}

TEST(TestPoint3D,Compare)
{
  const auto p = Point3D(1.0, 2.0, 3.0);
  const auto q = Point3D(1.0, 2.0, 3.0);
  const auto r = Point3D(10.0, 2.0, 3.0);
  const auto s = Point3D(1.0, 20.0, 3.0);
  const auto t = Point3D(1.0, 2.0, 30.0);
  const auto u = Point3D(10.0, 20.0, 3.0);
  const auto v = Point3D(1.0, 20.0, 30.0);
  const auto w = Point3D(10.0, 2.0, 30.0);
  const auto x = Point3D(10.0, 2.0, 30.0);

  EXPECT_TRUE(p == q);
  EXPECT_FALSE(p == r);
  EXPECT_FALSE(p == s);
  EXPECT_FALSE(p == t);
  EXPECT_FALSE(p == u);
  EXPECT_FALSE(p == v);
  EXPECT_FALSE(p == w);
  EXPECT_FALSE(p == x);
}

TEST(TestPoint3D, CompareWithTolerance)
{
  const auto p = Point3D(1.0, 2.0, 3.0);
  const auto q = Point3D(1.01, 2.0, 3.0);

  auto comparator1 = Point3DComparator(0.1, 1.0);
  EXPECT_TRUE(comparator1.isEqual(p, q));
  auto comparator2 = Point3DComparator(0.001, 1.0);
  EXPECT_FALSE(comparator2.isEqual(p, q));

  const auto r = Point3D(100.0, 200.0, 300.0);
  const auto s = Point3D(100.01, 200.0, 300.0);

  auto comparator3 = Point3DComparator(0.1, 1e-3);
  EXPECT_TRUE(comparator3.isEqual(r, s));
  auto comparator4 = Point3DComparator(0.1, 1e-7);
  EXPECT_FALSE(comparator4.isEqual(r, s));
}
