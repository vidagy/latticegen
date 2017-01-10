#include <gtest/gtest.h>
#include <gmock/gmock.h>

#include <Geometry/SymmetryOperations.h>

namespace
{
  static const double pi = 3.14159265358979323846;
}

using namespace Geometry;

TEST(TestRotation,CtorThrows)
{
  EXPECT_THROW(Rotation(Vector3D(0.0, 0.0, 0.0)), std::invalid_argument);
}

TEST(TestRotation,Identity)
{
  const Vector3D v = 2.0 * pi / sqrt(1.0 + 4.0 + 9.0) * Vector3D(1.0, 2.0, 3.0);
  Rotation rotation = Rotation(v);

  const Vector3D v100 = Vector3D(1.0,0.0,0.0);
  const Vector3D v010 = Vector3D(0.0,1.0,0.0);
  const Vector3D v001 = Vector3D(0.0,0.0,1.0);

  EXPECT_EQ(v100, rotation(v100));
  EXPECT_EQ(v010, rotation(v010));
  EXPECT_EQ(v001, rotation(v001));
}

TEST(TestRotation,Degrees90)
{
  const Vector3D x = pi / 2.0 * Vector3D(1.0, 0.0, 0.0);
  const Vector3D y = pi / 2.0 * Vector3D(0.0, 1.0, 0.0);
  const Vector3D z = pi / 2.0 * Vector3D(0.0, 0.0, 1.0);

  Rotation rotation_x = Rotation(x);
  Rotation rotation_y = Rotation(y);
  Rotation rotation_z = Rotation(z);

  // rotation with x
  EXPECT_EQ(x, rotation_x(x));
  EXPECT_EQ(z, rotation_x(y));
  EXPECT_EQ(-1.0 * y, rotation_x(z));

  // rotation with y
  EXPECT_EQ(-1.0 * z, rotation_y(x));
  EXPECT_EQ(y, rotation_y(y));
  EXPECT_EQ(x, rotation_y(z));

  // rotation with z
  EXPECT_EQ(y, rotation_z(x));
  EXPECT_EQ(-1.0 * x, rotation_z(y));
  EXPECT_EQ(z, rotation_z(z));
}

TEST(TestReflection,CtorThrows)
{
  EXPECT_THROW(Reflection(Vector3D(0.0, 0.0, 0.0)), std::invalid_argument);
  EXPECT_THROW(Reflection(Vector3D(1.0, 2.0, -3.0)), std::invalid_argument);
}

TEST(TestReflection,PlaneXY)
{
  Reflection reflection = Reflection(Vector3D{0.0, 0.0, 1.0});

  const Vector3D x = Vector3D(1.0,0.0,0.0);
  const Vector3D y = Vector3D(0.0,1.0,0.0);
  const Vector3D z = Vector3D(0.0,0.0,1.0);

  EXPECT_EQ(reflection * x, x );
  EXPECT_EQ(reflection(x), x );
  EXPECT_EQ(reflection * y, y );
  EXPECT_EQ(reflection(y), y );

  EXPECT_EQ(reflection * z, (-1.0)*z );
  EXPECT_EQ(reflection(z), (-1.0)*z );
}

TEST(TestReflection,Plane111)
{
  const Vector3D x = Vector3D(1.0,0.0,0.0);

  Vector3D v111 = Vector3D(1.0, 1.0, 1.0);
  v111 = (1.0/sqrt(3)) * v111;
  Reflection reflection = Reflection(v111);

  EXPECT_EQ(reflection * x, Vector3D(1.0/3.0, -2.0/3.0, -2.0/3.0) );
}