#include <gtest/gtest.h>
#include <gmock/gmock.h>

#include <Geometry/SymmetryOperations.h>

namespace
{
  static const double pi = 3.14159265358979323846;
}

using namespace Geometry;

TEST(TestSymmetyOperations,RotationCtorThrows)
{
  EXPECT_THROW(Rotation(Vector3D(0.0, 0.0, 0.0)), std::invalid_argument);
}

TEST(TestSymmetyOperations,RotationIdentity)
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

TEST(TestSymmetyOperations,RotationDegrees90)
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

TEST(TestSymmetyOperations,ReflectionCtorThrows)
{
  EXPECT_THROW(Reflection(Vector3D(0.0, 0.0, 0.0)), std::invalid_argument);
  EXPECT_THROW(Reflection(Vector3D(1.0, 2.0, -3.0)), std::invalid_argument);
}

TEST(TestSymmetyOperations,ReflectionPlaneXY)
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

TEST(TestSymmetyOperations,ReflectionPlane111)
{
  const Vector3D x = Vector3D(1.0,0.0,0.0);

  Vector3D v111 = Vector3D(1.0, 1.0, 1.0);
  v111 = (1.0/sqrt(3)) * v111;
  Reflection reflection = Reflection(v111);

  EXPECT_EQ(reflection * x, Vector3D(1.0/3.0, -2.0/3.0, -2.0/3.0) );
}

TEST(TestSymmetyOperations,ImproperRotationRotateAndReflect)
{
  const Vector3D xz = Vector3D(1.0,0.0,1.0);

  Vector3D z90 = Vector3D(0.0, 0.0, pi/2.0);
  ImproperRotation improper_rotation = ImproperRotation(z90);

  EXPECT_EQ(improper_rotation * xz, Vector3D(0.0, 1.0, -1.0) );
}

namespace 
{
  void test_inversion(const Vector3D& vector)
  {
    Rotation rotation = Rotation(vector/vector.length()*pi);
    Reflection reflection = Reflection(vector/vector.length());

    const auto rotref = rotation * reflection;
    const auto refrot = reflection * rotation;

    for (auto x = -2.0; x <= 2.0; x += 0.6)
    {
      for (auto y = -2.0; y <= 2.0; y += 0.7)
      {
        for (auto z = -2.0; z <= 2.0; z += 0.8)
        {
          Vector3D v = Vector3D(x,y,z);
          EXPECT_EQ(rotref * v, Inversion() * v);
          EXPECT_EQ(refrot * v, Inversion() * v);
        }
      }
    }
  }
}

TEST(TestSymmetyOperations,InversionRotate180AndReflectIsInversion)
{
  test_inversion(Vector3D(1.0,0.0,0.0));
  test_inversion(Vector3D(0.0,1.0,0.0));
  test_inversion(Vector3D(0.0,0.0,1.0));

  test_inversion(Vector3D( 1.0, 1.0, 0.0));
  test_inversion(Vector3D( 1.0, 0.0, 1.0));
  test_inversion(Vector3D( 0.0, 1.0, 1.0));

  test_inversion(Vector3D(-1.0, 1.0, 0.0));
  test_inversion(Vector3D(-1.0, 0.0, 1.0));
  test_inversion(Vector3D( 0.0,-1.0, 1.0));

  test_inversion(Vector3D( 1.0,-1.0, 0.0));
  test_inversion(Vector3D( 1.0, 0.0,-1.0));
  test_inversion(Vector3D( 0.0, 1.0,-1.0));

  test_inversion(Vector3D(-1.0,-1.0, 0.0));
  test_inversion(Vector3D(-1.0, 0.0,-1.0));
  test_inversion(Vector3D( 0.0,-1.0,-1.0));

  test_inversion(Vector3D( 1.0, 1.0, 1.0));
  test_inversion(Vector3D(-1.0, 1.0, 1.0));
  test_inversion(Vector3D( 1.0,-1.0, 1.0));
  test_inversion(Vector3D( 1.0, 1.0,-1.0));
  test_inversion(Vector3D(-1.0,-1.0, 1.0));
  test_inversion(Vector3D(-1.0, 1.0,-1.0));
  test_inversion(Vector3D( 1.0,-1.0,-1.0));
  test_inversion(Vector3D(-1.0,-1.0,-1.0));

  test_inversion(Vector3D(1.0,2.0,3.0));
}

namespace
{
  void test_identity(const Vector3D& vector)
  {
    Reflection reflection = Reflection(vector/vector.length());
    Rotation rotation = Rotation(vector*pi*0.783);
    Rotation inv_rotation = Rotation(-1.0 * vector*pi*0.783);
    ImproperRotation improper_rotation = ImproperRotation(vector * 1.123);
    ImproperRotation inv_improper_rotation = ImproperRotation(-1.0 * vector * 1.123);

    const auto reflect2 = reflection * reflection;
    const auto rotinvrot = rotation * inv_rotation;
    const auto invrotrot = inv_rotation * rotation;
    const auto irotinvirot = improper_rotation * inv_improper_rotation ;
    const auto invirotirot = inv_improper_rotation * improper_rotation ;

    for (auto x = -2.0; x <= 2.0; x += 0.6)
    {
      for (auto y = -2.0; y <= 2.0; y += 0.7)
      {
        for (auto z = -2.0; z <= 2.0; z += 0.8)
        {
          Vector3D v = Vector3D(x,y,z);
          EXPECT_EQ(reflect2 * v, Identity() * v);

          EXPECT_EQ(rotinvrot * v, Identity() * v);
          EXPECT_EQ(invrotrot * v, Identity() * v);

          EXPECT_EQ(irotinvirot * v, Identity() * v);
          EXPECT_EQ(invirotirot * v, Identity() * v);
        }
      }
    }
  }
}

TEST(TestSymmetyOperations,Identity)
{
  test_identity(Vector3D(1.0,0.0,0.0));
  test_identity(Vector3D(0.0,1.0,0.0));
  test_identity(Vector3D(0.0,0.0,1.0));

  test_identity(Vector3D( 1.0, 1.0, 0.0));
  test_identity(Vector3D( 1.0, 0.0, 1.0));
  test_identity(Vector3D( 0.0, 1.0, 1.0));

  test_identity(Vector3D(-1.0, 1.0, 0.0));
  test_identity(Vector3D(-1.0, 0.0, 1.0));
  test_identity(Vector3D( 0.0,-1.0, 1.0));

  test_identity(Vector3D( 1.0,-1.0, 0.0));
  test_identity(Vector3D( 1.0, 0.0,-1.0));
  test_identity(Vector3D( 0.0, 1.0,-1.0));

  test_identity(Vector3D(-1.0,-1.0, 0.0));
  test_identity(Vector3D(-1.0, 0.0,-1.0));
  test_identity(Vector3D( 0.0,-1.0,-1.0));

  test_identity(Vector3D( 1.0, 1.0, 1.0));
  test_identity(Vector3D(-1.0, 1.0, 1.0));
  test_identity(Vector3D( 1.0,-1.0, 1.0));
  test_identity(Vector3D( 1.0, 1.0,-1.0));
  test_identity(Vector3D(-1.0,-1.0, 1.0));
  test_identity(Vector3D(-1.0, 1.0,-1.0));
  test_identity(Vector3D( 1.0,-1.0,-1.0));
  test_identity(Vector3D(-1.0,-1.0,-1.0));

  test_identity(Vector3D(1.0,2.0,3.0));
}
