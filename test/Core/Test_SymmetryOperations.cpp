#include <gtest/gtest.h>
#include <gmock/gmock.h>

#include <Core/SymmetryOperations.h>

namespace
{
  static const double pi = 3.14159265358979323846;
}

using namespace Core::Geometry;

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

  EXPECT_NEAR(v100.x, rotation(v100).x, 1e-15 );
  EXPECT_NEAR(v100.y, rotation(v100).y, 1e-15 );
  EXPECT_NEAR(v100.z, rotation(v100).z, 1e-15 );

  EXPECT_NEAR(v010.x, rotation(v010).x, 1e-15 );
  EXPECT_NEAR(v010.y, rotation(v010).y, 1e-15 );
  EXPECT_NEAR(v010.z, rotation(v010).z, 1e-15 );

  EXPECT_NEAR(v001.x, rotation(v001).x, 1e-15 );
  EXPECT_NEAR(v001.y, rotation(v001).y, 1e-15 );
  EXPECT_NEAR(v001.z, rotation(v001).z, 1e-15 );
}

TEST(TestRotation,90Degrees)
{
  const Vector3D x = pi / 2.0 * Vector3D(1.0, 0.0, 0.0);
  const Vector3D y = pi / 2.0 * Vector3D(0.0, 1.0, 0.0);
  const Vector3D z = pi / 2.0 * Vector3D(0.0, 0.0, 1.0);

  Rotation rotation_x = Rotation(x);
  Rotation rotation_y = Rotation(y);
  Rotation rotation_z = Rotation(z);

  // rotation with x
  EXPECT_NEAR(x.x, rotation_x(x).x, 1e-15 );
  EXPECT_NEAR(x.y, rotation_x(x).y, 1e-15 );
  EXPECT_NEAR(x.z, rotation_x(x).z, 1e-15 );

  EXPECT_NEAR(z.x, rotation_x(y).x, 1e-15 );
  EXPECT_NEAR(z.y, rotation_x(y).y, 1e-15 );
  EXPECT_NEAR(z.z, rotation_x(y).z, 1e-15 );
  
  EXPECT_NEAR(-y.x, rotation_x(z).x, 1e-15 );
  EXPECT_NEAR(-y.y, rotation_x(z).y, 1e-15 );
  EXPECT_NEAR(-y.z, rotation_x(z).z, 1e-15 );

  // rotation with y
  EXPECT_NEAR(-z.x, rotation_y(x).x, 1e-15 );
  EXPECT_NEAR(-z.y, rotation_y(x).y, 1e-15 );
  EXPECT_NEAR(-z.z, rotation_y(x).z, 1e-15 );

  EXPECT_NEAR(y.x, rotation_y(y).x, 1e-15 );
  EXPECT_NEAR(y.y, rotation_y(y).y, 1e-15 );
  EXPECT_NEAR(y.z, rotation_y(y).z, 1e-15 );
  
  EXPECT_NEAR(x.x, rotation_y(z).x, 1e-15 );
  EXPECT_NEAR(x.y, rotation_y(z).y, 1e-15 );
  EXPECT_NEAR(x.z, rotation_y(z).z, 1e-15 );

  // rotation with z
  EXPECT_NEAR(y.x, rotation_z(x).x, 1e-15 );
  EXPECT_NEAR(y.y, rotation_z(x).y, 1e-15 );
  EXPECT_NEAR(y.z, rotation_z(x).z, 1e-15 );

  EXPECT_NEAR(-x.x, rotation_z(y).x, 1e-15 );
  EXPECT_NEAR(-x.y, rotation_z(y).y, 1e-15 );
  EXPECT_NEAR(-x.z, rotation_z(y).z, 1e-15 );
  
  EXPECT_NEAR(z.x, rotation_z(z).x, 1e-15 );
  EXPECT_NEAR(z.y, rotation_z(z).y, 1e-15 );
  EXPECT_NEAR(z.z, rotation_z(z).z, 1e-15 );
}