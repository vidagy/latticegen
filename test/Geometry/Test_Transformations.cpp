#include <TestUtils/base.h>

#include <Geometry/Transformations.h>

namespace
{
  static const double pi = 3.14159265358979323846;
}

using namespace Geometry;

TEST(TestTransformations,RotationCtorThrows)
{
  EXPECT_THROW(Rotation(Vector3D(0.0, 0.0, 0.0)), std::invalid_argument);
}

TEST(TestTransformations,RotationIdentity)
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

TEST(TestTransformations,RotationDegrees90)
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

TEST(TestTransformations,ReflectionCtorThrows)
{
  EXPECT_THROW(Reflection(Vector3D(0.0, 0.0, 0.0)), std::invalid_argument);
  EXPECT_THROW(Reflection(Vector3D(1.0, 2.0, -3.0)), std::invalid_argument);
}

TEST(TestTransformations,ReflectionPlaneXY)
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

TEST(TestTransformations,ReflectionPlane111)
{
  const Vector3D x = Vector3D(1.0,0.0,0.0);

  Vector3D v111 = Vector3D(1.0, 1.0, 1.0);
  v111 = (1.0/sqrt(3)) * v111;
  Reflection reflection = Reflection(v111);

  EXPECT_EQ(reflection * x, Vector3D(1.0/3.0, -2.0/3.0, -2.0/3.0) );
}

TEST(TestTransformations,ImproperRotationRotateAndReflect)
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

TEST(TestTransformations,InversionRotate180AndReflectIsInversion)
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
    const auto reflection = Reflection(vector / vector.length());
    const auto rotation = Rotation(vector * pi * 0.783);
    const auto inv_rotation = Rotation(-1.0 * vector * pi * 0.783);
    const auto improper_rotation = ImproperRotation(vector * 1.123);
    const auto inv_improper_rotation = ImproperRotation(-1.0 * vector * 1.123);

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
          const auto v = Vector3D(x, y, z);
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

TEST(TestTransformations,Identity)
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

TEST(TestTransformations, Multiplication)
{
  auto vector = Vector3D{1.0, 2.0, 3.0};
  vector /= vector.length();

  const auto identity = Identity();
  const auto inversion = Inversion();
  const auto reflection = Reflection(vector);
  const auto rotation = Rotation(vector * pi * 0.783);
  const auto inv_rotation = Rotation(-1.0 * vector * pi * 0.783);
  const auto improper_rotation = ImproperRotation(vector * 1.123);
  const auto inv_improper_rotation = ImproperRotation(-1.0 * vector * 1.123);

  // multiplication with identity
  EXPECT_EQ(identity * identity, identity);
  EXPECT_EQ(identity * inversion, inversion);
  EXPECT_EQ(inversion * identity, inversion);
  EXPECT_EQ(identity * rotation, rotation);
  EXPECT_EQ(rotation * identity, rotation);
  EXPECT_EQ(identity * reflection, reflection);
  EXPECT_EQ(reflection * identity, reflection);
  EXPECT_EQ(identity * improper_rotation, improper_rotation);
  EXPECT_EQ(improper_rotation * identity, improper_rotation);

  // multiplication with inverse
  EXPECT_EQ(inversion * inversion, identity);
  EXPECT_EQ(reflection * reflection, identity);
  EXPECT_EQ(rotation * inv_rotation, identity);
  EXPECT_EQ(improper_rotation * inv_improper_rotation, identity);

  // with same kind
  EXPECT_EQ(rotation * rotation, Rotation(2.0 * vector * pi * 0.783));
  EXPECT_EQ(improper_rotation * improper_rotation, Rotation(2.0 * vector * 1.123));

  // mixed
  EXPECT_EQ(rotation * improper_rotation, ImproperRotation(vector * (1.123 + pi * 0.783)));
  EXPECT_EQ((reflection * improper_rotation).type, Transformation::Rotation);
}

TEST(TestTransformations, RotationsAsEulerMultiplied)
{
  const double step = pi / 12.0;
  for (double alpha = step; alpha <= 2.0 * pi; alpha += step) {
    for (double beta = step; beta <= 2.0 * pi; beta += step) {
      for (double gamma = step; gamma <= 2.0 * pi; gamma += step) {
        auto rot_euler = Rotation(alpha, beta, gamma);

        auto rot_alpha = Rotation(alpha * Point3D{0.0, 0.0, 1.0});
        auto rot_beta = Rotation(beta * Point3D{0.0, 1.0, 0.0});
        auto rot_gamma = Rotation(gamma * Point3D{0.0, 0.0, 1.0});

        auto rot_components = rot_gamma * rot_beta * rot_alpha;
        EXPECT_TRUE(rot_components.transformation_matrix == rot_euler.transformation_matrix)
                << " for alpha = " << alpha / pi << " beta = " << beta / pi << " gamma = " << gamma / pi
                << "\nrot_euler = " << rot_euler << "\nrot_components = " << rot_components;
      }
    }
  }
}


TEST(TestTransformations, RotationEulerMatrices)
{
  const double step = pi / 12.0;
  for (double alpha = 0.0 * pi; alpha <= 2.0 * pi; alpha += step) {
    for (double beta = 0.0 * pi; beta <= 2.0 * pi; beta += step) {
      for (double gamma = 0.0 * pi; gamma <= 2.0 * pi; gamma += step) {
        auto rot = Rotation(alpha, beta, gamma);
        auto euler = rot.get_euler_angles();
        auto rot2 = Rotation(euler.alpha, euler.beta, euler.gamma);
        auto euler2 = rot2.get_euler_angles();

        EXPECT_TRUE(rot == rot2)
                << "Rot original = \n" << rot << "\nRot new = \n" << rot2
                << "\nOriginal euler angles = "
                << "alpha " << euler.alpha / pi << " beta " << euler.beta / pi << " gamma " << euler.gamma / pi
                << "\nNew euler angles = "
                << "alpha " << euler2.alpha / pi << " beta2 " << euler.beta / pi << " gamma2 " << euler.gamma / pi;
      }
    }
  }
}
