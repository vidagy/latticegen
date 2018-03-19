#include <TestUtils/base.h>

#include <Geometry/Transformations.h>

namespace
{
  static const double pi = 3.14159265358979323846;
}

using namespace Geometry;

TEST(TestTransformations,RotationCtorThrows)
{
  EXPECT_THROW(Rotation(Point3D(0.0, 0.0, 0.0)), std::invalid_argument);
}

TEST(TestTransformations,RotationIdentity)
{
  const auto v = 2.0 * pi / sqrt(1.0 + 4.0 + 9.0) * Point3D(1.0, 2.0, 3.0);
  auto rotation = Rotation(v);

  const auto v100 = Point3D(1.0, 0.0, 0.0);
  const auto v010 = Point3D(0.0, 1.0, 0.0);
  const auto v001 = Point3D(0.0, 0.0, 1.0);

  EXPECT_POINTS_CLOSE(v100, rotation(v100));
  EXPECT_POINTS_CLOSE(v010, rotation(v010));
  EXPECT_POINTS_CLOSE(v001, rotation(v001));
}

TEST(TestTransformations,RotationDegrees90)
{
  const auto x = pi / 2.0 * Point3D(1.0, 0.0, 0.0);
  const auto z = pi / 2.0 * Point3D(0.0, 0.0, 1.0);
  const auto y = pi / 2.0 * Point3D(0.0, 1.0, 0.0);

  Rotation rotation_x = Rotation(x);
  Rotation rotation_y = Rotation(y);
  Rotation rotation_z = Rotation(z);

  // rotation with x
  EXPECT_POINTS_CLOSE(x, rotation_x(x));
  EXPECT_POINTS_CLOSE(z, rotation_x(y));
  EXPECT_POINTS_CLOSE(-1.0 * y, rotation_x(z));

  // rotation with y
  EXPECT_POINTS_CLOSE(-1.0 * z, rotation_y(x));
  EXPECT_POINTS_CLOSE(y, rotation_y(y));
  EXPECT_POINTS_CLOSE(x, rotation_y(z));

  // rotation with z
  EXPECT_POINTS_CLOSE(y, rotation_z(x));
  EXPECT_POINTS_CLOSE(-1.0 * x, rotation_z(y));
  EXPECT_POINTS_CLOSE(z, rotation_z(z));
}

TEST(TestTransformations,ReflectionCtorThrows)
{
  EXPECT_THROW(Reflection(Point3D(0.0, 0.0, 0.0)), std::invalid_argument);
  EXPECT_THROW(Reflection(Point3D(1.0, 2.0, -3.0)), std::invalid_argument);
}

TEST(TestTransformations,ReflectionPlaneXY)
{
  Reflection reflection = Reflection(Point3D{0.0, 0.0, 1.0});

  const auto x = Point3D(1.0, 0.0, 0.0);
  const auto y = Point3D(0.0, 1.0, 0.0);
  const auto z = Point3D(0.0, 0.0, 1.0);

  EXPECT_POINTS_CLOSE(reflection(x), x);
  EXPECT_POINTS_CLOSE(reflection(y), y);

  EXPECT_POINTS_CLOSE(reflection(z), (-1.0) * z);
}

TEST(TestTransformations,ReflectionPlane111)
{
  const auto x = Point3D(1.0, 0.0, 0.0);

  auto v111 = Point3D(1.0, 1.0, 1.0);
  v111 = (1.0/sqrt(3)) * v111;
  Reflection reflection = Reflection(v111);

  EXPECT_POINTS_CLOSE(reflection(x), Point3D(1.0 / 3.0, -2.0 / 3.0, -2.0 / 3.0));
}

TEST(TestTransformations,ImproperRotationRotateAndReflect)
{
  const auto xz = Point3D(1.0, 0.0, 1.0);

  auto z90 = Point3D(0.0, 0.0, pi / 2.0);
  ImproperRotation improper_rotation = ImproperRotation(z90);

  EXPECT_POINTS_CLOSE(improper_rotation * xz, Point3D(0.0, 1.0, -1.0));
}

namespace 
{
  void test_inversion(const Point3D &vector) {
    auto rotation = Rotation(vector / vector.norm() * pi);
    auto reflection = Reflection(vector / vector.norm());

    const auto rotref = rotation * reflection;
    const auto refrot = reflection * rotation;

    for (auto x = -2.0; x <= 2.0; x += 0.6)
    {
      for (auto y = -2.0; y <= 2.0; y += 0.7)
      {
        for (auto z = -2.0; z <= 2.0; z += 0.8)
        {
          auto v = Point3D(x, y, z);
          EXPECT_POINTS_CLOSE(rotref(v), Inversion()(v));
          EXPECT_POINTS_CLOSE(refrot(v), Inversion()(v));
        }
      }
    }
  }
}

TEST(TestTransformations,InversionRotate180AndReflectIsInversion)
{
  test_inversion(Point3D(1.0, 0.0, 0.0));
  test_inversion(Point3D(0.0, 1.0, 0.0));
  test_inversion(Point3D(0.0, 0.0, 1.0));

  test_inversion(Point3D(1.0, 1.0, 0.0));
  test_inversion(Point3D(1.0, 0.0, 1.0));
  test_inversion(Point3D(0.0, 1.0, 1.0));

  test_inversion(Point3D(-1.0, 1.0, 0.0));
  test_inversion(Point3D(-1.0, 0.0, 1.0));
  test_inversion(Point3D(0.0, -1.0, 1.0));

  test_inversion(Point3D(1.0, -1.0, 0.0));
  test_inversion(Point3D(1.0, 0.0, -1.0));
  test_inversion(Point3D(0.0, 1.0, -1.0));

  test_inversion(Point3D(-1.0, -1.0, 0.0));
  test_inversion(Point3D(-1.0, 0.0, -1.0));
  test_inversion(Point3D(0.0, -1.0, -1.0));

  test_inversion(Point3D(1.0, 1.0, 1.0));
  test_inversion(Point3D(-1.0, 1.0, 1.0));
  test_inversion(Point3D(1.0, -1.0, 1.0));
  test_inversion(Point3D(1.0, 1.0, -1.0));
  test_inversion(Point3D(-1.0, -1.0, 1.0));
  test_inversion(Point3D(-1.0, 1.0, -1.0));
  test_inversion(Point3D(1.0, -1.0, -1.0));
  test_inversion(Point3D(-1.0, -1.0, -1.0));

  test_inversion(Point3D(1.0, 2.0, 3.0));
}

namespace
{
  void test_identity(const Point3D &vector) {
    const auto reflection = Reflection(vector / vector.norm());
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
          const auto v = Point3D(x, y, z);
          EXPECT_POINTS_CLOSE(reflect2(v), Identity()(v));

          EXPECT_POINTS_CLOSE(rotinvrot(v), Identity()(v));
          EXPECT_POINTS_CLOSE(invrotrot(v), Identity()(v));

          EXPECT_POINTS_CLOSE(irotinvirot(v), Identity()(v));
          EXPECT_POINTS_CLOSE(invirotirot(v), Identity()(v));
        }
      }
    }
  }
}

TEST(TestTransformations,Identity)
{
  test_identity(Point3D(1.0, 0.0, 0.0));
  test_identity(Point3D(0.0, 1.0, 0.0));
  test_identity(Point3D(0.0, 0.0, 1.0));

  test_identity(Point3D(1.0, 1.0, 0.0));
  test_identity(Point3D(1.0, 0.0, 1.0));
  test_identity(Point3D(0.0, 1.0, 1.0));

  test_identity(Point3D(-1.0, 1.0, 0.0));
  test_identity(Point3D(-1.0, 0.0, 1.0));
  test_identity(Point3D(0.0, -1.0, 1.0));

  test_identity(Point3D(1.0, -1.0, 0.0));
  test_identity(Point3D(1.0, 0.0, -1.0));
  test_identity(Point3D(0.0, 1.0, -1.0));

  test_identity(Point3D(-1.0, -1.0, 0.0));
  test_identity(Point3D(-1.0, 0.0, -1.0));
  test_identity(Point3D(0.0, -1.0, -1.0));

  test_identity(Point3D(1.0, 1.0, 1.0));
  test_identity(Point3D(-1.0, 1.0, 1.0));
  test_identity(Point3D(1.0, -1.0, 1.0));
  test_identity(Point3D(1.0, 1.0, -1.0));
  test_identity(Point3D(-1.0, -1.0, 1.0));
  test_identity(Point3D(-1.0, 1.0, -1.0));
  test_identity(Point3D(1.0, -1.0, -1.0));
  test_identity(Point3D(-1.0, -1.0, -1.0));

  test_identity(Point3D(1.0, 2.0, 3.0));
}

TEST(TestTransformations, Multiplication)
{
  auto vector = Point3D{1.0, 2.0, 3.0};
  vector.normalize();

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
  //const double step = pi / 12.0;
  //for (double alpha = step; alpha <= 2.0 * pi; alpha += step) {
  //  for (double beta = step; beta <= 2.0 * pi; beta += step) {
  //    for (double gamma = step; gamma <= 2.0 * pi; gamma += step) {
  double alpha = 1.0 * pi;
  double beta = 2.0 * pi;
  double gamma = 1.0 * pi;
        auto rot_euler = Rotation(alpha, beta, gamma);

        auto rot_alpha = Rotation(alpha * Point3D{0.0, 0.0, 1.0});
        auto rot_beta = Rotation(beta * Point3D{0.0, 1.0, 0.0});
        auto rot_gamma = Rotation(gamma * Point3D{0.0, 0.0, 1.0});

        auto rot_components = rot_gamma * rot_beta * rot_alpha;
  EXPECT_TRUE(rot_components == rot_euler)
                << " for alpha = " << alpha / pi << " beta = " << beta / pi << " gamma = " << gamma / pi
                << "\nrot_euler = " << rot_euler << "\nrot_components = " << rot_components;
  //   }
  // }
  //}
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
