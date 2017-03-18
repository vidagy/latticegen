#include "Transformations.h"
#include <Core/Exceptions.h>

using namespace Core;
using namespace Geometry;

namespace
{
  static const Matrix3D identity_matrix = {
    1.0, 0.0, 0.0,
    0.0, 1.0, 0.0,
    0.0, 0.0, 1.0
  };
}

Identity::Identity()
  : Transformation(Transformation::Identity, identity_matrix)
{}

namespace
{
  Matrix3D get_rotation_matrix(const Vector3D& rotation_vector)
  {
    const double angle = rotation_vector.length();

    if (nearlyZero(angle))
      THROW_INVALID_ARGUMENT("Rotation can not be created from null vector");

    const Vector3D axis = rotation_vector * (1.0 / angle);

    Matrix3D result = {
        cos(angle)+pow(axis.x,2)*(1.0-cos(angle)) , 
        axis.x*axis.y*(1.0-cos(angle))-axis.z*sin(angle) ,
        axis.x*axis.z*(1.0-cos(angle))+axis.y*sin(angle) ,

        axis.y*axis.x*(1.0-cos(angle))+axis.z*sin(angle) ,
        cos(angle)+pow(axis.y,2)*(1.0-cos(angle)) , 
        axis.y*axis.z*(1.0-cos(angle))-axis.x*sin(angle) ,

        axis.z*axis.x*(1.0-cos(angle))-axis.y*sin(angle) ,
        axis.z*axis.y*(1.0-cos(angle))+axis.x*sin(angle) ,
        cos(angle)+pow(axis.z,2)*(1.0-cos(angle)) 
    };

    return result;
  }
}

Rotation::Rotation(const Vector3D& rotation_vector)
  : Transformation(Transformation::Rotation, get_rotation_matrix(rotation_vector))
{}

namespace
{
  Matrix3D get_diad(const Vector3D& vector)
  {
    return Matrix3D {
      vector.x *vector.x, vector.x *vector.y, vector.x *vector.z,
      vector.y *vector.x, vector.y *vector.y, vector.y *vector.z,
      vector.z *vector.x, vector.z *vector.y, vector.z *vector.z
    };
  }

  Matrix3D get_reflection_matrix(const Vector3D& vector)
  {
    const Matrix3D unit_matrix = Identity().transformation_matrix;
    return unit_matrix -2.0 * get_diad(vector);
  }
}

Reflection::Reflection(const Vector3D& reflection_plane)
  : Transformation(Transformation::Reflection, get_reflection_matrix(reflection_plane))
{
  if (! equalsWithTolerance(reflection_plane.length(), 1.0))
    THROW_INVALID_ARGUMENT("Non-unit vector on input of Reflection: " + std::to_string(reflection_plane));
}

ImproperRotation::ImproperRotation(const Vector3D& rotation_vector)
  : Transformation(Transformation::ImproperRotation,
      Geometry::Reflection(rotation_vector/rotation_vector.length()).transformation_matrix *
      Geometry::Rotation(rotation_vector).transformation_matrix
    )
{}

namespace
{
  static const Matrix3D inversion_matrix = {
    -1.0, 0.0, 0.0,
    0.0, -1.0, 0.0,
    0.0, 0.0, -1.0
  };
}

Inversion::Inversion()
  : Transformation(Transformation::Inversion, inversion_matrix)
{}
