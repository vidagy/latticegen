#include "SymmetryOperations.h" 

using namespace Core;
using namespace Geometry;

namespace
{
  Matrix3D get_rotation_matrix(const Vector3D& rotation_vector)
  {
    const double angle = rotation_vector.length();
    
    if (Core::nearlyZero(angle))
      throw std::invalid_argument("Rotation can not be created from null vector");

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
  : rotation_matrix(get_rotation_matrix(rotation_vector))
{}

Vector3D Rotation::operator()(const Vector3D& vector) const
{
  return rotation_matrix * vector;
}

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
    Matrix3D unit_matrix = {1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0};
    return unit_matrix -2.0 * get_diad(vector);
  }
}

Reflection::Reflection(const Vector3D& reflection_plane)
  : reflection_matrix(get_reflection_matrix(reflection_plane))
{
  if (! equalsWithTolerance(reflection_plane.length(), 1.0))
    throw std::invalid_argument("Non-unit vector on input of Reflection: " + std::to_string(reflection_plane));
}

Vector3D Reflection::operator()(const Vector3D& vector) const
{
  return reflection_matrix * vector;
}

ImproperRotation::ImproperRotation(const Vector3D& rotation_vector)
  : transformation_matrix(
    Reflection(rotation_vector/rotation_vector.length()).reflection_matrix *
    Rotation(rotation_vector).rotation_matrix )
{
}

Vector3D ImproperRotation::operator()(const Vector3D& vector) const
{
  return transformation_matrix * vector;
}

const Matrix3D Identity::identity_matrix = Matrix3D{
  1.0, 0.0, 0.0,
  0.0, 1.0, 0.0,
  0.0, 0.0, 1.0
};

const Matrix3D Inversion::inversion_matrix = Matrix3D{
  -1.0, 0.0, 0.0,
  0.0, -1.0, 0.0,
  0.0, 0.0, -1.0
};
