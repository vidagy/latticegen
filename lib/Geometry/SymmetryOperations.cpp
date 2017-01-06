#include "SymmetryOperations.h" 

#include <Core/ComparisonHelpers.h>

#include <stdexcept>

using namespace Core;
using namespace Geometry;

namespace
{
  Rotation::RotationMatrix get_rotation_matrix(const Vector3D& rotation_vector)
  {
    const double angle = rotation_vector.getLength();
    
    if (Core::nearlyZero(angle))
      throw std::invalid_argument("Rotation can not be created from null vector");

    const Vector3D axis = rotation_vector * (1.0 / angle);

    Rotation::RotationMatrix result = {
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
  return Vector3D(
      rotation_matrix[0][0] * vector.x + rotation_matrix[0][1] * vector.y + rotation_matrix[0][2] * vector.z,
      rotation_matrix[1][0] * vector.x + rotation_matrix[1][1] * vector.y + rotation_matrix[1][2] * vector.z,
      rotation_matrix[2][0] * vector.x + rotation_matrix[2][1] * vector.y + rotation_matrix[2][2] * vector.z
    );
}
