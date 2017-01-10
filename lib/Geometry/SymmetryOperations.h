#ifndef LATTICEGEN_SYMMETRYOPERATIONS_H_
#define LATTICEGEN_SYMMETRYOPERATIONS_H_

#include "Point3D.h"
#include "Matrix3D.h"
#include <array>

namespace Geometry
{
  class Rotation
  {
  public:
    Rotation(const Vector3D& rotation_vector);
    Vector3D operator()(const Vector3D& vector) const;
  private:
    Matrix3D rotation_matrix;
  };

  inline Vector3D operator*(const Rotation& rotation, const Vector3D& vector)
  {
    return rotation(vector);
  }

  class Reflection
  {
  public:
    Reflection(const Vector3D& reflection_plane);
    Vector3D operator()(const Vector3D& vector) const;
  private:
    Matrix3D reflection_matrix;
  };

  inline Vector3D operator*(const Reflection& reflection, const Vector3D& vector)
  {
    return reflection(vector);
  }
}
 
#endif // LATTICEGEN_SYMMETRYOPERATIONS_H_