#ifndef LATTICEGEN_SYMMETRYOPERATIONS_H_
#define LATTICEGEN_SYMMETRYOPERATIONS_H_

#include "Point3D.h"
#include "Matrix3D.h"
#include <array>

namespace Geometry
{
  class ImproperRotation;

  class Rotation
  {
  public:
    Rotation(const Vector3D& rotation_vector);
    Vector3D operator()(const Vector3D& vector) const;

    operator const Matrix3D&() const
    {
      return rotation_matrix;
    }

  private:
    Matrix3D rotation_matrix;
    friend class ImproperRotation;
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

    operator const Matrix3D&() const
    {
      return reflection_matrix;
    }
  private:
    Matrix3D reflection_matrix;
    friend class ImproperRotation;
  };

  inline Vector3D operator*(const Reflection& reflection, const Vector3D& vector)
  {
    return reflection(vector);
  }

  class ImproperRotation
  {
  public:
    ImproperRotation(const Vector3D& rotation_vector);
    Vector3D operator()(const Vector3D& vector) const;

    operator const Matrix3D&() const
    {
      return transformation_matrix;
    }
  private:
    Matrix3D transformation_matrix;
  };

  inline Vector3D operator*(const ImproperRotation& improper_rotation, const Vector3D& vector)
  {
    return improper_rotation(vector);
  }

  class Inversion
  {
  public:
    Vector3D operator()(const Vector3D& vector) const
    {
      return inversion_matrix * vector;
    }
    operator const Matrix3D&() const
    {
      return inversion_matrix;
    }
    static Matrix3D inversion_matrix;
  };

  inline Vector3D operator*(const Inversion& inversion, const Vector3D& vector)
  {
    return inversion(vector);
  }
}
 
#endif // LATTICEGEN_SYMMETRYOPERATIONS_H_