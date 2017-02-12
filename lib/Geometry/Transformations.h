#ifndef LATTICEGEN_TRANSFORMATIONS_H_
#define LATTICEGEN_TRANSFORMATIONS_H_

#include "Point3D.h"
#include "Matrix3D.h"
#include <array>

#include <Core/ComparisonHelpers.h>

namespace Geometry
{
  class Transformation
  {
  public:
    enum Type
    {
      Identity,
      Rotation,
      Reflection,
      ImproperRotation,
      Inversion
    };

    operator const Matrix3D&() const
    {
      return transformation_matrix;
    }
    Vector3D operator()(const Vector3D& vector) const
    {
      return transformation_matrix * vector;
    }

    const Type type;
    const Matrix3D transformation_matrix;

  protected:
    Transformation(const Type type_, const Matrix3D& matrix_)
      : type(type_), transformation_matrix(matrix_)
    {}
  };

  inline Vector3D operator*(const Transformation& symmetry_element, const Vector3D& vector)
  {
    return symmetry_element.transformation_matrix * vector;
  }
  inline bool operator==(const Transformation& lhs, const Transformation& rhs)
  {
    return lhs.type == rhs.type && lhs.transformation_matrix == rhs.transformation_matrix;
  }

  class Identity : public Transformation
  {
  public:
    Identity();
  };

  class Rotation : public Transformation
  {
  public:
    explicit Rotation(const Vector3D& rotation_vector);
  };


  class Reflection : public Transformation
  {
  public:
    explicit Reflection(const Vector3D& reflection_plane);
  };


  class ImproperRotation : public Transformation
  {
  public:
    explicit ImproperRotation(const Vector3D& rotation_vector);
  };

  class Inversion : public Transformation
  {
  public:
    Inversion();
  };
}
 
#endif // LATTICEGEN_TRANSFORMATIONS_H_
