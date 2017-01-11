#ifndef LATTICEGEN_SYMMETRYELEMENTS_H_
#define LATTICEGEN_SYMMETRYELEMENTS_H_

#include "Point3D.h"
#include "Matrix3D.h"
#include <array>

namespace Geometry
{
  class SymmetryElement
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

    virtual Type get_type() const { return type; }

    operator const Matrix3D&() const
    {
      return transformation_matrix;
    }
    Vector3D operator()(const Vector3D& vector) const
    {
      return transformation_matrix * vector;
    }

    const Matrix3D transformation_matrix;

  protected:
    SymmetryElement(const Type type_, const Matrix3D& matrix_)
      : transformation_matrix(matrix_), type(type_)
    {}
  private:
    const Type type;
  };

  inline Vector3D operator*(const SymmetryElement& symmetry_element, const Vector3D& vector)
  {
    return symmetry_element.transformation_matrix * vector;
  }

  class Identity : public SymmetryElement
  {
  public:
    Identity()
      : SymmetryElement(SymmetryElement::Identity, Matrix3D{
        1.0, 0.0, 0.0,
        0.0, 1.0, 0.0,
        0.0, 0.0, 1.0
      })
    {}
  };

  class Rotation : public SymmetryElement
  {
  public:
    Rotation(const Vector3D& rotation_vector);
  };


  class Reflection : public SymmetryElement
  {
  public:
    Reflection(const Vector3D& reflection_plane);
  };


  class ImproperRotation : public SymmetryElement
  {
  public:
    ImproperRotation(const Vector3D& rotation_vector);
  };

  class Inversion : public SymmetryElement
  {
  public:
    Inversion()
      : SymmetryElement(SymmetryElement::Inversion, Matrix3D{
        -1.0, 0.0, 0.0,
        0.0, -1.0, 0.0,
        0.0, 0.0, -1.0
      })
    {}
  };

}
 
#endif // LATTICEGEN_SYMMETRYELEMENTS_H_