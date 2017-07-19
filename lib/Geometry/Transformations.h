#ifndef LATTICEGEN_TRANSFORMATIONS_H_
#define LATTICEGEN_TRANSFORMATIONS_H_

#include <Core/ComparisonHelpers.h>
#include <Core/Exceptions.h>
#include <Core/Point3D.h>
#include "Matrix3D.h"
#include <array>

using namespace Core;

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

    friend Transformation Geometry::operator*(const Transformation &, const Transformation &);
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

  inline Transformation operator*(const Transformation &lhs, const Transformation &rhs)
  {
    // identity shortcuts
    if (lhs.type == Transformation::Identity)
      return rhs;
    if (rhs.type == Transformation::Identity)
      return lhs;

    auto res = lhs.transformation_matrix * rhs.transformation_matrix;
    auto tr = trace(res);
    if (equalsWithTolerance(tr, 3.0))
      return Geometry::Identity();
    if (equalsWithTolerance(tr, -3.0))
      return Geometry::Inversion();

    auto det = determinant(res);
    if (equalsWithTolerance(det, 1.0))
      return Transformation(Transformation::Rotation, res);
    else if (equalsWithTolerance(det, -1.0)) {
      if (equalsWithTolerance(tr, 1.0))
        return Transformation(Transformation::Reflection, res);
      else
        return Transformation(Transformation::ImproperRotation, res);
    } else {
      THROW_LOGIC_ERROR(
        "invalid Transformation types: lhs = " + std::to_string(lhs.type) + " rhs = " + std::to_string(rhs.type) +
        " determinant is not +- 1.0 but " + std::to_string(det)
      );
    }
  }
}

namespace std
{
  using namespace Geometry;

  inline std::string to_string(const Geometry::Transformation &transformation)
  {
    std::stringstream ss;
    switch (transformation.type) {
      case Geometry::Transformation::Identity:
        ss << "Identity";
        break;
      case Geometry::Transformation::Rotation:
        ss << "Rotation";
        break;
      case Geometry::Transformation::Reflection:
        ss << "Reflection";
        break;
      case Geometry::Transformation::ImproperRotation:
        ss << "ImproperRotation";
        break;
      case Geometry::Transformation::Inversion:
        ss << "Inversion";
        break;
      default:
        THROW_LOGIC_ERROR("unhandled transformation type " + std::to_string(transformation.type));
    }
    ss << std::to_string(transformation.transformation_matrix);
    return ss.str();
  }

  inline ostream &operator<<(ostream &o, const Geometry::Transformation &transformation)
  {
    o << std::to_string(transformation);
    return o;
  }
}
 
#endif // LATTICEGEN_TRANSFORMATIONS_H_
