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
    // identity
    if (lhs.type == Transformation::Identity)
      return rhs;
    if (rhs.type == Transformation::Identity)
      return lhs;
    // inversion
    if (rhs.type == Transformation::Inversion || lhs.type == Transformation::Inversion) {
      if (rhs.type == Transformation::Inversion && lhs.type == Transformation::Inversion)
        return Geometry::Identity();
      auto res = lhs.transformation_matrix * rhs.transformation_matrix;
      if (rhs.type == Transformation::Reflection || lhs.type == Transformation::Reflection)
        return Geometry::Transformation(Transformation::Rotation, res);
      if (rhs.type == Transformation::Rotation || lhs.type == Transformation::Rotation)
        return Geometry::Transformation(Transformation::ImproperRotation, res);
      if (rhs.type == Transformation::ImproperRotation || lhs.type == Transformation::ImproperRotation)
        return Geometry::Transformation(Transformation::Rotation, res);
    }
    // rotation
    if (rhs.type == Transformation::Rotation || lhs.type == Transformation::Rotation) {
      auto res = lhs.transformation_matrix * rhs.transformation_matrix;
      if (rhs.type == Transformation::Rotation && lhs.type == Transformation::Rotation) {
        if (equalsWithTolerance(trace(res), 3.0))
          return Geometry::Identity();
        else
          return Geometry::Transformation(Transformation::Rotation, res);
      }
      if (rhs.type == Transformation::Reflection || lhs.type == Transformation::Reflection) {
        return Geometry::Transformation(Transformation::ImproperRotation, res);
      }
      if (rhs.type == Transformation::ImproperRotation || lhs.type == Transformation::ImproperRotation) {
        if (equalsWithTolerance(trace(res), 1.0))
          return Geometry::Transformation(Transformation::Reflection, res);
        if (equalsWithTolerance(trace(res), -3.0))
          return Geometry::Inversion();
        else
          return Geometry::Transformation(Transformation::ImproperRotation, res);
      }
    }
    // reflection or improper rotation
    if ((rhs.type == Transformation::Reflection || rhs.type == Transformation::ImproperRotation) &&
        (lhs.type == Transformation::Reflection || lhs.type == Transformation::ImproperRotation)) {
      auto res = lhs.transformation_matrix * rhs.transformation_matrix;
      if (equalsWithTolerance(trace(res), 3.0))
        return Geometry::Identity();
      else
        return Geometry::Transformation(Transformation::Rotation, res);
    }

    THROW_LOGIC_ERROR(
      "invalid Transformation types: lhs = " + std::to_string(lhs.type) + " rhs = " + std::to_string(rhs.type)
    );
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
