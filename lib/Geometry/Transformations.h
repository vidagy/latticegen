#ifndef LATTICEGEN_TRANSFORMATIONS_H_
#define LATTICEGEN_TRANSFORMATIONS_H_

#include <Core/ComparisonHelpers.h>
#include <Core/Exceptions.h>
#include <Core/Point3D.h>
#include <Core/Matrix3D.h>
#include <array>

using namespace Core;

namespace Geometry
{
  class Transformation;

  Transformation operator*(const Transformation &lhs, const Transformation &rhs);

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

    Eigen::Vector3d operator()(const Point3DCRef &vector) const
    {
      return transformation_matrix * vector;
    }

    Type type;
    Eigen::Matrix3d transformation_matrix;

  protected:
    Transformation(const Type type_, Eigen::Matrix3d matrix_)
      : type(type_), transformation_matrix(std::move(matrix_))
    {}

    friend Transformation Geometry::operator*(const Transformation &, const Transformation &);
  };

  inline Eigen::Vector3d operator*(Transformation &symmetry_element, const Point3DCRef &vector)
  {
    return symmetry_element.transformation_matrix * vector;
  }
  inline bool operator==(const Transformation& lhs, const Transformation& rhs)
  {
    // TODO comparison cannot shortcut on the Transformation::Type, since Rotation(0,0,0) is equal to Identity()
    return lhs.transformation_matrix.isApprox(rhs.transformation_matrix);
  }

  class Identity : public Transformation
  {
  public:
    Identity();
  };

  class Rotation : public Transformation
  {
  public:
    explicit Rotation(const Point3DCRef &rotation_vector);

    /// create from euler angles using the z-y-z convention
    /// R = R_z(gamma) * R_z(beta) * R_z(alpha)
    Rotation(double alpha, double beta, double gamma);

    struct EulerAngles
    {
      EulerAngles(double alpha, double beta, double gamma)
        : alpha(alpha), beta(beta), gamma(gamma) {}

      double alpha, beta, gamma;
    };

    EulerAngles get_euler_angles() const;
  };


  class Reflection : public Transformation
  {
  public:
    explicit Reflection(const Point3DCRef &reflection_plane);
  };


  class ImproperRotation : public Transformation
  {
  public:
    explicit ImproperRotation(const Point3DCRef &rotation_vector);
  };

  class Inversion : public Transformation
  {
  public:
    Inversion();
  };
}

namespace std
{
  using namespace Geometry;

  inline std::string to_string(const Geometry::Transformation::Type type)
  {
    std::stringstream ss;
    switch (type) {
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
        THROW_LOGIC_ERROR("unhandled transformation type " + std::to_string(static_cast<int>(type)));
    }
    return ss.str();
  }

  inline std::string to_string(const Geometry::Transformation &transformation) {
    std::stringstream ss;
    ss << std::to_string(transformation.type);
    ss << transformation.transformation_matrix;
    return ss.str();
  }

  inline ostream &operator<<(ostream &o, const Geometry::Transformation &transformation)
  {
    o << std::to_string(transformation);
    return o;
  }
}
 
#endif // LATTICEGEN_TRANSFORMATIONS_H_
