#include "Transformations.h"
#include <Math/CommonFunctions.h>

using namespace Core;
using namespace Geometry;

namespace {
  const Matrix3D identity_matrix = {{
                                             {{1.0, 0.0, 0.0}},
                                             {{0.0, 1.0, 0.0}},
                                             {{0.0, 0.0, 1.0}}
                                           }};
}

Identity::Identity()
  : Transformation(Transformation::Identity, identity_matrix) {}

namespace {
  Matrix3D get_rotation_matrix(const Vector3D &rotation_vector) {
    const double angle = rotation_vector.length();

    if (nearlyZero(angle))
      THROW_INVALID_ARGUMENT("Rotation can not be created from null vector");

    const Vector3D axis = rotation_vector * (1.0 / angle);

    Matrix3D result = {{
                         {{cos(angle) + pow(axis.x, 2) * (1.0 - cos(angle)),
                            axis.x * axis.y * (1.0 - cos(angle)) - axis.z * sin(angle),
                            axis.x * axis.z * (1.0 - cos(angle)) + axis.y * sin(angle)}},

                         {{axis.y * axis.x * (1.0 - cos(angle)) + axis.z * sin(angle),
                            cos(angle) + pow(axis.y, 2) * (1.0 - cos(angle)),
                            axis.y * axis.z * (1.0 - cos(angle)) - axis.x * sin(angle)}},

                         {{axis.z * axis.x * (1.0 - cos(angle)) - axis.y * sin(angle),
                            axis.z * axis.y * (1.0 - cos(angle)) + axis.x * sin(angle),
                            cos(angle) + pow(axis.z, 2) * (1.0 - cos(angle))}}
                       }};

    return result;
  }
}

Rotation::Rotation(const Vector3D &rotation_vector)
  : Transformation(Transformation::Rotation, get_rotation_matrix(rotation_vector)) {}

namespace {
  Matrix3D get_rotation_from_euler(double alpha, double beta, double gamma) {
    const auto c1 = cos(alpha);
    const auto c2 = cos(beta);
    const auto c3 = cos(gamma);

    const auto s1 = sin(alpha);
    const auto s2 = sin(beta);
    const auto s3 = sin(gamma);

    return Matrix3D {{
                       {{c3 * c2 * c1 - s3 * s1, -c1 * s3 - c3 * c2 * s1, c3 * s2}},
                       {{c3 * s1 + c2 * c1 * s3, c3 * c1 - c2 * s3 * s1, s3 * s2}},
                       {{-c1 * s2, s2 * s1, c2}}
                     }};
  }
}

Rotation::Rotation(double alpha, double beta, double gamma)
  : Transformation(Transformation::Rotation, get_rotation_from_euler(alpha, beta, gamma)) {}

Rotation::EulerAngles Rotation::get_euler_angles() const {
  const auto &m = transformation_matrix;
  if (equalsWithTolerance(1.0, fabs(m[2][2]))) { // sin(beta) == 0 <==> cos(beta) == (+/-)1.0
    double beta;
    if (m[2][2] > 0.0)
      beta = 0;
    else
      beta = pi;

    // we will only calculate alpha and have gamma as zero.
    if (nearlyZero(m[0][0])) { // sin(beta) == 0; cos(alpha) == 0 <==> sin(alpha) == (+/-)1.0
      if (m[1][0] > 0.0)
        return EulerAngles(pi / 2.0, beta, 0.0);
      else
        return EulerAngles(-pi / 2.0, beta, 0.0);
    } else { // sin(beta) == 0; cos(alpha) != 0
      double alpha = atan(m[1][0] / m[1][1]);
      if (cos(alpha) * m[1][1] < 0) // they have the same sign
        alpha += pi;
      return EulerAngles(alpha, beta, 0.0);
    }
  } else { // sin(beta) != 0 <==> cos(beta) != (+/-)1.0
    double beta = acos(m[2][2]); // plus minus acos but we have the freedom to pick one, so we pick the positive
    double alpha, gamma;

    if (nearlyZero(m[2][0])) { // cos(alpha) == 0  => alpha= (+/-)pi/2.0
      alpha = (m[2][1] / sin(beta) > 0.0 ? 1.0 : -1.0) * pi / 2.0;
    } else { // cos(alpha) != 0
      alpha = atan(-m[2][1] / m[2][0]);
      auto cos_alpha_flip = m[2][0] / cos(alpha) * sin(beta) > 0.0;
      auto sin_alpha_flip = !nearlyZero(sin(alpha)) && m[2][1] / sin(alpha) * sin(beta) < 0.0;
      if (sin_alpha_flip && cos_alpha_flip) {
        alpha += pi;
      } else if (sin_alpha_flip) {
        alpha *= -1.0;
      } else if (cos_alpha_flip) {
        alpha = pi - alpha;
      }
    }
    if (nearlyZero(m[0][2])) { // cos(gamma) == 0  => gamma= (+/-)pi/2.0
      gamma = (m[1][2] / sin(beta) >= 0.0 ? 1.0 : -1.0) * pi / 2.0;
    } else { // cos(gamma) != 0
      gamma = atan(m[1][2] / m[0][2]);
      auto cos_gamma_flip = m[0][2] / cos(gamma) * sin(beta) < 0.0;
      auto sin_gamma_flip = !nearlyZero(sin(gamma)) && m[1][2] / sin(gamma) * sin(beta) < 0.0;
      if (sin_gamma_flip && cos_gamma_flip) {
        gamma += pi;
      } else if (sin_gamma_flip) {
        gamma *= -1.0;
      } else if (cos_gamma_flip) {
        gamma = pi - gamma;
      }
    }

    return EulerAngles(alpha, beta, gamma);
  }
}

namespace {
  Matrix3D get_diad(const Vector3D &vector) {
    return Matrix3D {{
                       {{vector.x * vector.x, vector.x * vector.y, vector.x * vector.z}},
                       {{vector.y * vector.x, vector.y * vector.y, vector.y * vector.z}},
                       {{vector.z * vector.x, vector.z * vector.y, vector.z * vector.z}}
                     }};
  }

  Matrix3D get_reflection_matrix(const Vector3D &vector) {
    const Matrix3D unit_matrix = Identity().transformation_matrix;
    return unit_matrix - 2.0 * get_diad(vector);
  }
}

Reflection::Reflection(const Vector3D &reflection_plane)
  : Transformation(Transformation::Reflection, get_reflection_matrix(reflection_plane)) {
  if (!equalsWithTolerance(reflection_plane.length(), 1.0))
    THROW_INVALID_ARGUMENT("Non-unit vector on input of Reflection: " + std::to_string(reflection_plane));
}

ImproperRotation::ImproperRotation(const Vector3D &rotation_vector)
  : Transformation(Transformation::ImproperRotation,
                   Geometry::Reflection(rotation_vector / rotation_vector.length()).transformation_matrix *
                   Geometry::Rotation(rotation_vector).transformation_matrix
) {}

namespace {
  const Matrix3D inversion_matrix = {{
                                              {{-1.0, 0.0, 0.0}},
                                              {{0.0, -1.0, 0.0}},
                                              {{0.0, 0.0, -1.0}}
                                            }};
}

Inversion::Inversion()
  : Transformation(Transformation::Inversion, inversion_matrix) {}


Transformation Geometry::operator*(const Transformation &lhs, const Transformation &rhs) {
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