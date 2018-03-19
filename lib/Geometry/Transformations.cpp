#include "Transformations.h"
#include <Math/CommonFunctions.h>

using namespace Core;
using namespace Geometry;

namespace {
  const Eigen::Matrix3d identity_matrix = Eigen::Matrix3d::Identity();
}

Identity::Identity()
  : Transformation(Transformation::Identity, identity_matrix) {}

namespace {
  Eigen::Matrix3d get_rotation_matrix(const Point3DCRef &rotation_vector) {
    const double angle = rotation_vector.norm();

    if (nearlyZero(angle))
      THROW_INVALID_ARGUMENT("Rotation can not be created from null vector");

    const Point3D axis = (1.0 / angle) * rotation_vector;

    Eigen::Matrix3d result;
    result <<
           cos(angle) + pow(axis(0), 2) * (1.0 - cos(angle)),
      axis(0) * axis(1) * (1.0 - cos(angle)) - axis(2) * sin(angle),
      axis(0) * axis(2) * (1.0 - cos(angle)) + axis(1) * sin(angle),

      axis(1) * axis(0) * (1.0 - cos(angle)) + axis(2) * sin(angle),
      cos(angle) + pow(axis(1), 2) * (1.0 - cos(angle)),
      axis(1) * axis(2) * (1.0 - cos(angle)) - axis(0) * sin(angle),

      axis(2) * axis(0) * (1.0 - cos(angle)) - axis(1) * sin(angle),
      axis(2) * axis(1) * (1.0 - cos(angle)) + axis(0) * sin(angle),
      cos(angle) + pow(axis(2), 2) * (1.0 - cos(angle));

    return result;
  }
}

Rotation::Rotation(const Point3DCRef &rotation_vector)
  : Transformation(Transformation::Rotation, get_rotation_matrix(rotation_vector)) {}

namespace {
  Eigen::Matrix3d get_rotation_from_euler(double alpha, double beta, double gamma) {
    const auto c1 = cos(alpha);
    const auto c2 = cos(beta);
    const auto c3 = cos(gamma);

    const auto s1 = sin(alpha);
    const auto s2 = sin(beta);
    const auto s3 = sin(gamma);

    Eigen::Matrix3d res;

    res <<
        c3 * c2 * c1 - s3 * s1, -c1 * s3 - c3 * c2 * s1, c3 * s2,
      c3 * s1 + c2 * c1 * s3, c3 * c1 - c2 * s3 * s1, s3 * s2,
      -c1 * s2, s2 * s1, c2;

    return res;
  }
}

Rotation::Rotation(double alpha, double beta, double gamma)
  : Transformation(Transformation::Rotation, get_rotation_from_euler(alpha, beta, gamma)) {}

Rotation::EulerAngles Rotation::get_euler_angles() const {
  const auto &m = transformation_matrix;
  if (equalsWithTolerance(1.0, fabs(m(2, 2)))) { // sin(beta) == 0 <==> cos(beta) == (+/-)1.0
    double beta;
    if (m(2, 2) > 0.0)
      beta = 0;
    else
      beta = pi;

    // we will only calculate alpha and have gamma as zero.
    if (nearlyZero(m(0, 0))) { // sin(beta) == 0; cos(alpha) == 0 <==> sin(alpha) == (+/-)1.0
      if (m(1, 0) > 0.0)
        return EulerAngles(pi / 2.0, beta, 0.0);
      else
        return EulerAngles(-pi / 2.0, beta, 0.0);
    } else { // sin(beta) == 0; cos(alpha) != 0
      double alpha = atan(m(1, 0) / m(1, 1));
      if (cos(alpha) * m(1, 1) < 0) // they have the same sign
        alpha += pi;
      return EulerAngles(alpha, beta, 0.0);
    }
  } else { // sin(beta) != 0 <==> cos(beta) != (+/-)1.0
    double beta = acos(m(2, 2)); // plus minus acos but we have the freedom to pick one, so we pick the positive
    double alpha, gamma;

    if (nearlyZero(m(2, 0))) { // cos(alpha) == 0  => alpha= (+/-)pi/2.0
      alpha = (m(2, 1) / sin(beta) > 0.0 ? 1.0 : -1.0) * pi / 2.0;
    } else { // cos(alpha) != 0
      alpha = atan(-m(2, 1) / m(2, 0));
      auto cos_alpha_flip = m(2, 0) / cos(alpha) * sin(beta) > 0.0;
      auto sin_alpha_flip = !nearlyZero(sin(alpha)) && m(2, 1) / sin(alpha) * sin(beta) < 0.0;
      if (sin_alpha_flip && cos_alpha_flip) {
        alpha += pi;
      } else if (sin_alpha_flip) {
        alpha *= -1.0;
      } else if (cos_alpha_flip) {
        alpha = pi - alpha;
      }
    }
    if (nearlyZero(m(0, 2))) { // cos(gamma) == 0  => gamma= (+/-)pi/2.0
      gamma = (m(1, 2) / sin(beta) >= 0.0 ? 1.0 : -1.0) * pi / 2.0;
    } else { // cos(gamma) != 0
      gamma = atan(m(1, 2) / m(0, 2));
      auto cos_gamma_flip = m(0, 2) / cos(gamma) * sin(beta) < 0.0;
      auto sin_gamma_flip = !nearlyZero(sin(gamma)) && m(1, 2) / sin(gamma) * sin(beta) < 0.0;
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
  Eigen::Matrix3d get_diad(const Point3DCRef &vector) {
    return vector * (vector.transpose());
  }

  Eigen::Matrix3d get_reflection_matrix(const Point3D &vector) {
    const Eigen::Matrix3d unit_matrix = Identity().transformation_matrix;
    return unit_matrix - 2.0 * get_diad(vector);
  }
}

Reflection::Reflection(const Point3DCRef &reflection_plane)
  : Transformation(Transformation::Reflection, get_reflection_matrix(reflection_plane)) {
  if (!equalsWithTolerance(reflection_plane.norm(), 1.0))
    THROW_INVALID_ARGUMENT("Non-unit vector on input of Reflection: " + std::to_string(reflection_plane));
}

ImproperRotation::ImproperRotation(const Point3DCRef &rotation_vector)
  : Transformation(Transformation::ImproperRotation,
                   Geometry::Reflection(rotation_vector / rotation_vector.norm()).transformation_matrix *
                   Geometry::Rotation(rotation_vector).transformation_matrix
) {}

namespace {
  const Eigen::Matrix3d inversion_matrix = -1.0 * Eigen::Matrix3d::Identity();
}

Inversion::Inversion()
  : Transformation(Transformation::Inversion, inversion_matrix) {}


Transformation Geometry::operator*(const Transformation &lhs, const Transformation &rhs) {
  const auto tolerance_multiplier = 100.0;
  // identity shortcuts
  if (lhs.type == Transformation::Identity)
    return rhs;
  if (rhs.type == Transformation::Identity)
    return lhs;

  auto res = lhs.transformation_matrix * rhs.transformation_matrix;
  auto tr = res.trace();
  if (equalsWithTolerance(tr, 3.0, tolerance_multiplier * default_abs_tolerance))
    return Geometry::Identity();
  if (equalsWithTolerance(tr, -3.0, tolerance_multiplier * default_abs_tolerance))
    return Geometry::Inversion();

  auto det = res.determinant();
  if (equalsWithTolerance(det, 1.0, tolerance_multiplier * default_abs_tolerance))
    return Transformation(Transformation::Rotation, res);
  else if (equalsWithTolerance(det, -1.0, tolerance_multiplier * default_abs_tolerance)) {
    if (equalsWithTolerance(tr, 1.0, tolerance_multiplier * default_abs_tolerance))
      return Transformation(Transformation::Reflection, res);
    else
      return Transformation(Transformation::ImproperRotation, res);
  } else {
    THROW_LOGIC_ERROR(
      "invalid Transformation type: lhs = " + std::to_string(lhs.type) + " rhs = " + std::to_string(rhs.type) +
      " determinant is not +- 1.0 but " + std::to_string(det * 1e10)
    );
  }
}