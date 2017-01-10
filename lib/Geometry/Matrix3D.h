#ifndef LATTICEGEN_MATRIX3D_H
#define LATTICEGEN_MATRIX3D_H

#include <Core/ComparisonHelpers.h>

#include "Point3D.h"

using namespace Core;

#include <array>
#include <limits>
#include <iostream>
#include <iomanip>
#include <sstream>

namespace Geometry
{
  typedef std::array< std::array<double,3>, 3> Matrix3D;

  inline Matrix3D operator*(const double alpha, const Matrix3D& matrix)
  {
    return Matrix3D {
      matrix[0][0] * alpha, matrix[0][1] * alpha, matrix[0][2] * alpha,
      matrix[1][0] * alpha, matrix[1][1] * alpha, matrix[1][2] * alpha,
      matrix[2][0] * alpha, matrix[2][1] * alpha, matrix[2][2] * alpha
    };
  }

  inline Point3D operator*(const Matrix3D& matrix, const Vector3D& vector)
  {
    return Vector3D(
    matrix[0][0] * vector.x + matrix[0][1] * vector.y + matrix[0][2] * vector.z,
    matrix[1][0] * vector.x + matrix[1][1] * vector.y + matrix[1][2] * vector.z,
    matrix[2][0] * vector.x + matrix[2][1] * vector.y + matrix[2][2] * vector.z
    );
  }

  inline Matrix3D operator*(const Matrix3D& lhs, const Matrix3D& rhs)
  {
    return Matrix3D {
    lhs[0][0] * rhs[0][0] + lhs[0][1] * rhs[1][0] + lhs[0][2] * rhs[2][0],
    lhs[0][0] * rhs[0][1] + lhs[0][1] * rhs[1][1] + lhs[0][2] * rhs[2][1],
    lhs[0][0] * rhs[0][2] + lhs[0][1] * rhs[1][2] + lhs[0][2] * rhs[2][2],

    lhs[1][0] * rhs[0][0] + lhs[1][1] * rhs[1][0] + lhs[1][2] * rhs[2][0],
    lhs[1][0] * rhs[0][1] + lhs[1][1] * rhs[1][1] + lhs[1][2] * rhs[2][1],
    lhs[1][0] * rhs[0][2] + lhs[1][1] * rhs[1][2] + lhs[1][2] * rhs[2][2],

    lhs[2][0] * rhs[0][0] + lhs[2][1] * rhs[1][0] + lhs[2][2] * rhs[2][0],
    lhs[2][0] * rhs[0][1] + lhs[2][1] * rhs[1][1] + lhs[2][2] * rhs[2][1],
    lhs[2][0] * rhs[0][2] + lhs[2][1] * rhs[1][2] + lhs[2][2] * rhs[2][2]
    };
  }

  inline Matrix3D operator+(const Matrix3D& lhs, const Matrix3D& rhs)
  {
    return Matrix3D {
    lhs[0][0] + rhs[0][0],  lhs[0][1] + rhs[0][1],  lhs[0][2] + rhs[0][2],
    lhs[1][0] + rhs[1][0],  lhs[1][1] + rhs[1][1],  lhs[1][2] + rhs[1][2],
    lhs[2][0] + rhs[2][0],  lhs[2][1] + rhs[2][1],  lhs[2][2] + rhs[2][2]
    };
  }

  inline Matrix3D operator-(const Matrix3D& lhs, const Matrix3D& rhs)
  {
    return Matrix3D {
    lhs[0][0] - rhs[0][0],  lhs[0][1] - rhs[0][1],  lhs[0][2] - rhs[0][2],
    lhs[1][0] - rhs[1][0],  lhs[1][1] - rhs[1][1],  lhs[1][2] - rhs[1][2],
    lhs[2][0] - rhs[2][0],  lhs[2][1] - rhs[2][1],  lhs[2][2] - rhs[2][2]
    };
  }
}

namespace std
{
  using namespace Geometry;

  inline std::string to_string(const Geometry::Matrix3D& matrix)
  {
    std::stringstream ss;
    ss << std::scientific << std::setprecision(std::numeric_limits<double>::max_digits10)
       << "[" << matrix[0][0] << " , " << matrix[0][1] << " , " << matrix[0][2] << std::endl
       << matrix[1][0] << " , " << matrix[1][1] << " , " << matrix[1][2] << std::endl
       << matrix[2][0] << " , " << matrix[2][1] << " , " << matrix[2][2] << " ] " << std::endl;
    return ss.str();
  }
  inline ostream& operator<<(ostream& o, const Geometry::Matrix3D& matrix)
  {
    o << std::to_string(matrix);
    return o;
  }
}
#endif //LATTICEGEN_MATRIX3D_H
