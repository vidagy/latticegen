#ifndef LATTICEGEN_MATRIX3D_H
#define LATTICEGEN_MATRIX3D_H

#include <Core/ComparisonHelpers.h>

LATTICEGEN_MUTE_BEGIN
LATTICEGEN_MUTE_EIGEN
#include <Eigen/Dense>

LATTICEGEN_MUTE_END

#include <limits>
#include <iostream>
#include <iomanip>
#include <sstream>

namespace Core {
  typedef Eigen::Ref<Eigen::Matrix3d> Matrix3dRef;
  typedef Eigen::Ref<const Eigen::Matrix3d> Matrix3dCRef;
}

#endif //LATTICEGEN_MATRIX3D_H
