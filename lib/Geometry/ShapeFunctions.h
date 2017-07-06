#ifndef LATTICEGEN_SHAPEFUNCTIONS_H
#define LATTICEGEN_SHAPEFUNCTIONS_H

#include "UnitCell3D.h"
#include <Core/lm_vector.h>
#include <complex>

namespace Geometry
{
  using namespace Core;

  /// @brief theta and phi mesh configuration for spherical integration
  ///
  /// We use uniform mesh for theta. For phi, we pick phi_res number of points on the equator of the sphere and
  /// proportionally less as we approach the pole to keep the density of samples constant.
  struct ShapeFunctionsConfig
  {
    ShapeFunctionsConfig(
      int theta_res_ = default_theta_res, int phi_res_ = default_phi_res,
      unsigned int bracketing_max_iter_ = default_bracketing_max_iter
    )
      : theta_res(theta_res_), phi_res(phi_res_), bracketing_max_iter(bracketing_max_iter_)
    {
      if (theta_res < 2)
        THROW_INVALID_ARGUMENT("theta_res must be grater than 1, now theta_res = " + std::to_string(theta_res));
      if (phi_res < 2)
        THROW_INVALID_ARGUMENT("phi_res must be grater than 1, now theta_res = " + std::to_string(phi_res));
      if (bracketing_max_iter < 2u)
        THROW_INVALID_ARGUMENT("bracketing_max_iter must be grater than 1, now bracketing_max_iter = "
                               + std::to_string(bracketing_max_iter));
    }

    const static int default_theta_res = 200;
    const static int default_phi_res = 200;
    const static unsigned int default_bracketing_max_iter = 100;

    const int theta_res;
    const int phi_res;
    const unsigned int bracketing_max_iter;
  };

  class ShapeFunctions
  {
  public:
    ShapeFunctions(
      const UnitCell3D &unit_cell_, unsigned int l_max, const std::vector<double> &r_points_,
      const ShapeFunctionsConfig &config_ = ShapeFunctionsConfig());

    const UnitCell3D unit_cell;
    const lm_vector<std::vector<std::complex<double>>> shape_functions;
    const std::vector<double> r_points;
  };
}


#endif //LATTICEGEN_SHAPEFUNCTIONS_H
