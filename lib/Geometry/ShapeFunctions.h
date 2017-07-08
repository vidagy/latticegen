#ifndef LATTICEGEN_SHAPEFUNCTIONS_H
#define LATTICEGEN_SHAPEFUNCTIONS_H

#include "UnitCell3D.h"
#include <Math/LebedevQuadrature.h>
#include <Core/lm_vector.h>
#include <complex>

namespace Geometry
{
  using namespace Core;
  using namespace Math;

  /// @brief theta and phi mesh configuration for spherical integration
  ///
  /// We use uniform mesh for theta. For phi, we pick phi_res number of points on the equator of the sphere and
  /// proportionally less as we approach the pole to keep the density of samples constant.
  struct ShapeFunctionsConfig
  {
    ShapeFunctionsConfig(
      LebedevQuadrature::Order lebedev_order_ = default_lebedev_order,
      unsigned int bracketing_max_iter_ = default_bracketing_max_iter
    )
      : lebedev_order(lebedev_order_), bracketing_max_iter(bracketing_max_iter_)
    {
      if (bracketing_max_iter < 2u)
        THROW_INVALID_ARGUMENT("bracketing_max_iter must be grater than 1, now bracketing_max_iter = "
                               + std::to_string(bracketing_max_iter));
    }

    const static LebedevQuadrature::Order default_lebedev_order = LebedevQuadrature::Order::LD5810;
    const static unsigned int default_bracketing_max_iter = 100;

    const LebedevQuadrature::Order lebedev_order;
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
