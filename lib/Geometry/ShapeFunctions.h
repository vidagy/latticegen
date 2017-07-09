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

  struct ShapeFunctionsConfig
  {
    ShapeFunctionsConfig(
      LebedevQuadrature::Order lebedev_order_ = default_lebedev_order,
      unsigned int r_ws_bits_ = default_r_ws_bits,
      unsigned int bracketing_max_iter_ = default_bracketing_max_iter
    )
      : lebedev_order(lebedev_order_), r_ws_bits(r_ws_bits_), bracketing_max_iter(bracketing_max_iter_)
    {
      if (bracketing_max_iter < 2u)
        THROW_INVALID_ARGUMENT("bracketing_max_iter must be grater than 1, now bracketing_max_iter = "
                               + std::to_string(bracketing_max_iter));
    }

    const static LebedevQuadrature::Order default_lebedev_order = LebedevQuadrature::Order::LD5810;
    const static unsigned int default_r_ws_bits = 12u;
    const static unsigned int default_bracketing_max_iter = 100u;

    const LebedevQuadrature::Order lebedev_order;
    const unsigned int r_ws_bits;
    const unsigned int bracketing_max_iter;
  };

  class ShapeFunctions
  {
  public:
    // TODO we should implement interpolation from going to r_WS mesh to custom input r mesh
    ShapeFunctions(
      const UnitCell3D &unit_cell_, unsigned int l_max, const ShapeFunctionsConfig &config_ = ShapeFunctionsConfig());

    const UnitCell3D unit_cell;
    const lm_vector<std::vector<std::complex<double>>> shape_functions;
  };
}


#endif //LATTICEGEN_SHAPEFUNCTIONS_H
