#ifndef LATTICEGEN_SHAPEFUNCTIONS_H
#define LATTICEGEN_SHAPEFUNCTIONS_H

#include "UnitCell3D.h"
#include <Core/lm_vector.h>
#include <complex>
#include <Core/RadialMesh.h>

namespace Geometry
{
  using namespace Core;

  struct ShapeFunctionsConfig
  {
    ShapeFunctionsConfig(
      unsigned int quadrature_order_ = default_qadrature_order,
      unsigned int r_ws_bits_ = default_r_ws_bits,
      unsigned int bracketing_max_iter_ = default_bracketing_max_iter
    )
      : quadrature_order(quadrature_order_), r_ws_bits(r_ws_bits_), bracketing_max_iter(bracketing_max_iter_)
    {
      if (bracketing_max_iter < 2u)
        THROW_INVALID_ARGUMENT("bracketing_max_iter must be grater than 1, now bracketing_max_iter = "
                               + std::to_string(bracketing_max_iter));
    }

    const static unsigned int default_qadrature_order = 8u;
    const static unsigned int default_r_ws_bits = 12u;
    const static unsigned int default_bracketing_max_iter = 100u;

    const unsigned int quadrature_order;
    const unsigned int r_ws_bits;
    const unsigned int bracketing_max_iter;
  };

  class ShapeFunctions
  {
  public:
    ShapeFunctions(
      const UnitCell3D &unit_cell_, unsigned int l_max, const std::shared_ptr<RadialMesh> &mesh_,
      const ShapeFunctionsConfig &config_ = ShapeFunctionsConfig()
    );

    const UnitCell3D unit_cell;
    const unsigned int l_max;
    const std::shared_ptr<RadialMesh> mesh;
    const lm_vector<std::vector<std::complex<double>>> shape_functions;
  };
}


#endif //LATTICEGEN_SHAPEFUNCTIONS_H
