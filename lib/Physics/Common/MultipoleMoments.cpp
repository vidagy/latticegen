#include "MultipoleMoments.h"
#include <Math/CommonFunctions.h>
#include <Math/Integrator.h>

using namespace Math;
using namespace Physics::Common;

namespace
{
  std::complex<double> integral(
    unsigned int l,
    const ExponentialMesh &mesh,
    std::vector<std::complex<double>> shape_truncated_charge_density
  )
  {
    const auto &r = mesh.get_points();
    for (auto i = 0u; i < r.size(); ++i) {
      shape_truncated_charge_density[i] *= shape_truncated_charge_density[i] * Math::pow(r[i], 2 + l);
    }

    return IntegratorExponential::simpson_alt(shape_truncated_charge_density, mesh);
  }
}

MultipoleMoments::MultipoleMoments(
  const ShapeTruncatedChargeDensity &shape_truncated_charge_density
) : moments(shape_truncated_charge_density.density.l_max)
{
  const auto &mesh = *shape_truncated_charge_density.mesh;

  for (auto l = 0u; l <= shape_truncated_charge_density.density.l_max; ++l) {
    for (auto m = 0; m <= ((int) l); ++m) {
      moments.at(l, m) =
        sqrt(4.0 * pi) / (2.0 * l + 1.0) * integral(l, mesh, shape_truncated_charge_density.density.at(l, m));
      if (m != 0)
        moments.at(l, -m) = sign(m) * std::conj(moments.at(l, m));
    }
  }
}
