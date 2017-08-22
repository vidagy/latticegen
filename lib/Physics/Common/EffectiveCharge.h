#ifndef LATTICEGEN_EFFECTIVECHARGE_H
#define LATTICEGEN_EFFECTIVECHARGE_H

#include <Core/RadialMesh.h>

using namespace Core;

namespace Physics
{
  namespace Common
  {
    /// @brief The generalized Coulomb potential in atomic units is V(r) = - Z(r)/r.
    /// @brief EffectiveCharge is this Z(r) function
    class EffectiveCharge
    {
    public:
      EffectiveCharge(const std::vector<double> &z_, const std::shared_ptr<const ExponentialMesh> &r_)
        : z(z_), r(r_)
      {
        if (z.size() != r->get_points().size())
          THROW_INVALID_ARGUMENT("in EffectiveCharge z.size()=" + std::to_string(z.size()) +
                                 " while r->points.size()=" + std::to_string(r->get_points().size()));
      }

      const std::vector<double> z;
      const std::shared_ptr<const ExponentialMesh> r;
    };
  }
}

#endif //LATTICEGEN_EFFECTIVECHARGE_H
