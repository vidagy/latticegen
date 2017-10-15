#ifndef LATTICEGEN_CHARGEDENSITY_H
#define LATTICEGEN_CHARGEDENSITY_H

#include <boost/shared_ptr.hpp>
#include <Core/RadialMesh.h>
#include <Core/lm_vector.h>
#include <Math/ClebschGordan.h>
#include <Geometry/ShapeFunctions.h>

namespace Physics
{
  namespace Common
  {
    using namespace Core;
    using namespace Math;
    using namespace Geometry;

    struct ChargeDensity
    {
      ChargeDensity(
        const std::shared_ptr<ExponentialMesh> &mesh_, const lm_vector<std::vector<std::complex<double>>> &density_
      ) : mesh(mesh_), density(density_)
      {
        for (auto l = 0u; l <= density_.l_max; ++l) {
          for (auto m = -((int) l); m <= ((int) l); ++m) {
            if (mesh->get_points().size() != density.at(l, m).size())
              THROW_LOGIC_ERROR(
                "mesh has " + std::to_string(mesh->get_points().size()) +
                " while density at l = " + std::to_string(l) + " m = " + std::to_string(m) +
                " has " + std::to_string(density.at(l, m).size())
              );
          }
        }
      }

      std::shared_ptr<ExponentialMesh> mesh;
      lm_vector<std::vector<std::complex<double>>> density;
    };

    /// @brief as defined in Zabloudil et al (19.15)
    struct ShapeTruncatedChargeDensity
    {
      ShapeTruncatedChargeDensity(
        const ChargeDensity &charge_density, const ShapeFunctions &shape_functions
      ) : mesh(charge_density.mesh), density(charge_density.density.l_max)
      {
        if (!std::dynamic_pointer_cast<ExponentialMesh>(shape_functions.mesh))
          THROW_INVALID_ARGUMENT("shape functions must be defined over an exponential mesh");

        if (charge_density.mesh != shape_functions.mesh) { // pointer equality
          if (charge_density.mesh->get_points().size() != shape_functions.mesh->get_points().size())
            THROW_INVALID_ARGUMENT(
              " charge_density has to be created with the same radial mesh for ChargeDensity and shape_functions. "
                "Now size is different: "
                "for charge_density " + std::to_string(charge_density.mesh->get_points().size()) +
              " and for shape_functions " + std::to_string(shape_functions.mesh->get_points().size()));
          for (
            auto it1 = charge_density.mesh->get_points().cbegin(), it2 = shape_functions.mesh->get_points().cbegin();
            it1 != charge_density.mesh->get_points().cend() && it2 != shape_functions.mesh->get_points().cend();
            ++it1, ++it2
            ) {
            if (!equalsWithTolerance(*it1, *it2))
              THROW_INVALID_ARGUMENT(
                " charge_density has to be created with the same radial mesh for ChargeDensity and shape_functions. "
                  "Now elements are different, first diff is : " + std::to_string(fabs(*it1 - *it2))
              );
          }
        }

        if (shape_functions.l_max < charge_density.density.l_max)
          THROW_INVALID_ARGUMENT(
            "shape_functions.l_max = " + std::to_string(shape_functions.l_max) +
            " is less than charge_density.l_max = " + std::to_string(charge_density.density.l_max)
          );

        for (auto l = 0u; l <= density.l_max; ++l) {
          for (auto m = 0; m <= ((int) l); ++m) {
            density.at(l, m) = std::vector<std::complex<double>>(charge_density.mesh->get_points().size(), 0.0);
          }
        }

        for (auto l = 0u; l <= charge_density.density.l_max; ++l) {
          for (auto m = 0; m <= ((int) l); ++m) {
            auto &d = density.at(l, m);
            for (auto lp = 0u; lp <= charge_density.density.l_max; ++lp) {
              for (auto mp = -((int) lp); mp <= ((int) lp); ++mp) {
                for (auto lpp = 0u; lpp <= shape_functions.l_max; ++lpp) {
                  for (auto mpp = -((int) lpp); mpp <= ((int) lpp); ++mpp) {
                    auto cg = ClebschGordan::calculate(l, m, lpp, mpp, lp, mp);
                    if (!nearlyZero(cg))
                      for (auto i = 0; i < charge_density.mesh->get_points().size(); ++i)
                        d[i] +=
                          cg * charge_density.density.at(lp, mp)[i] * shape_functions.shape_functions.at(lpp, mpp)[i];
                  }
                }
              }
            }
          }
        }
      }

      std::shared_ptr<ExponentialMesh> mesh;
      lm_vector<std::vector<std::complex<double>>> density;
    };
  }
}

#endif //LATTICEGEN_CHARGEDENSITY_H
