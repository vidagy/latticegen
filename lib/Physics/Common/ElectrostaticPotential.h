#ifndef LATTICEGEN_ELECTROSTATICPOTENTIAL_H
#define LATTICEGEN_ELECTROSTATICPOTENTIAL_H

#include <boost/shared_ptr.hpp>
#include <utility>
#include <Core/RadialMesh.h>
#include <Core/lm_vector.h>

namespace Physics
{
  namespace Common
  {
    using namespace Core;

    /// @brief Electrostatic potential in atomic units a.k.a. V(r)
    struct ElectrostaticPotential
    {
      ElectrostaticPotential(const Core::lm_vector<std::vector<std::complex<double>>> &v_,
                             std::shared_ptr<ExponentialMesh> mesh_
      ) : mesh(std::move(mesh_)), v(v_) {
        for (auto l = 0u; l <= v_.l_max; ++l) {
          for (auto m = -((int) l); m <= ((int) l); ++m) {
            if (mesh->get_points().size() != v.at(l, m).size())
              THROW_LOGIC_ERROR(
                "mesh has " + std::to_string(mesh->get_points().size()) +
                " while v at l = " + std::to_string(l) + " m = " + std::to_string(m) +
                " has " + std::to_string(v.at(l, m).size())
              );
          }
        }
      }

      std::shared_ptr<ExponentialMesh> mesh;
      lm_vector<std::vector<std::complex<double>>> v;
    };
  }
}

#endif //LATTICEGEN_ELECTROSTATICPOTENTIAL_H
