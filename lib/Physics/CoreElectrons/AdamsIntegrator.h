#ifndef LATTICEGEN_ADAMSINTEGRATOR_H
#define LATTICEGEN_ADAMSINTEGRATOR_H

#include <vector>
#include <stdexcept>

#include <Core/ExponentialMesh.h>
#include <memory>

using namespace Core;

namespace Physics
{
  namespace CoreElectrons
  {
    class AdamsIntegrator
    {
    public:
      AdamsIntegrator(
        const std::shared_ptr<const ExponentialMesh> &r_, const std::vector<double> &z_, double energy_,
        unsigned int l_,
        unsigned long practical_infinity_, unsigned long classical_turning_point_,
        const unsigned int quadrature_,
        int number_of_iter_
      )
        : r(r_), z(z_), energy(energy_), l(l_), practical_infinity(practical_infinity_),
          classical_turning_point(classical_turning_point_), quadrature(quadrature_), number_of_iter(number_of_iter_)
      {
        if (!(quadrature_ > 0u && quadrature_ < 9u))
          throw std::invalid_argument(
            "invalid quadrature in AdamsIntegrator, quadrature = " + std::to_string(quadrature_));
        if (energy_ > 0.0)
          throw std::invalid_argument("invalid energy in AdamsIntegrator, energy = " + std::to_string(energy_));
      }

      void integrate(std::vector<double> &R, std::vector<double> &dR_dr);

      double get_old_dR_dr_scaled() const { return old_dR_dr_scaled; }

    private:
      void adams_moulton_method(std::vector<double> &R, std::vector<double> &dR_dr, unsigned long from,
                                unsigned long to) const;

      void integrate_inward(std::vector<double> &R, std::vector<double> &dR_dr) const;

      void start_inward(std::vector<double> &R, std::vector<double> &dR_dr) const;

      void integrate_outward(std::vector<double> &R, std::vector<double> &dR_dr) const;

      void start_outward(std::vector<double> &R, std::vector<double> &dR_dr) const;

      void match_solutions(std::vector<double> &R, std::vector<double> &dR_dr, double old_R_at_ctp) const;

      const std::shared_ptr<const ExponentialMesh> &r;
      const std::vector<double> &z;
      const double energy;
      const unsigned int l;
      const unsigned long practical_infinity;
      const unsigned long classical_turning_point;
      unsigned int quadrature;
      int number_of_iter;

      double old_dR_dr_scaled = 0.0;

      friend class TestAccessor;
    };
  }
}

#endif //LATTICEGEN_ADAMSINTEGRATOR_H
