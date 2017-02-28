#ifndef LATTICEGEN_ADAMSINTEGRATOR_H
#define LATTICEGEN_ADAMSINTEGRATOR_H

#include <vector>

namespace Physics
{
  namespace CoreElectrons
  {
    class AdamsIntegrator
    {
    public:
      AdamsIntegrator(
        const std::vector<double> &r_, const std::vector<double> &z_, double energy_, unsigned int l_,
        unsigned long practical_infinity_, unsigned long classical_turning_point_
      )
        : r(r_), z(z_), energy(energy_), l(l_), practical_infinity(practical_infinity_),
          classical_turning_point(classical_turning_point_) {}

      void integrate(std::vector<double> &R, std::vector<double> &dR_dr);

      double get_old_dR_dr_scaled() const { return old_dR_dr_scaled; }

    private:
      void adams_integrator(std::vector<double> &R, std::vector<double> &dR_dr) const;

      void integrate_inward(std::vector<double> &R, std::vector<double> &dR_dr) const;

      void integrate_outward(std::vector<double> &R, std::vector<double> &dR_dr) const;

      void match_solutions(
        std::vector<double> &R, std::vector<double> &dR_dr, double old_R_at_ctp, double old_dR_dr_at_ctp
      ) const;

      const std::vector<double> &r;
      const std::vector<double> &z;
      const double energy;
      const unsigned int l;
      const unsigned long practical_infinity;
      const unsigned long classical_turning_point;

      double old_dR_dr_scaled = 0.0;
    };
  }
}

#endif //LATTICEGEN_ADAMSINTEGRATOR_H
