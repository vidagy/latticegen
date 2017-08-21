#ifndef LATTICEGEN_ADAMSINTEGRATOR_H
#define LATTICEGEN_ADAMSINTEGRATOR_H

#include <vector>
#include <memory>

#include <Core/RadialMesh.h>
#include <Core/Exceptions.h>

using namespace Core;

namespace Physics
{
  namespace NonRelativistic
  {
    namespace CoreElectrons
    {
      struct AdamsIntegratorConfig
      {
#define MAX_ADAMS_MOULTON_QUADRATURE_ORDER 16
#define MAX_INWARD_ASYMPTOTIC_EXPANSION_ORDER 32

        AdamsIntegratorConfig(
          int adams_moulton_quadrature_order_ = default_adams_moulton_quadrature_order,
          int inward_asymptotic_expansion_order_ = default_inward_asymptotic_expansion_order,
          double inward_asymptotic_expansion_cutoff_ = default_inward_asymptotic_expansion_cutoff,
          int outward_quadrature_order_ = default_outward_quadrature_order,
          int outward_scheme_repetition_ = default_outward_scheme_repetition
        ) : adams_moulton_quadrature_order(adams_moulton_quadrature_order_),
            inward_asymptotic_expansion_order(inward_asymptotic_expansion_order_),
            inward_asymptotic_expansion_cutoff(inward_asymptotic_expansion_cutoff_),
            outward_quadrature_order(outward_quadrature_order_),
            outward_scheme_repetition(outward_scheme_repetition_)
        {
          if (adams_moulton_quadrature_order < 1 || adams_moulton_quadrature_order > MAX_ADAMS_MOULTON_QUADRATURE_ORDER)
            THROW_INVALID_ARGUMENT("adams_moulton_quadrature_order must be in range [1,16]");
          if (inward_asymptotic_expansion_order < 1
              || inward_asymptotic_expansion_order > MAX_INWARD_ASYMPTOTIC_EXPANSION_ORDER)
            THROW_INVALID_ARGUMENT("inward_asymptotic_expansion_order must be positive");
          if (inward_asymptotic_expansion_cutoff < 0.0)
            THROW_INVALID_ARGUMENT("inward_asymptotic_expansion_cutoff must be positive");
          if (outward_quadrature_order < 1)
            THROW_INVALID_ARGUMENT("outward_quadrature_order must be positive");
          if (outward_scheme_repetition < 1)
            THROW_INVALID_ARGUMENT("outward_scheme_repetition must be positive");
        }

        static const int default_adams_moulton_quadrature_order = 8;
        static const int default_inward_asymptotic_expansion_order = 15;
        constexpr static const double default_inward_asymptotic_expansion_cutoff = 1e-11;
        static const int default_outward_quadrature_order = 8;
        static const int default_outward_scheme_repetition = 1;

        int adams_moulton_quadrature_order;
        int inward_asymptotic_expansion_order;
        double inward_asymptotic_expansion_cutoff;
        int outward_quadrature_order;
        int outward_scheme_repetition;
      };

      class AdamsIntegrator
      {
      public:
        AdamsIntegrator(
          const std::shared_ptr<const ExponentialMesh> &r_, const std::vector<double> &z_, double energy_,
          int l_, int practical_infinity_, int classical_turning_point_, const AdamsIntegratorConfig &config_
        )
          : r(r_), z(z_), energy(energy_), l(l_), practical_infinity(practical_infinity_),
            classical_turning_point(classical_turning_point_), config(config_)
        {
          if (energy_ > 0.0)
            THROW_INVALID_ARGUMENT("invalid energy in AdamsIntegrator, energy = " + std::to_string(energy_));
          if (l_ < 0)
            THROW_INVALID_ARGUMENT("l must be non-negative l = " + std::to_string(l));
        }

        void integrate(std::vector<double> &R, std::vector<double> &dR_dr);

        double get_old_dR_dr_scaled() const { return old_dR_dr_scaled; }

      private:
        void adams_moulton_method(std::vector<double> &R, std::vector<double> &dR_dr, int from, int to) const;

        void integrate_inward(std::vector<double> &R, std::vector<double> &dR_dr) const;

        void start_inward(std::vector<double> &R, std::vector<double> &dR_dr) const;

        void integrate_outward(std::vector<double> &R, std::vector<double> &dR_dr) const;

        void start_outward(std::vector<double> &R, std::vector<double> &dR_dr) const;

        void match_solutions(std::vector<double> &R, std::vector<double> &dR_dr, double old_R_at_ctp) const;

        static std::vector<double> get_adams_parameters(int quadrature);

        const std::shared_ptr<const ExponentialMesh> &r;
        const std::vector<double> &z;
        const double energy;
        const int l;
        const int practical_infinity;
        const int classical_turning_point;
        const AdamsIntegratorConfig config;

        double old_dR_dr_scaled = 0.0;

        friend class TestAccessor;
      };
    }
  }
}

#endif //LATTICEGEN_ADAMSINTEGRATOR_H
