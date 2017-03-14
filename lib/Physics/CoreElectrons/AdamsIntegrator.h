#ifndef LATTICEGEN_ADAMSINTEGRATOR_H
#define LATTICEGEN_ADAMSINTEGRATOR_H

#include <vector>
#include <stdexcept>
#include <memory>

#include <Core/ExponentialMesh.h>

using namespace Core;

namespace Physics
{
  namespace CoreElectrons
  {
    struct AdamsIntegratorConfig
    {
      AdamsIntegratorConfig(
        int adams_moulton_quadrature_order_ = default_adams_moulton_quadrature_order,
        int inward_asymptotic_expansion_order_ = default_inward_asymptotic_expansion_order,
        double inward_asymptotic_expansion_cutoff_ = default_inward_asymptotic_expansion_cutoff,
        int outward_quadrature_order_ = default_outward_quadrature_order,
        int outward_scheme_repetition_ = default_outward_scheme_repetition
      ) : adams_moulton_quadrature_order(adams_moulton_quadrature_order_),
          inward_asymptotic_expansion_order(inward_asymptotic_expansion_order_),
          inward_asymptotic_expansion_cutoff(inward_asymptotic_expansion_cutoff_),
          outward_quadrature_order(outward_quadrature_order_), outward_scheme_repetition(outward_scheme_repetition_)
      {
        if (adams_moulton_quadrature_order < 1 || adams_moulton_quadrature_order > 8)
          throw std::invalid_argument("adams_moulton_quadrature_order must be int the range [1,8]");
        if (inward_asymptotic_expansion_order < 1)
          throw std::invalid_argument("inward_asymptotic_expansion_order must be positive");
        if (inward_asymptotic_expansion_cutoff < 0.0)
          throw std::invalid_argument("inward_asymptotic_expansion_cutoff must be positive");
        if (outward_quadrature_order < 1)
          throw std::invalid_argument("outward_quadrature_order must be positive");
        if (outward_scheme_repetition < 1)
          throw std::invalid_argument("outward_scheme_repetition must be positive");
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
          throw std::invalid_argument("invalid energy in AdamsIntegrator, energy = " + std::to_string(energy_));
        if (l_ < 0)
          throw std::invalid_argument("l must be non-negative l = " + std::to_string(l));
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

#endif //LATTICEGEN_ADAMSINTEGRATOR_H
