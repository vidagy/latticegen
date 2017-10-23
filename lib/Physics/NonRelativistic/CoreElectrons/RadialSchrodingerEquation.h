#ifndef LATTICEGEN_RADIALSCHRODINGEREQUATION_H
#define LATTICEGEN_RADIALSCHRODINGEREQUATION_H

#include <vector>
#include <memory>

#include <Physics/Common/EffectiveCharge.h>
#include "AdamsIntegrator.h"

using namespace Core;
using namespace Physics::Common;

namespace Physics
{
  namespace NonRelativistic
  {
    namespace CoreElectrons
    {
      ///////////////////////////////////////////////////////////////////////////
      // The numeric solution of the Schrodinger Equation is based on chapter  //
      //          2.3 of Atomic Structure Theory by Walter R. Johnson          //
      ///////////////////////////////////////////////////////////////////////////

      class RadialSolution
      {
      public:
        RadialSolution(
          int n_, int l_,
          const std::shared_ptr<const ExponentialMesh> &r_, std::vector<double> R_, std::vector<double> dR_dr_,
          double E_, int practical_infinity_, int number_of_iteration_
        )
          : n(n_), l(l_),
            r(r_), R(R_), dR_dr(dR_dr_),
            E(E_), practical_infinity(practical_infinity_), number_of_iteration(number_of_iteration_)
        {
          if (R.size() != r->get_points().size())
            THROW_INVALID_ARGUMENT("in RadialSolution R.size()=" + std::to_string(R.size()) +
                                   " while r->get_points().size()=" + std::to_string(r->get_points().size()));
          if (R.size() != dR_dr.size())
            THROW_INVALID_ARGUMENT("in RadialSolution R.size()=" + std::to_string(R.size()) +
                                   " while dR_dr.size()=" + std::to_string(dR_dr.size()));
          if (practical_infinity_ >= static_cast<int>(r->get_points().size()))
            THROW_INVALID_ARGUMENT("in RadialSolution practical_infinity=" + std::to_string(practical_infinity) +
                                   " while r->get_points().size()=" + std::to_string(r->get_points().size()));
        }

        const int n;                                    /// principal quantum number
        const int l;                                    /// angular momentum quantum number
        const std::shared_ptr<const ExponentialMesh> r; /// radial points
        const std::vector<double> R;                    /// radial component of the wavefunction
        const std::vector<double> dR_dr;                /// derivative of the radial component of the wavefunction
        const double E;                                 /// eigenvalue
        const int practical_infinity;                   /// practical infinity for the wavefunction (index in r)
        const int number_of_iteration;                  /// number of iteration in the solver
      };

      class RadialSchrodingerEquation
      {
      public:
        explicit RadialSchrodingerEquation(const EffectiveCharge &effective_charge_)
          : effective_charge(effective_charge_) {}

        constexpr static const double default_energy_tolerance = 1e-11;
        static const int default_max_iter = 50;

        RadialSolution solve(
          int n, int l, double energy_guess,
          const AdamsIntegratorConfig &adams_integrator_config = AdamsIntegratorConfig(),
          double energy_tolerance_ = default_energy_tolerance,
          int max_iter_ = default_max_iter
        ) const;


        const EffectiveCharge effective_charge;

      private:
        static int get_practical_infinity(
          const std::vector<double> &r, const std::vector<double> &z, double energy);

        static int get_classical_turning_point(
          const std::vector<double> &r, const std::vector<double> &z, double energy, int practical_infinity);

        static int get_number_of_nodes(const std::vector<double> &R, int practical_infinity);

        static double get_norm(
          const std::vector<double> &dr_di, double dx, const std::vector<double> &R, int practical_infinity);

        static void normalize_solution(
          std::vector<double> &R, std::vector<double> &dR_dr, double norm, int practical_infinity);

        friend class Accessor;
      };
    }
  }
}

#endif //LATTICEGEN_RADIALSCHRODINGEREQUATION_H
