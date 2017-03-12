#ifndef LATTICEGEN_RADIALSCHRODINGEREQUATION_H
#define LATTICEGEN_RADIALSCHRODINGEREQUATION_H

#include <vector>
#include <memory>
#include <stdexcept>

#include <Core/ExponentialMesh.h>

using namespace Core;

namespace Physics
{
  namespace CoreElectrons
  {
    ///////////////////////////////////////////////////////////////////////////
    // The numeric solution of the Schrodinger Equation is based on chapter  //
    //          2.3 of Atomic Structure Theory by Walter R. Johnson          // 
    ///////////////////////////////////////////////////////////////////////////

    /// @brief The generalized Coulomb potential in atomic units is V(r) = - Z(r)/r.
    /// @brief EffectiveCharge is this Z(r) function
    class EffectiveCharge
    {
    public:
      EffectiveCharge(const std::vector<double> &z_, const std::shared_ptr<const ExponentialMesh> &r_)
        : z(z_), r(r_)
      {
        if (z.size() != r->points.size())
          throw std::invalid_argument("in EffectiveCharge z.size()=" + std::to_string(z.size()) +
                                      " while r->points.size()=" + std::to_string(r->points.size()));
      }

      const std::vector<double> z;
      const std::shared_ptr<const ExponentialMesh> r;
    };

    class RadialSolution
    {
    public:
      RadialSolution(
        const unsigned int n_, const unsigned int l_,
        const std::shared_ptr<const ExponentialMesh> &r_, std::vector<double> R_, std::vector<double> dR_dr_,
        double E_, unsigned long practical_infinity_, unsigned int number_of_iteration_
      )
        : n(n_), l(l_),
          r(r_), R(R_), dR_dr(dR_dr_),
          E(E_), practical_infinity(practical_infinity_), number_of_iteration(number_of_iteration_)
      {
        if (R.size() != r->points.size())
          throw std::invalid_argument("in RadialSolution R.size()=" + std::to_string(R.size()) +
                                      " while r->points.size()=" + std::to_string(r->points.size()));
        if (R.size() != dR_dr.size())
          throw std::invalid_argument("in RadialSolution R.size()=" + std::to_string(R.size()) +
                                      " while dR_dr.size()=" + std::to_string(dR_dr.size()));
        if (practical_infinity_ >= r->points.size())
          throw std::invalid_argument("in RadialSolution practical_infinity=" + std::to_string(practical_infinity) +
                                      " while r->points.size()=" + std::to_string(r->points.size()));
      }

      const unsigned int n;                           /// principal quantum number
      const unsigned int l;                           /// angular momentum quantum number
      const std::shared_ptr<const ExponentialMesh> r; /// radial points
      const std::vector<double> R;                    /// radial component of the wavefunction
      const std::vector<double> dR_dr;                /// derivative of the radial component of the wavefunction
      const double E;                                 /// eigenvalue
      const unsigned long practical_infinity;         /// practical infinity for the wavefunction (index in r)
      const unsigned int number_of_iteration;         /// number of iteration in the solver
    };

    class RadialSchrodingerEquation
    {
    public:
      RadialSchrodingerEquation(
        const EffectiveCharge &effective_charge_, double energy_tolerance_ = 1e-11, unsigned int max_iter_ = 50
      )
        : effective_charge(effective_charge_), energy_tolerance(energy_tolerance_), max_iter(max_iter_) {}

      RadialSolution solve(unsigned int n, unsigned int l, double energy_guess, unsigned int quadrature) const;

      const EffectiveCharge effective_charge;
      const double energy_tolerance;
      const unsigned int max_iter;

    private:
      static unsigned long get_practical_infinity(
        const std::vector<double> &r, const std::vector<double> &z, double energy);

      static unsigned long get_classical_turning_point(
        const std::vector<double> &r, const std::vector<double> &z, double energy, unsigned long practical_infinity);

      static unsigned int get_number_of_nodes(const std::vector<double> &R, unsigned long practical_infinity);

      static double get_norm(
        const std::vector<double> &r, double dx, const std::vector<double> &R, unsigned long practical_infinity);

      static void normalize_solution(
        std::vector<double> &R, std::vector<double> &dR_dr, double norm, unsigned long practical_infinity);

      friend class Accessor;
    };
  }
}

#endif //LATTICEGEN_RADIALSCHRODINGEREQUATION_H
