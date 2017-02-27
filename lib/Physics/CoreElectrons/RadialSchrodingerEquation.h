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
        double E_, unsigned long r_max_, unsigned int number_of_iteration_
      )
        : n(n_), l(l_),
          r(r_), R(R_), dR_dr(dR_dr_),
          E(E_), r_max(r_max_), number_of_iteration(number_of_iteration_)
      {
        if (R.size() != r->points.size())
          throw std::invalid_argument("in EffectiveCharge R.size()=" + std::to_string(R.size()) +
                                      " while r->points.size()=" + std::to_string(r->points.size()));
        if (R.size() != dR_dr.size())
          throw std::invalid_argument("in EffectiveCharge R.size()=" + std::to_string(R.size()) +
                                      " while dR_dr.size()=" + std::to_string(dR_dr.size()));
        if (r_max_ >= r->points.size())
          throw std::invalid_argument("in EffectiveCharge r_max=" + std::to_string(r_max) +
                                      " while r->points.size()=" + std::to_string(r->points.size()));
      }

      const unsigned int n;                           /// principal quantum number
      const unsigned int l;                           /// angular momentum quantum number
      const std::shared_ptr<const ExponentialMesh> r; /// radial points
      const std::vector<double> R;                    /// radial component of the wavefunction
      const std::vector<double> dR_dr;                /// derivative of the radial component of the wavefunction
      const double E;                                 /// eigenvalue
      const unsigned long r_max;                      /// practical infinity for the wavefunction (index in r)
      const unsigned int number_of_iteration;         /// number of iteration in the solver
    };

    class RadialSchrodingerEquation
    {
    public:
      RadialSchrodingerEquation(const EffectiveCharge &effective_charge_, unsigned int max_iter_ = 50)
        : effective_charge(effective_charge_), max_iter(max_iter_) {}

      RadialSolution solve(unsigned int n, unsigned int l, double energy_guess);

    private:
      const EffectiveCharge effective_charge;
      const unsigned int max_iter;

      unsigned long get_practical_infinity(const std::vector<double> &r, const std::vector<double> &z, double energy);

      unsigned long get_classical_turning_point(
        const std::vector<double> &r, const std::vector<double> &z, double energy, unsigned long practical_infinity);

      void integrate_inward(const std::vector<double> &r, const std::vector<double> &z, double energy,
                            unsigned int l, unsigned long r_max, unsigned long classical_turning_point,
                            std::vector<double> &R, std::vector<double> &dR_dr);

      void integrate_outward(const std::vector<double> &r, const std::vector<double> &z, double energy,
                             unsigned int l, unsigned long classical_turning_point,
                             std::vector<double> &R, std::vector<double> &dR_dr);

      void match_solutions(unsigned long classical_turning_point, std::vector<double> &R, std::vector<double> &dr);

      unsigned int get_number_of_nodes(const std::vector<double> &R);

      void update_energy_coarse(double &energy, double &lower, double &upper);

      double get_norm(const std::vector<double> &r, const std::vector<double> &R);

      bool update_energy_fine(double &energy, double &lower, double &upper, double norm);

      void normalize_solution(std::vector<double> &R, std::vector<double> &dr, double norm);

      friend class Accessor;
    };
  }
}

#endif //LATTICEGEN_RADIALSCHRODINGEREQUATION_H
