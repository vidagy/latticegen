#include "RadialSchrodingerEquation.h"

using namespace Physics::CoreElectrons;

RadialSolution Physics::CoreElectrons::RadialSchrodingerEquation::solve(
  unsigned int n, unsigned int l, double energy)
{
  const auto &r = effective_charge.r->points; // safe because class always has a strong reference on it
  const auto &z = effective_charge.z;
  const auto required_number_of_nodes = n - l - 1;

  std::vector<double> R(r.size(), 0.0);
  std::vector<double> dR_dr(r.size(), 0.0);
  auto e_lower = 0.0;
  auto e_upper = 0.0;
  auto practical_infinity = r.size() - 1;

  auto number_of_iteration = 0u;
  while (number_of_iteration++ < max_iter) {
    /// above r_max the wave function is practically zero. we don't calculate in this region
    practical_infinity = get_practical_infinity(r, z, energy);
    /// classical_turning_point is the boundary between the inward and the outward integration
    auto classical_turning_point = get_classical_turning_point(r, z, energy, practical_infinity);

    /// integrating the Schrodinger Equation from r_max to classical_turning_point
    /// this function will overwrite the corresponding region of R and dR_dr
    integrate_inward(r, z, energy, l, practical_infinity, classical_turning_point, R, dR_dr);
    /// integrating the Schrodinger Equation from 0 to classical_turning_point
    /// this function will overwrite the corresponding region of R and dR_dr
    integrate_outward(r, z, energy, l, classical_turning_point, R, dR_dr);
    /// make R and dR_dr continuous
    match_solutions(classical_turning_point, R, dR_dr);

    /// get the number of R = 0 (not counting the origin).
    auto number_of_nodes = get_number_of_nodes(R);
    if (required_number_of_nodes != number_of_nodes) {
      /// if it is not the required then we need a bigger energy step.
      update_energy_coarse(energy, e_lower, e_upper);
    } else {
      /// if the node number is OK, then we fine-tune the energy so that the dR_dr
      /// becomes continuous as well
      auto norm = get_norm(r, R);
      bool is_ready = update_energy_fine(energy, e_lower, e_upper, norm);
      if (is_ready)
        /// if dR_dr is continuous we normalize and exit
        normalize_solution(R, dR_dr, norm);
      break;
      /// if not then let's iterate again
    }
  }

  if (number_of_iteration > max_iter)
    throw std::logic_error("RadialSchrodingerEquation::solve didn't converge. number_of_iteration="
                           + std::to_string(number_of_iteration) + " while max_iter="
                           + std::to_string(max_iter));

  return RadialSolution(n, l, effective_charge.r, R, dR_dr, energy, practical_infinity, number_of_iteration);
}

unsigned long
RadialSchrodingerEquation::get_practical_infinity(
  const std::vector<double> &r, const std::vector<double> &z, double energy)
{
  auto practical_infinity = r.size() - 1;
  static const double asymptotic_region = 40.0 * 40.0 / 2.0;

  while (
    (((energy * r[practical_infinity] + z[practical_infinity]) * r[practical_infinity] + asymptotic_region) < 0.0) &&
    (practical_infinity > 0)
    ) {
    --practical_infinity;
  }

  if (practical_infinity == 0)
    throw std::logic_error("in get_practical_infinity, practical_infinity == 0");

  return practical_infinity;
}

unsigned long
RadialSchrodingerEquation::get_classical_turning_point(
  const std::vector<double> &r, const std::vector<double> &z, double energy, unsigned long practical_infinity)
{
  auto classical_turning_point = practical_infinity;

  while (
    ((energy * r[classical_turning_point] + z[classical_turning_point]) < 0.0) &&
    (classical_turning_point > 0)
    ) {
    --classical_turning_point;
  }

  if (classical_turning_point == 0)
    throw std::logic_error("in get_classical_turning_pointy, classical_turning_point == 0");

  return classical_turning_point;
}

void RadialSchrodingerEquation::integrate_inward(
  const std::vector<double> &r, const std::vector<double> &z, double energy,
  unsigned int l, unsigned long r_max, unsigned long classical_turning_point,
  std::vector<double> &R, std::vector<double> &dR_dr)
{

}

void
RadialSchrodingerEquation::integrate_outward(
  const std::vector<double> &r, const std::vector<double> &z, double energy,
  unsigned int l, unsigned long classical_turning_point,
  std::vector<double> &R, std::vector<double> &dR_dr)
{

}

void RadialSchrodingerEquation::match_solutions(
  unsigned long classical_turning_point, std::vector<double> &R, std::vector<double> &dr
)
{

}

unsigned int RadialSchrodingerEquation::get_number_of_nodes(const std::vector<double> &R)
{
  return 0;
}

void RadialSchrodingerEquation::update_energy_coarse(double &energy, double &lower, double &upper)
{

}

double RadialSchrodingerEquation::get_norm(const std::vector<double> &r, const std::vector<double> &R)
{
  return 0;
}

bool RadialSchrodingerEquation::update_energy_fine(double &energy, double &lower, double &upper, double norm)
{
  return false;
}

void RadialSchrodingerEquation::normalize_solution(std::vector<double> &R, std::vector<double> &dr, double norm)
{

}
