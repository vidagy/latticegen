#include <Math/Integrator.h>
#include "RadialSchrodingerEquation.h"
#include "EnergyUpdate.h"
#include "AdamsIntegrator.h"

using namespace Physics::CoreElectrons;

RadialSolution Physics::CoreElectrons::RadialSchrodingerEquation::solve(
  unsigned int n, unsigned int l, double energy) const
{
  const auto &r = effective_charge.r->points; // safe because class always has a strong reference on it
  const auto dx = effective_charge.r->dx;
  const auto &z = effective_charge.z;
  const auto required_number_of_nodes = n - l - 1;
  auto update_energy = EnergyUpdate(required_number_of_nodes, energy_tolerance);

  std::vector<double> R(r.size(), 0.0);
  std::vector<double> dR_dr(r.size(), 0.0);
  auto practical_infinity = r.size() - 1;

  auto number_of_iteration = 0u;
  while (number_of_iteration++ < max_iter) {
    /// above r_max the wave function is practically zero. we don't calculate in this region
    practical_infinity = get_practical_infinity(r, z, energy);
    /// classical_turning_point is the boundary between the inward and the outward integration
    auto classical_turning_point = get_classical_turning_point(r, z, energy, practical_infinity);

    auto integrator = AdamsIntegrator(r, z, energy, l, practical_infinity, classical_turning_point);
    integrator.integrate(R, dR_dr);

    /// get the number of R = 0 (not counting the origin).
    auto number_of_nodes = get_number_of_nodes(R, practical_infinity);

    if (required_number_of_nodes != number_of_nodes) {
      /// if it is not the required then we need a bigger energy step.
      energy = update_energy.coarse(number_of_nodes, energy);
    } else {
      /// if the node number is OK, then we fine-tune the energy so that the dR_dr
      /// becomes continuous as well
      auto norm = get_norm(r, dx, R, practical_infinity);
      auto new_R = R[classical_turning_point];
      auto new_dR_dr = dR_dr[classical_turning_point];

      bool finished;
      std::tie(energy, finished) =
        update_energy.fine(energy, norm, new_R, new_dR_dr, integrator.get_old_dR_dr_scaled());
      if (finished) {
        /// if energy diff is small
        normalize_solution(R, dR_dr, norm, practical_infinity);
        break;
      }
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


unsigned int RadialSchrodingerEquation::get_number_of_nodes(
  const std::vector<double> &R, unsigned long practical_infinity)
{
  auto number_of_nodes = 0u;
  auto sign = R[2] > 0.0;

  for (auto i = 3u; i < practical_infinity; ++i) {
    auto new_sign = R[i] > 0.0;
    if (sign != new_sign) {
      ++number_of_nodes;
      sign = new_sign;
    }
  }
  return number_of_nodes;
}

double RadialSchrodingerEquation::get_norm(
  const std::vector<double> &r, double dx, const std::vector<double> &R, unsigned long practical_infinity)
{
  std::vector<double> f;
  f.reserve(practical_infinity);
  for (auto i = 0u; i < practical_infinity; ++i) {
    f.push_back(R[i] * R[i] * r[i]);
  }

  return Math::IntegratorEquidistant::simpson_alt(f, dx);
}

void RadialSchrodingerEquation::normalize_solution(
  std::vector<double> &R, std::vector<double> &dR_dr, double norm, unsigned long practical_infinity)
{
  auto factor = 1.0 / sqrt(norm);
  for (auto i = 0u; i < practical_infinity; ++i) {
    R[i] *= factor;
  }
  for (auto i = 0u; i < practical_infinity; ++i) {
    dR_dr[i] *= factor;
  }
}