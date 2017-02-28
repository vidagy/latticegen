#include "AdamsIntegrator.h"

using namespace Physics::CoreElectrons;

void AdamsIntegrator::integrate(std::vector<double> &R, std::vector<double> &dR_dr)
{
  /// integrating the Schrodinger Equation from r_max to classical_turning_point
  /// this function will overwrite the corresponding region of R and dR_dr
  integrate_inward(R, dR_dr);

  auto old_R_at_ctp = R[classical_turning_point];
  auto old_dR_dr_at_ctp = dR_dr[classical_turning_point];

  /// integrating the Schrodinger Equation from 0 to classical_turning_point
  /// this function will overwrite the corresponding region of R and dR_dr
  integrate_outward(R, dR_dr);

  /// make R and dR_dr continuous
  match_solutions(R, dR_dr, old_R_at_ctp, old_dR_dr_at_ctp);

  old_dR_dr_scaled = R[classical_turning_point] / old_R_at_ctp * old_dR_dr_at_ctp;
}

void AdamsIntegrator::adams_integrator(std::vector<double> &R, std::vector<double> &dR_dr) const
{

}

void AdamsIntegrator::integrate_inward(std::vector<double> &R, std::vector<double> &dR_dr) const
{

}

void AdamsIntegrator::integrate_outward(std::vector<double> &R, std::vector<double> &dR_dr) const
{

}

void AdamsIntegrator::match_solutions(
  std::vector<double> &R, std::vector<double> &dR_dr, double old_R_at_ctp, double old_dR_dr_at_ctp
) const
{
  auto R_ratio = R[classical_turning_point] / old_R_at_ctp;
  auto dR_dr_ratio = R_ratio * old_dR_dr_at_ctp;

  for (auto i = classical_turning_point; i < practical_infinity; ++i) {
    R[i] *= R_ratio;
    dR_dr[i] *= dR_dr_ratio;
  }
}
