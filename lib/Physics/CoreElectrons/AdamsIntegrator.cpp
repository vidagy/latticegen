#include <algorithm>
#include <Math/Derivator.h>
#include <Math/LapackWrapper.h>
#include <Math/Factorial.h>
#include "AdamsIntegrator.h"

using namespace Physics::CoreElectrons;

void AdamsIntegrator::integrate(std::vector<double> &R, std::vector<double> &dR_dr)
{
  /// integrating the Schrodinger Equation from practical_infinity to classical_turning_point
  /// this function will overwrite the corresponding region of R and dR_dr
  integrate_inward(R, dR_dr);

  auto old_R_at_ctp = R[classical_turning_point];
  auto old_dR_dr_at_ctp = dR_dr[classical_turning_point];

  /// integrating the Schrodinger Equation from 0 to classical_turning_point
  /// this function will overwrite the corresponding region of R and dR_dr
  integrate_outward(R, dR_dr);

  /// make R and dR_dr continuous
  match_solutions(R, dR_dr, old_R_at_ctp);

  /// the mismatch between the two regions will be used in fine energy update
  old_dR_dr_scaled = R[classical_turning_point] / old_R_at_ctp * old_dR_dr_at_ctp;

}

#define MAX_LENGTH_OF_PARAMS 16

void AdamsIntegrator::adams_moulton_method(
  std::vector<double> &R, std::vector<double> &dR_dr, int from, int to) const
{
  const auto &quadrature = config.adams_moulton_quadrature_order;
  // TODO revisit assertions in AdamsIntegrator (for eg. this one makes no sense)
  if (abs(from - to) <= quadrature)
    THROW_LOGIC_ERROR("in AdamsIntegrator::adams_moulton_method from and to too close: from = " +
                      std::to_string(from) + " to = " + std::to_string(to) + " quadrature = " +
                      std::to_string(quadrature));
  if (quadrature >= MAX_LENGTH_OF_PARAMS)
    THROW_LOGIC_ERROR("quadrature must be smaller than " + std::to_string(MAX_LENGTH_OF_PARAMS));

  const auto &r_points = r->points;
  const auto &dr_di = r->d_points;

  auto diff = from > to ? -1 : 1;
  const auto ang = 0.5 * l * (l + 1);

  /// set up params to contain coefficients multiplied by dx
  auto params = get_adams_parameters(quadrature);

  for (auto i = 0; i <= quadrature; ++i) {
    params[i] *= r->dx;
  }

  /// set up f with initial values
  /// we read quadrature number of values from R and dR_dr ahead of from
  double f_R[MAX_LENGTH_OF_PARAMS];
  double f_dR_dr[MAX_LENGTH_OF_PARAMS];
  for (auto k = 0, i = from - diff * quadrature; k < quadrature; ++k, i += diff) {
    f_R[k] = diff * dr_di[i] * dR_dr[i];
    f_dR_dr[k] = -2.0 * diff * (energy + (z[i] - ang / r_points[i]) / r_points[i]) * dr_di[i] * R[i];
  }

  /// we are using the Adams method on the inclusive range between from and to
  for (auto i = from, k = 0; i != to + diff; i += diff, ++k) {
    auto b = diff * dr_di[i];
    auto c = -2.0 * diff * (energy + (z[i] - ang / r_points[i]) / r_points[i]) * dr_di[i];
    auto param_b = params[quadrature] * b;
    auto param_c = params[quadrature] * c;
    auto determinant = 1.0 - param_b * param_c;
    auto new_R = R[i - diff];
    auto new_dR_dr = dR_dr[i - diff];
    for (auto ii = 0; ii < quadrature; ++ii) {
      new_R += params[ii] * f_R[(ii + k) % quadrature];
      new_dR_dr += params[ii] * f_dR_dr[(ii + k) % quadrature];
    }
    R[i] = (new_R + param_b * new_dR_dr) / determinant;
    dR_dr[i] = (param_c * new_R + new_dR_dr) / determinant;

    /// update f with new values
    f_R[k % quadrature] = b * dR_dr[i];
    f_dR_dr[k % quadrature] = c * R[i];
  }
}

#define MAX_EXPANSION_ORDER 32

void AdamsIntegrator::start_inward(std::vector<double> &R, std::vector<double> &dR_dr) const
{
  const auto &inward_expansion = config.inward_asymptotic_expansion_order;
  if (inward_expansion >= MAX_EXPANSION_ORDER)
    THROW_INVALID_ARGUMENT("inward expansion order " + std::to_string(inward_expansion)
                           + " is not less than" + std::to_string(MAX_EXPANSION_ORDER));

  const auto &cutoff = config.inward_asymptotic_expansion_cutoff;
  const auto &r_points = r->points;

  double alam = sqrt(-2.0 * energy);
  double sigma = z[practical_infinity] / alam;
  double ang = l * (l + 1);

  /// set up coefficients
  double coeff_R[MAX_EXPANSION_ORDER];
  double coeff_dR_dr[MAX_EXPANSION_ORDER];
  coeff_R[0] = 1.0;
  coeff_dR_dr[0] = -alam;
  for (auto i = 1; i < inward_expansion; ++i) {
    coeff_R[i] = (ang - (sigma - i + 1) * (sigma - i)) * coeff_R[i - 1] / (2.0 * i * alam);
    coeff_dR_dr[i] = ((sigma - i + 1) * (sigma + i) - ang) * coeff_R[i - 1] / (2.0 * i);
  }

  /// fill in the first few elements using the expansion formula
  for (auto i = practical_infinity; i > practical_infinity - config.adams_moulton_quadrature_order; --i) {
    auto rfac = pow(r_points[i], sigma) * exp(-alam * r_points[i]);
    auto sum_R = coeff_R[0];
    auto sum_dR_dr = coeff_dR_dr[0];
    auto r_inverses = 1.0;
    for (auto k = 1; k < inward_expansion; ++k) {
      r_inverses /= r_points[i];
      auto R_part = coeff_R[k] * r_inverses;
      auto dR_dr_part = coeff_dR_dr[k] * r_inverses;
      sum_R += R_part;
      sum_dR_dr += dR_dr_part;
      if (std::max(fabs(R_part), fabs(dR_dr_part)) < cutoff)
        break;
    }

    R[i] = rfac * sum_R;
    dR_dr[i] = rfac * sum_dR_dr;
  }
}

void AdamsIntegrator::integrate_inward(std::vector<double> &R, std::vector<double> &dR_dr) const
{
  /// Fill in the first quadrature number of points in R and dR_dr using
  /// the asymptotic expansion of the solution of the Schrodinger Equation
  /// for Coulomb potential. This is valid since far away the potential is Coulomb like.
  start_inward(R, dR_dr);

  /// integrate the rest using the Adams--Moulton formula
  adams_moulton_method(R, dR_dr, practical_infinity - config.adams_moulton_quadrature_order, classical_turning_point);
}

void AdamsIntegrator::start_outward(std::vector<double> &R, std::vector<double> &dR_dr) const
{
  using namespace Math;

  // TODO this is 2 orders less accurate in case of l=0 and n!=1 for dR_dr then for the other cases.
  // Most probably the reason is related to the fact that in those cases dR_dr is finite while R is zero.
  // There must some loss or precision but I'm not sure where is it coming from.
  // A possible solution would be to refine the mesh for the start of the outward integration.
  // For that we would need interpolation for Z as well.

  const auto &outward_scheme = config.outward_quadrature_order;
  const auto &r_points = r->points;
  const auto &dr_di = r->d_points;
  const auto &dx = r->dx;

  double u0 = 1.0;
  double v0 = -z[0] / (l + 1.0);

  auto lagrangian_quadrature = Derivator::lagrange_quadrature(outward_scheme + 1);

  for (auto scheme_repetition = 0; scheme_repetition < config.outward_scheme_repetition; ++scheme_repetition) {
    auto starting_index = scheme_repetition * outward_scheme;

    /// Preparation set up m matrix
    std::vector<double> b(outward_scheme, 0.0);
    std::vector<double> c(outward_scheme, 0.0);
    std::vector<double> d(outward_scheme, 0.0);
    std::vector<double> m(outward_scheme * outward_scheme, 0.0);
    for (auto i = 0; i < outward_scheme; ++i) {
      b[i] = dx * dr_di[i + starting_index];
      c[i] = -2.0 * dx * (energy + z[i + starting_index] / r_points[i + starting_index]) * dr_di[i + starting_index];
      d[i] = -2.0 * dx * (l + 1.0) * dr_di[i + starting_index] / r_points[i + starting_index];

      for (auto j = 0; j < outward_scheme; ++j) {
        if (i == j)
          m[i * outward_scheme + j] = lagrangian_quadrature[i + 1][j + 1] - d[i];
        else
          m[i * outward_scheme + j] = lagrangian_quadrature[i + 1][j + 1];
      }
    }

    /// invert (lagrangian_quadrature - d * id) matrix
    LapackWrapper::invert_matrix(m);

    /// solve equation for R
    /// set up matrix and rhs
    std::vector<double> fm(outward_scheme * outward_scheme, 0.0);
    std::vector<double> rhs(outward_scheme, 0.0);
    for (auto i = 0; i < outward_scheme; ++i) {
      rhs[i] = -lagrangian_quadrature[i + 1][0] * u0;
      for (auto j = 0; j < outward_scheme; ++j) {
        fm[i * outward_scheme + j] = lagrangian_quadrature[i + 1][j + 1] - b[i] * m[i * outward_scheme + j] * c[j];
        rhs[i] += -b[i] * m[i * outward_scheme + j] * lagrangian_quadrature[j + 1][0] * v0;
      }
    }

    /// invert matrix
    LapackWrapper::invert_matrix(fm);

    /// multiply rhs;
    std::vector<double> R_solution(outward_scheme, 0.0);
    for (auto i = 0; i < outward_scheme; ++i) {
      for (auto j = 0; j < outward_scheme; ++j) {
        R_solution[i] += fm[i * outward_scheme + j] * rhs[j];
      }
    }

    /// solve for dR_dr
    std::vector<double> dR_dr_solution(outward_scheme, 0.0);
    for (auto i = 0; i < outward_scheme; ++i) {
      for (auto j = 0; j < outward_scheme; ++j) {
        dR_dr_solution[i] +=
          m[i * outward_scheme + j] * (c[j] * R_solution[j] - lagrangian_quadrature[j + 1][0] * v0);
      }
    }

    /// finalize solutions
    for (auto i = 0; i < outward_scheme; ++i) {
      R[i + starting_index] = Math::pow(r_points[i + starting_index], l + 1) * R_solution[i];
      dR_dr[i + starting_index] = Math::pow(r_points[i + starting_index], l) *
                                  (r_points[i + starting_index] * dR_dr_solution[i] + (l + 1.0) * R_solution[i]);
    }
    u0 = R_solution.back();
    v0 = dR_dr_solution.back();
  }
}

void AdamsIntegrator::integrate_outward(std::vector<double> &R, std::vector<double> &dR_dr) const
{
  /// use Lagrangian differentiation formula to obtain the first quadrature number of points
  start_outward(R, dR_dr);

  /// integrate the rest using the Adams--Moulton formula
  adams_moulton_method(R, dR_dr, config.adams_moulton_quadrature_order, classical_turning_point);
}

void AdamsIntegrator::match_solutions(std::vector<double> &R, std::vector<double> &dR_dr, double old_R_at_ctp) const
{
  auto R_ratio = R[classical_turning_point] / old_R_at_ctp;

  for (auto i = classical_turning_point + 1; i < practical_infinity; ++i) {
    R[i] *= R_ratio;
    dR_dr[i] *= R_ratio;
  }
}

namespace
{
  double get_adams_parameter(int n, int i)
  {
    // i = 0 ... n
    double prefactor = (((n - i) & 1) ? -1.0 : 1.0) / Math::factorial(i) / Math::factorial(n - i);

    std::vector<double> factors(n, 0.0);
    factors[0] = 1.0;
    for (auto j = 0; j <= n; ++j) {
      if (j != n - i) {
        double current = j - 1;

        for (auto k = n - 1; k > 0; --k) {
          factors[k] = current * factors[k] + factors[k - 1];
        }
        factors[0] *= current;
      }
    }
    double integral = 0.0;
    for (int j = 0; j < n; ++j) {
      integral += factors[j] / (j + 1.0);
    }
    integral += 1.0 / (n + 1.0);

    return prefactor * integral;
  }
}

std::vector<double> AdamsIntegrator::get_adams_parameters(int quadrature)
{
  if (quadrature < 1)
    THROW_INVALID_ARGUMENT("quadrature must be positive");

  std::vector<double> res;
  res.reserve(quadrature + 1);
  for (auto i = 0; i <= quadrature; ++i)
    res.push_back(get_adams_parameter(quadrature, i));
  return res;
}
