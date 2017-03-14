#include <algorithm>
#include <Math/Derivator.h>
#include <Math/LapackWrapper.h>
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

namespace
{
  // TODO write a generator for these coefficients

  static const double adams_params_1[1] =
    {0.5};
  static const double adams_params_2[2] =
    {-1.0 / 12.0, 8.0 / 12.0};
  static const double adams_params_3[3] =
    {1.0 / 24.0, -5.0 / 24.0, 19.0 / 24.0};
  static const double adams_params_4[4] =
    {-19.0 / 720.0, 106.0 / 720.0, -264.0 / 720.0, 646.0 / 720.0};
  static const double adams_params_5[5] =
    {27.0 / 1440.0, -173.0 / 1440.0, 482.0 / 1440.0, -798.0 / 1440.0, 1427.0 / 1440.0};
  static const double adams_params_6[6] =
    {-863.0 / 60480.0, 6312 / 60480.0, -20211.0 / 60480.0,
     37504.0 / 60480.0, -46461.0 / 60480.0, 65112.0 / 60480.0};
  static const double adams_params_7[7] =
    {1375.0 / 120960.0, -11351.0 / 120960.0, 41499.0 / 120960.0,
     -88547.0 / 120960.0, 123133.0 / 120960.0, -121797.0 / 120960.0, 139849.0 / 120960.0};
  static const double adams_params_8[8] =
    {-33953.0 / 3628800.0, 312874.0 / 3628800.0, -1291214.0 / 3628800.0, 3146338.0 / 3628800.0,
     -5033120.0 / 3628800.0, 5595358.0 / 3628800.0, -4604594.0 / 3628800.0, 4467094.0 / 3628800.0};

  static const double *adams_params[8] =
    {adams_params_1, adams_params_2, adams_params_3, adams_params_4,
     adams_params_5, adams_params_6, adams_params_7, adams_params_8};

  static const double adams_param[8] =
    {0.5, 5.0 / 12.0, 9.0 / 24.0, 251.0 / 720.0,
     475.0 / 1440.0, 19087.0 / 60480.0, 36799.0 / 120960.0, 1070017.0 / 3628800.0};
}
#define MAX_LENGTH_OF_PARAMS 8

void AdamsIntegrator::adams_moulton_method(
  std::vector<double> &R, std::vector<double> &dR_dr, unsigned long from, unsigned long to) const
{
  if (abs(static_cast<int>(from) - static_cast<int>(to)) <= static_cast<int>(quadrature))
    throw std::logic_error("in AdamsIntegrator::adams_moulton_method from and to too close: from = " +
                           std::to_string(from) + " to = " + std::to_string(to) + " quadrature = " +
                           std::to_string(quadrature));

  const auto &r_points = r->points;
  auto diff = from > to ? -1l : 1l;
  const auto ang = 0.5 * l * (l + 1);

  /// set up param and params to contain coefficients multiplied by dx
  const auto &dx = r->dx;
  const auto &param = dx * adams_param[quadrature - 1];

  double params[MAX_LENGTH_OF_PARAMS];
  for (auto i = 0u; i < quadrature; ++i) {
    params[i] = dx * adams_params[quadrature - 1][i];
  }

  /// set up f with initial values
  /// we read quadrature number of values from R and dR_dr ahead of from
  double f_R[MAX_LENGTH_OF_PARAMS];
  double f_dR_dr[MAX_LENGTH_OF_PARAMS];
  for (auto k = 0ul, i = from - diff * quadrature; k < quadrature; ++k, i += diff) {
    f_R[k] = static_cast<double>(diff) * r_points[i] * dR_dr[i];
    f_dR_dr[k] = -2.0 * static_cast<double>(diff) * (energy * r_points[i] + z[i] - ang / r_points[i]) * R[i];
  }

  /// we are using the Adams method on the inclusive range between from and to
  for (auto i = from, k = 0ul; i != to + diff; i += diff, ++k) {
    auto b = static_cast<double>(diff) * r_points[i];
    auto c = -2.0 * static_cast<double>(diff) * (energy * r_points[i] + z[i] - ang / r_points[i]);
    auto param_b = param * b;
    auto param_c = param * c;
    auto determinant = 1.0 - param_b * param_c;
    auto new_R = R[i - diff];
    auto new_dR_dr = dR_dr[i - diff];
    for (auto ii = 0u; ii < quadrature; ++ii) {
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

void AdamsIntegrator::start_inward(std::vector<double> &R, std::vector<double> &dR_dr) const
{
#define INWARD_ASYMPTOTIC_EXPANSION 15
#define TOLERANCE 1e-14

  const auto &r_points = r->points;

  double alam = sqrt(-2.0 * energy);
  double sigma = z[practical_infinity] / alam;
  double ang = l * (l + 1);

  /// set up coefficients
  double coeff_R[INWARD_ASYMPTOTIC_EXPANSION];
  double coeff_dR_dr[INWARD_ASYMPTOTIC_EXPANSION];
  coeff_R[0] = 1.0;
  coeff_dR_dr[0] = -alam;
  for (auto i = 1; i < INWARD_ASYMPTOTIC_EXPANSION; ++i) {
    coeff_R[i] = (ang - (sigma - i + 1) * (sigma - i)) * coeff_R[i - 1] / (2.0 * i * alam);
    coeff_dR_dr[i] = ((sigma - i + 1) * (sigma + i) - ang) * coeff_R[i - 1] / (2.0 * i);
  }

  /// fill in the first few elements using the expansion formula
  for (auto i = practical_infinity; i > practical_infinity - quadrature; --i) {
    auto rfac = pow(r_points[i], sigma) * exp(-alam * r_points[i]);
    auto sum_R = coeff_R[0];
    auto sum_dR_dr = coeff_dR_dr[0];
    auto r_inverses = 1.0;
    for (auto k = 1; k < INWARD_ASYMPTOTIC_EXPANSION; ++k) {
      r_inverses /= r_points[i];
      auto R_part = coeff_R[k] * r_inverses;
      auto dR_dr_part = coeff_dR_dr[k] * r_inverses;
      sum_R += R_part;
      sum_dR_dr += dR_dr_part;
      if (std::max(fabs(R_part), fabs(dR_dr_part)) < TOLERANCE)
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
  adams_moulton_method(R, dR_dr, practical_infinity - quadrature, classical_turning_point);
}

void AdamsIntegrator::start_outward(std::vector<double> &R, std::vector<double> &dR_dr) const
{
  using namespace Math;
#define OUTWARD_SCHEME_REPETITION 1
#define OUTWARD_SCHEME 8

  const auto &r_points = r->points;
  const auto &dx = r->dx;

  double u0 = 1.0;
  double v0 = -z[0] / (static_cast<double>(l) + 1.0);

  auto lagrangian_quadrature = Derivator::lagrange_quadrature(OUTWARD_SCHEME + 1);

  for (auto scheme_repetition = 0; scheme_repetition < OUTWARD_SCHEME_REPETITION; ++scheme_repetition) {
    auto starting_index = scheme_repetition * OUTWARD_SCHEME;

    /// Preparation set up m matrix
    std::vector<double> b(OUTWARD_SCHEME, 0.0);
    std::vector<double> c(OUTWARD_SCHEME, 0.0);
    std::vector<double> d(OUTWARD_SCHEME, 0.0);
    std::vector<double> m(OUTWARD_SCHEME * OUTWARD_SCHEME, 0.0);
    for (auto i = 0; i < OUTWARD_SCHEME; ++i) {
      b[i] = dx * r_points[i + starting_index];
      c[i] = -2.0 * dx * (energy * r_points[i + starting_index] + z[i + starting_index]);
      d[i] = -2.0 * dx * (static_cast<double>(l) + 1.0);

      for (auto j = 0; j < OUTWARD_SCHEME; ++j) {
        if (i == j)
          m[i * OUTWARD_SCHEME + j] = lagrangian_quadrature[i + 1][j + 1] - d[i];
        else
          m[i * OUTWARD_SCHEME + j] = lagrangian_quadrature[i + 1][j + 1];
      }
    }

    /// invert (lagrangian_quadrature - d * id) matrix
    LapackWrapper::invert_matrix(m);

    /// solve equation for R
    /// set up matrix and rhs
    std::vector<double> fm(OUTWARD_SCHEME * OUTWARD_SCHEME, 0.0);
    std::vector<double> rhs(OUTWARD_SCHEME, 0.0);
    for (auto i = 0; i < OUTWARD_SCHEME; ++i) {
      rhs[i] = -lagrangian_quadrature[i + 1][0] * u0;
      for (auto j = 0; j < OUTWARD_SCHEME; ++j) {
        fm[i * OUTWARD_SCHEME + j] = lagrangian_quadrature[i + 1][j + 1] - b[i] * m[i * OUTWARD_SCHEME + j] * c[j];
        rhs[i] += -b[i] * m[i * OUTWARD_SCHEME + j] * lagrangian_quadrature[j + 1][0] * v0;
      }
    }

    /// invert matrix
    LapackWrapper::invert_matrix(fm);

    /// multiply rhs;
    std::vector<double> R_solution(OUTWARD_SCHEME, 0.0);
    for (auto i = 0; i < OUTWARD_SCHEME; ++i) {
      for (auto j = 0; j < OUTWARD_SCHEME; ++j) {
        R_solution[i] += fm[i * OUTWARD_SCHEME + j] * rhs[j];
      }
    }

    /// solve for dR_dr
    std::vector<double> dR_dr_solution(OUTWARD_SCHEME, 0.0);
    for (auto i = 0; i < OUTWARD_SCHEME; ++i) {
      for (auto j = 0; j < OUTWARD_SCHEME; ++j) {
        dR_dr_solution[i] +=
          m[i * OUTWARD_SCHEME + j] * (c[j] * R_solution[j] - lagrangian_quadrature[j + 1][0] * v0);
      }
    }

    /// finalize solutions
    for (auto i = 0; i < OUTWARD_SCHEME; ++i) {
      R[i + starting_index] = pow(r_points[i + starting_index], l + 1) * R_solution[i];
      dR_dr[i + starting_index] = pow(r_points[i + starting_index], l) *
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
  adams_moulton_method(R, dR_dr, quadrature, classical_turning_point);
}

void AdamsIntegrator::match_solutions(std::vector<double> &R, std::vector<double> &dR_dr, double old_R_at_ctp) const
{
  auto R_ratio = R[classical_turning_point] / old_R_at_ctp;

  for (auto i = classical_turning_point + 1; i < practical_infinity; ++i) {
    R[i] *= R_ratio;
    dR_dr[i] *= R_ratio;
  }
}
