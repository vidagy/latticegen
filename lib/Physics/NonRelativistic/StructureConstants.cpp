#include <Math/SphericalHarmonics.h>
#include <Math/Factorial.h>
#include <Math/Bessel.h>
#include <Math/ClebschGordan.h>
#include <Geometry/Mesh.h>
#include "StructureConstants.h"

using namespace Physics::NonRelativistic;
using namespace std::complex_literals;
using namespace Math;

///@brief Zabloudil et al (15.6)
std::complex<double>
RealStructureConstants::calculate(unsigned int l, unsigned int m, unsigned int lprime, unsigned int mprime,
                                  Coordinates3D n, Coordinates3D nprime, const std::complex<double> &z) const
{
  if (m > l)
    THROW_INVALID_ARGUMENT(" m = " + std::to_string(m) + " is grater than l = " + std::to_string(l));
  if (mprime > lprime)
    THROW_INVALID_ARGUMENT(
      " mprime = " + std::to_string(mprime) + " is grater than lprime = " + std::to_string(lprime));
  if (n == nprime)
    THROW_INVALID_ARGUMENT(
      " n and nprime is the same = " + std::to_string(n));

  auto p = sqrt(z);
  Vector3D R = (nprime.a - n.a) * unit_cell.v1
               + (nprime.b - n.b) * unit_cell.v2
               + (nprime.c - n.c) * unit_cell.v3;

  int mpp = mprime - m;
  unsigned int lpp_min = static_cast<unsigned int>(abs(lprime - l));
  unsigned int lpp_max = lprime + l;

  auto res = 0.0i;
  for (auto lpp = lpp_min; lpp <= lpp_max; ++lpp) {
    res += ipow(l - lprime - lpp)
           * hankel_1(lpp, p * R.length())
           * Complex::spherical_harmonic(lpp, mpp, R)
           * Gaunt::calculate(l, m, lpp, mpp, lprime, mprime);

  }
  return -4 * pi * p * 1.0i * res;
}

///@brief Zabloudil et al (15.75)
std::complex<double>
ReciprocalStructureConstantsCalculator::calculate(unsigned int l, int m, unsigned int lprime, int mprime,
                                                  const Vector3D &k) const
{
  if (static_cast<unsigned int>(abs(m)) > l)
    THROW_INVALID_ARGUMENT(" m = " + std::to_string(m) + " is grater than l = " + std::to_string(l));
  if (static_cast<unsigned int>(abs(mprime)) > lprime)
    THROW_INVALID_ARGUMENT(
      " mprime = " + std::to_string(mprime) + " is grater than lprime = " + std::to_string(lprime));
  if (l > integral_cache.l_max)
    THROW_INVALID_ARGUMENT(
      " l = " + std::to_string(l) + " is grater than l_max = " + std::to_string(integral_cache.l_max)
      + ". Please create a calculator with greater l_max");
  if (lprime > integral_cache.l_max)
    THROW_INVALID_ARGUMENT(
      " lprime = " + std::to_string(lprime) + " is grater than l_max = " + std::to_string(integral_cache.l_max)
      + ". Please create a calculator with greater l_max");

  auto mpp = m + mprime;
  auto lpp_min = static_cast<unsigned int>(abs(lprime - l));
  auto lpp_max = lprime + l;

  auto res = 0.0i;
  for (auto lpp = lpp_min; lpp <= lpp_max; ++lpp) {
    res += ipow(l - lprime - lpp)
           * Gaunt::calculate(l, m, lprime, mprime, lpp, mpp)
           * (D1(lpp, mpp, k) + D2(lpp, mpp, k) + D3(lpp, mpp));
  }
  return 4.0 * pi * res;
}

///@brief Zabloudil et al (15.77)
std::complex<double>
ReciprocalStructureConstantsCalculator::D1(unsigned int l, int m, const Vector3D &k) const
{
  if (nearlyZero(k.length()))
    THROW_INVALID_ARGUMENT("k cannot be zero, k = " + std::to_string(k));

  auto p = std::sqrt(z);
  const auto &ewald_param = config->ewald_param;

  auto sum = 0.0i;
  // calculating the sum over reciprocal lattice positions
  for (auto shell: *reciprocal_shells) {
    auto shell_diff = 0.0i;
    for (auto K: shell.points) {
      auto K_plus_k = K + k;
      auto K_plus_k_squared = K_plus_k * K_plus_k;
      auto K_plus_k_squared_minus_z = K_plus_k_squared - z;
      if (nearlyZero(abs(K_plus_k_squared_minus_z)))
        THROW_LOGIC_ERROR(
          "division by zero: K = " + std::to_string(K)
          + " k = " + std::to_string(k)
          + " z = " + std::to_string(z.real()) + " (+) " + std::to_string(z.imag()) + "i");

      shell_diff += Math::pow(K_plus_k.length(), l)
                    * exp(-K_plus_k_squared / ewald_param)
                    / (K_plus_k_squared_minus_z)
                    * std::conj(Complex::spherical_harmonic(l, m, K_plus_k));
    }
    sum += shell_diff;
    if (abs(shell_diff / sum) < config->integral_tolerance) {
      break;
    }
  }

  return -4.0 * pi
         / std::abs(unit_cell->v1 * cross_product(unit_cell->v2, unit_cell->v3))
         * ipow(l)
         / Math::pow(p, l)
         * std::exp(z / ewald_param)
         * sum;
}

namespace
{
  std::complex<double> integral(int l, const StructureConstantsConfig &config, double R_squared, std::complex<double> z)
  {
    const double tol = config.integral_tolerance;
    const int max_iter = config.max_step_count;
    const auto exponent = z / 4.0 - R_squared;
    const double step_size = sqrt(std::abs(1.0 / exponent)) / config.steps_per_unit;
    auto res = 0.0i;
    auto diff = 0.0i;
    auto xi_0 = sqrt(config.ewald_param) / 2.0;
    auto iter = 0;
    do {
      auto xi = xi_0 + iter * step_size;
      diff = Math::pow(xi, 2 * l) * std::exp(exponent * xi * xi) * step_size;
      res += diff;
      ++iter;
    } while ((std::abs(diff) > tol) && iter != max_iter);
    if (iter == max_iter)
      THROW_LOGIC_ERROR("reached max iteration, last diff = "
                        + std::to_string(diff.real()) + " (+) " + std::to_string(diff.imag())
                        + " xi = " + std::to_string(xi_0 + (iter - 1) * step_size));
    // Logger::log(Logger::Debug, "D2 integral iter: " + std::to_string(iter));
    return res;
  }
}

///@brief Zabloudil et al (15.78)
std::complex<double>
ReciprocalStructureConstantsCalculator::D2(unsigned int l, int m, const Vector3D &k) const
{
  auto p = std::sqrt(z);

  auto sum = 0.0i;
  // calculating the sum over lattice positions
  auto i = 0u;
  for (auto shell: *direct_shells) {
    auto length = shell.r();
    if (!nearlyZero(length)) {
      auto integral = integral_cache.get(l, i++);
      auto shell_diff = 0.0i;
      for (auto R: shell.points) {
        shell_diff += Math::pow(length, l)
                      * std::exp(1.0i * (R * k))
                      * std::conj(Complex::spherical_harmonic(l, m, R))
                      * integral;
      }
      sum += shell_diff;
      if (abs(shell_diff / sum) < config->integral_tolerance) {
        break;
      }
    }
  }

  return (((l + 1) & 1) == 0 ? 1.0 : -1.0) * (1 << (l + 1)) / sqrt(pi) / Math::pow(p, l) * sum;
}

///@brief Zabloudil et al (15.79)
std::complex<double>
ReciprocalStructureConstantsCalculator::D3(unsigned int l, int m) const
{
  if ((l != 0) || (m != 0))
    return 0.0i;

  auto ewald_param = config->ewald_param;
  auto tol = config->integral_tolerance;
  auto res = 0.0i;
  auto diff = 1.0;
  auto n = 0;
  if (std::abs(z) < tol)
    res = -1;
  else {
    while (diff > tol) {
      auto d = pow(z / ewald_param, n) / (2.0 * n - 1.0) / factorial(n);
      res += d;
      diff = std::abs(d) / std::abs(res);
      ++n;
    }
  }
  return -sqrt(ewald_param) / 2.0 / pi * res;
}

namespace
{
  std::shared_ptr<const std::vector<Shell>> get_direct_shells(const Cell3D &cell, double cutoff_scale)
  {
    auto scaled_cutoff = cutoff_scale * std::max(
      std::max(
        cell.v1.length(),
        cell.v2.length()),
      cell.v3.length()
    );

    auto mesh = LatticeMesh(cell).generate(CutoffSphere(scaled_cutoff));
    auto transformations = SymmetryTransformationFactory::generate(cell);

    return std::make_shared<const std::vector<Shell>>(Shell::get_shells(transformations, mesh));
  }

  std::shared_ptr<const std::vector<Shell>> get_reciprocal_shells(const UnitCell3D &unit_cell, double cutoff_scale)
  {
    auto reciprocal_unit_cell = ReciprocalUnitCell3D(unit_cell);
    return get_direct_shells(reciprocal_unit_cell, cutoff_scale);
  }
}

ReciprocalStructureConstants::ReciprocalStructureConstants(
  const UnitCell3D &unit_cell_, const StructureConstantsConfig config_
) : unit_cell(std::make_shared<const UnitCell3D>(unit_cell_)),
    direct_shells(get_direct_shells(unit_cell_, config_.lattice_cutoff_scale)),
    reciprocal_shells(get_reciprocal_shells(unit_cell_, config_.lattice_cutoff_scale)),
    config(std::make_shared<const StructureConstantsConfig>(config_)) {}

ReciprocalStructureConstantsCalculator
ReciprocalStructureConstants::get_calculator(unsigned int l_max, const std::complex<double> &z) const
{
  auto integral_cache = IntegralCache(direct_shells, config, 2 * l_max, z);
  return ReciprocalStructureConstantsCalculator(unit_cell, direct_shells, reciprocal_shells, config, z, integral_cache);
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

namespace
{
  std::vector<std::complex<double>>
  populate_integral_cache(const std::shared_ptr<const std::vector<Shell>> &direct_shells,
                          const std::shared_ptr<const StructureConstantsConfig> &config,
                          const int l_max,
                          const std::complex<double> z)
  {
    auto cache = std::vector<std::complex<double>>();
    cache.reserve(l_max * direct_shells->size());

    for (auto l = 0; l <= l_max; ++l) {
      for (auto shell: *direct_shells) {
        auto r = shell.r();
        // the integral is not defined for r == 0, so we leave that "shell" out.
        // this means that the indexing of the shells are shifted by -1
        if (!nearlyZero(r)) {
          cache.push_back(integral(l, *config, r * r, z));
        }
      }
    }
    return cache;
  }
}

IntegralCache::IntegralCache(const std::shared_ptr<const std::vector<Shell>> &direct_shells,
                             const std::shared_ptr<const StructureConstantsConfig> &config,
                             unsigned int l_max_, std::complex<double> z
) : n_max((const unsigned int) direct_shells->size()), l_max(l_max_),
    cache(populate_integral_cache(direct_shells, config, l_max, z))
{
}

std::complex<double> IntegralCache::get(unsigned int l, unsigned int n) const
{
  if (l > l_max)
    THROW_INVALID_ARGUMENT("l = " + std::to_string(l) + "while l_max = " + std::to_string(l_max));
  if (n > n_max)
    THROW_INVALID_ARGUMENT("n = " + std::to_string(n) + "while n_max = " + std::to_string(n_max));

  return cache[l * n_max + n];
}
