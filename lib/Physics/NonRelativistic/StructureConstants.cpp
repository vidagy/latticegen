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
StructureConstants::calculate_real_space(unsigned int l, unsigned int m, unsigned int lprime, unsigned int mprime,
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

  auto res = 0i;
  for (auto lpp = lpp_min; lpp <= lpp_max; ++lpp) {
    res += ipow(l - lprime - lpp)
           * hankel_1(lpp, p * R.length())
           * Complex::spherical_harmonic(lpp, mpp, R)
           * Gaunt::calculate(l, m, lpp, mpp, lprime, mprime);

  }
  return -4 * pi * p * 1i * res;
}

///@brief Zabloudil et al (15.75)
std::complex<double>
StructureConstants::calculate_reciprocal_space(unsigned int l, unsigned int m, unsigned int lprime, unsigned int mprime,
                                               const Vector3D &k, const std::complex<double> &z) const
{
  if (m > l)
    THROW_INVALID_ARGUMENT(" m = " + std::to_string(m) + " is grater than l = " + std::to_string(l));
  if (mprime > lprime)
    THROW_INVALID_ARGUMENT(
      " mprime = " + std::to_string(mprime) + " is grater than lprime = " + std::to_string(lprime));

  int mpp = m + mprime;
  unsigned int lpp_min = static_cast<unsigned int>(abs(lprime - l));
  unsigned int lpp_max = lprime + l;

  auto res = 0i;
  for (auto lpp = lpp_min; lpp <= lpp_max; ++lpp) {
    res += ipow(l - lprime - lpp)
           * Gaunt::calculate(l, m, lprime, mprime, lpp, mpp)
           * (D1(lpp, mpp, k, z) + D2(lpp, mpp, k, z) + D3(lpp, mpp, k, z));
  }
  return 4 * pi * res;
}

///@brief Zabloudil et al (15.77)
std::complex<double>
StructureConstants::D1(unsigned int l, int m, const Vector3D &k, const std::complex<double> &z) const
{
  auto p = std::sqrt(z);
  const auto &ewald_param = config.ewald_param;

  auto sum = 0i;
  // calculating the sum over reciprocal lattice positions
  for (auto K: reciprocal_mesh) {
    auto K_plus_k = K + k;
    auto K_plus_k_squared = K_plus_k * K_plus_k;
    if (nearlyZero(abs(K_plus_k_squared - z)))
      THROW_LOGIC_ERROR(
        "division by zero: K = " + std::to_string(K)
        + " k = " + std::to_string(k)
        + " z = " + std::to_string(z.real()) + " (+) " + std::to_string(z.imag()) + "i");

    sum += Math::pow(K_plus_k.length(), l)
           * exp(-K_plus_k_squared / ewald_param)
           / (K_plus_k_squared - z)
           * std::conj(Complex::spherical_harmonic(l, m, K_plus_k));
  }

  return -4 * pi
         / std::abs(unit_cell.v1 * cross_product(unit_cell.v2, unit_cell.v3))
         * ipow(l)
         * std::pow(p, -l)
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
    const double step_size = sqrt(std::abs(1.0 / exponent)) / 200.0;
    auto res = 0i;
    auto diff = 0i;
    auto xi_0 = sqrt(config.ewald_param) / 2.0;
    auto iter = 0;
    do {
      auto xi = xi_0 + iter * step_size;
      diff = Math::pow(xi, 2 * l) * std::exp(exponent * xi * xi);
      res += diff;
      ++iter;
    } while ((std::abs(diff) > tol) && iter != max_iter);
    if (iter == max_iter)
      THROW_LOGIC_ERROR("reached max iteration, last diff = "
                        + std::to_string(diff.real()) + " (+) " + std::to_string(diff.imag())
                        + " xi = " + std::to_string(xi_0 + (iter - 1) * step_size));
    //std::cout << "iter = " << iter << std::endl;
    return res;
  }
}

///@brief Zabloudil et al (15.78)
std::complex<double>
StructureConstants::D2(unsigned int l, int m, const Vector3D &k, const std::complex<double> &z) const
{
  auto p = std::sqrt(z);

  auto sum = 0i;
  // calculating the sum over lattice positions
  for (auto R: direct_mesh) {
    auto R_R = R * R;
    auto length = sqrt(R_R);
    if (!nearlyZero(length)) {
      auto diff = Math::pow(length, l)
                  * std::exp(1i * (R * k))
                  * std::conj(Complex::spherical_harmonic(l, m, R))
                  * integral(l, config, R_R, z);
      std::cout << " diff " << diff << std::endl;
      sum += diff;
    }
  }

  return (((l + 1) & 1) == 0 ? 1.0 : -1.0) * (1 << (l + 1)) / sqrt(pi) * std::pow(p, -l) * sum;
}

///@brief Zabloudil et al (15.79)
std::complex<double>
StructureConstants::D3(unsigned int l, int m, const Vector3D &k, const std::complex<double> &z) const
{
  if (l != 0)
    return 0i;
  if (m != 0)
    return 0i;

  auto &ewald_param = config.ewald_param;
  auto &tol = config.integral_tolerance;
  auto res = 0i;
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
  std::vector<Point3D> get_direct_mesh(const Cell3D &cell)
  {
    auto cutoff = 5 * std::max(
      std::max(
        cell.v1.length(),
        cell.v2.length()),
      cell.v3.length()
    );

    return LatticeMesh(cell).generate(CutoffSphere(cutoff));
  }

  std::vector<Point3D> get_reciprocal_mesh(const UnitCell3D &unit_cell)
  {
    auto reciprocal_unit_cell = ReciprocalUnitCell3D(unit_cell);
    return get_direct_mesh(reciprocal_unit_cell);
  }
}

StructureConstants::StructureConstants(const UnitCell3D &unit_cell_, const StructureConstantsConfig config_)
  : unit_cell(unit_cell_), direct_mesh(get_direct_mesh(unit_cell)), reciprocal_mesh(get_reciprocal_mesh(unit_cell)),
    config(config_) {}
