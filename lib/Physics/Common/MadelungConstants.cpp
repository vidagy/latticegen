#include <Math/SphericalHarmonics.h>
#include <Math/ClebschGordan.h>
#include "MadelungConstants.h"

#include <boost/math/special_functions/gamma.hpp>

using namespace Math;
using namespace Physics::Common;
using namespace std::complex_literals;

std::complex<double> RealMadelungConstants::calculate(
  unsigned int l, int m, unsigned int lprime, int mprime, Coordinates3D n, Coordinates3D nprime
) const
{
  if (n == nprime)
    THROW_INVALID_ARGUMENT(
      "n and nprime must be different n = " + std::to_string(n) + " nprime = " + std::to_string(nprime));

  if (l == 0)
    THROW_INVALID_ARGUMENT("l cannot be zero, most probably expression is incorrect");
  if (lprime == 0)
    THROW_INVALID_ARGUMENT("lprime cannot be zero, most probably expression is incorrect");

  auto R_minus_Rprime = n * unit_cell - nprime * unit_cell;

  // TODO this expression might be wrong! (A.4) has different signs then the expression in (19.28), needs to be checked
  auto res =
    Math::sign(l)
    * 4.0 * pi * double_factorial(((int) (2 * (l + lprime))) - 1)
    / double_factorial(((int) (2 * l)) - 1)
    / double_factorial(((int) (2 * lprime)) - 1)
    * ClebschGordan::calculate(l, m, l + lprime, mprime - m, lprime, mprime)
    * std::conj(Complex::spherical_harmonic(l + lprime, mprime - m, R_minus_Rprime))
    / Math::pow(R_minus_Rprime.length(), l + lprime + 1);
  return res;
}

std::complex<double> RealMadelungConstants::calculateReduced(
  unsigned int l, int m, Coordinates3D n, Coordinates3D nprime
) const
{
  if (n == nprime)
    THROW_INVALID_ARGUMENT(
      "n and nprime must be different n = " + std::to_string(n) + " nprime = " + std::to_string(nprime));

  auto R_minus_Rprime = n * unit_cell - nprime * unit_cell;

  auto res =
    Math::sign(l) * 4.0 * pi / (2.0 * l + 1)
    * std::conj(Complex::spherical_harmonic(l, m, R_minus_Rprime))
    / Math::pow(R_minus_Rprime.length(), l + 1);
  return res;
}

namespace
{
  std::vector<Shell> get_direct_shells(const Cell3D &cell, double cutoff_scale)
  {
    auto scaled_cutoff = cutoff_scale * std::max(
      std::max(
        cell.v1.length(),
        cell.v2.length()),
      cell.v3.length()
    );

    auto mesh = LatticeMesh(cell).generate(CutoffSphere(scaled_cutoff));
    auto transformations = SymmetryTransformationFactory::generate(cell);

    return Shell::get_shells(transformations, mesh);
  }

  std::vector<Shell> get_reciprocal_shells(const UnitCell3D &unit_cell, double cutoff_scale)
  {
    auto reciprocal_unit_cell = ReciprocalUnitCell3D(unit_cell);
    return get_direct_shells(reciprocal_unit_cell, cutoff_scale);
  }

  class MadelungConstantCalculator
  {
  public:
    MadelungConstantCalculator(const UnitCell3D &unit_cell_, const MadelungConstantsConfig &config_)
      : unit_cell(unit_cell_), config(config_),
        direct_shells(get_direct_shells(unit_cell, config.lattice_cutoff_scale)),
        reciprocal_shells(get_reciprocal_shells(unit_cell, config.lattice_cutoff_scale)) {}

    lm_vector<std::complex<double>> calculate(unsigned int l_max) const
    {
      auto res = lm_vector<std::complex<double>>(l_max);
      for (auto l = 0u; l <= l_max; ++l) {
        for (auto m = -((int) l); m <= (int) l; ++m) {
          res.at(l, m) = G(l, m);
        }
      }

      return res;
    }

  private:
    ///@brief Zabloudil et al (19.89)
    std::complex<double> G(unsigned int l, int m) const
    {
      return D1(l, m) + D2(l, m) + D3(l, m);
    }

    ///@brief Zabloudil et al (19.82)
    std::complex<double> D1(unsigned int l, int m) const
    {
      auto sum = 0i;
      for (const auto &shell: reciprocal_shells) {
        auto length = shell.r();
        if (!nearlyZero(length)) {
          for (const auto &G: shell.points) {
            sum +=
              std::conj(Complex::spherical_harmonic(l, m, G))
              * Math::pow(length, l) / length / length
              * exp(-length * length * config.ewald_param * config.ewald_param);
          }
        }
      }

      auto V = unit_cell.volume();
      auto last = 0i;
      if (l == 0u && m == 0)
        last = 4.0 * pi * sqrt(4.0 * pi) * config.ewald_param * config.ewald_param / V;
      return sum * 16.0 * pi * pi * sqrt(pi) * Math::ipow(l) / V / static_cast<double>(1ull << (l + 1)) /
             Math::gamma_plus_half(l + 1) +
             last;
    }

    ///@brief Zabloudil et al (19.87)
    std::complex<double> D2(unsigned int l, int m) const
    {
      auto sum = 0i;
      for (const auto &shell: direct_shells) {
        auto length = shell.r();
        if (!nearlyZero(length)) {
          for (const auto &t: shell.points) {
            sum +=
              std::conj(Complex::spherical_harmonic(l, m, -1.0 * t / length))
              / Math::pow(length, l + 1)
              * boost::math::tgamma_lower(l + 0.5, length * length / 4.0 / (config.ewald_param * config.ewald_param));
          }
        }
      }
      return sum * 2.0 * pi * Math::sign(l) / Math::gamma_plus_half(l + 1);
    }

    ///@brief Zabloudil et al (19.88)
    std::complex<double> D3(unsigned int l, int m) const
    {
      if (l == 0u && m == 0)
        return -2.0 / config.ewald_param;
      else
        return 0.0;
    }

    const UnitCell3D unit_cell;
    const MadelungConstantsConfig config;
    const std::vector<Shell> direct_shells;
    const std::vector<Shell> reciprocal_shells;
  };
}

ReducedMadelungConstants::ReducedMadelungConstants(
  const UnitCell3D &unit_cell_, const unsigned int l_max, const MadelungConstantsConfig &config_)
  : unit_cell(unit_cell_), madelung_constants(MadelungConstantCalculator(unit_cell_, config_).calculate(l_max)),
    config(config_) {}
