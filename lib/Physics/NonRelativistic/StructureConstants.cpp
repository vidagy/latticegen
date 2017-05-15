#include <Math/SphericalHarmonics.h>
#include <Math/Factorial.h>
#include <Math/Bessel.h>
#include <Math/ClebschGordan.h>
#include "StructureConstants.h"

using namespace Physics::NonRelativistic;
using namespace std::complex_literals;

std::complex<double>
RealSpaceStructureConstants::calculate(
  unsigned int l, unsigned int m, unsigned int lprime, unsigned int mprime,
  Coordinates3D n, Coordinates3D nprime,
  const std::complex<double> &z
) const
{
  if (m > l)
    THROW_INVALID_ARGUMENT(" m = " + std::to_string(m) + " is grater than l = " + std::to_string(l));
  if (mprime > lprime)
    THROW_INVALID_ARGUMENT(
      " mprime = " + std::to_string(mprime) + " is grater than lprime = " + std::to_string(lprime));

  auto p = sqrt(z);
  Vector3D R = (std::get<0>(nprime) - std::get<0>(n)) * unit_cell.a
               + (std::get<1>(nprime) - std::get<1>(n)) * unit_cell.b
               + (std::get<2>(nprime) - std::get<2>(n)) * unit_cell.c;

  unsigned int mpp = mprime - m;
  unsigned int lpp_min = static_cast<unsigned int>(abs(lprime - l));
  unsigned int lpp_max = lprime + l;

  auto res = 0i;
  for (auto lpp = lpp_min; lpp <= lpp_max; ++lpp) {
    res += Math::ipow(l - lprime - lpp)
           * Math::hankel_1(lpp, p * R.length())
           * Math::Complex::spherical_harmonic(lpp, mpp, R)
           * Math::Gaunt::calculate(l, m, lpp, mpp, lprime, mprime);

  }
  return -4 * pi * p * 1i * res;
}
