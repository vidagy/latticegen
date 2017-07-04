#include <cmath>
#include <cstdlib>
#include <algorithm>
#include <Core/Exceptions.h>
#include "ClebschGordan.h"
#include "CommonFunctions.h"
#include "SphericalHarmonics.h"

using namespace Math;

double ClebschGordan::calculate(double j1, double m1, double j2, double m2, double J, double M)
{
  int double_j1 = static_cast<int>(lround(2.0 * j1));
  int double_j2 = static_cast<int>(lround(2.0 * j2));
  int double_m1 = static_cast<int>(lround(2.0 * m1));
  int double_m2 = static_cast<int>(lround(2.0 * m2));
  int double_J = static_cast<int>(lround(2.0 * J));
  int double_M = static_cast<int>(lround(2.0 * M));

  double tol = 10 * std::numeric_limits<double>::epsilon();

  if (fabs(2.0 * j1 - double_j1) > tol)
    THROW_INVALID_ARGUMENT("j1 = " + std::to_string(j1) + " is not half or integer");
  if (fabs(2.0 * j2 - double_j2) > tol)
    THROW_INVALID_ARGUMENT("j2 = " + std::to_string(j2) + " is not half or integer");
  if (fabs(2.0 * m1 - double_m1) > tol)
    THROW_INVALID_ARGUMENT("m1 = " + std::to_string(m1) + " is not half or integer");
  if (fabs(2.0 * m2 - double_m2) > tol)
    THROW_INVALID_ARGUMENT("m2 = " + std::to_string(m2) + " is not half or integer");
  if (fabs(2.0 * J - double_J) > tol)
    THROW_INVALID_ARGUMENT("J = " + std::to_string(J) + " is not half or integer");
  if (fabs(2.0 * M - double_M) > tol)
    THROW_INVALID_ARGUMENT("M = " + std::to_string(M) + " is not half or integer");

  if ((double_j1 < 0) || (double_j2 < 0) || (double_J < 0))
    return 0.0;
  if ((abs(double_m1) > double_j1) || (abs(double_m2) > double_j2) || (abs(double_M) > double_J))
    return 0.0;
  if (M != m1 + m2)
    return 0.0;
  if ((double_J < abs(double_j1 - double_j2)) || (double_J > double_j1 + double_j2))
    return 0.0;
  if (((double_j1 + double_j2) & 1) != (double_J & 1))
    return 0.0;

  if (double_M < 0)
    return ((((double_J - double_j1 - double_j2) / 2) & 1) ? -1.0 : 1.0) * calculate(j1, -m1, j2, -m2, J, -M);
  if (double_j1 < double_j2)
    return ((((double_J - double_j1 - double_j2) / 2) & 1) ? -1.0 : 1.0) * calculate(j2, m2, j1, m1, J, M);

  double sum = 0.0;
  int k_min = 0;
  k_min = std::max(k_min, (-double_J + double_j2 - double_m1) / 2);
  k_min = std::max(k_min, (-double_J + double_j1 + double_m2) / 2);
  int k_max = (double_j1 + double_j2 - double_J) / 2;
  k_max = std::min(k_max, (double_j1 - double_m1) / 2);
  k_max = std::min(k_max, (double_j2 + double_m2) / 2);

  for (auto k = k_min; k <= k_max; ++k) {
    sum += ((k & 1) ? -1.0 : 1.0)
           / factorial(k)
           / factorial((double_j1 + double_j2 - double_J) / 2 - k)
           / factorial((double_j1 - double_m1) / 2 - k)
           / factorial((double_j2 + double_m2) / 2 - k)
           / factorial((double_J - double_j2 + double_m1) / 2 + k)
           / factorial((double_J - double_j1 - double_m2) / 2 + k);
  }

  return
    sqrt(
      (double_J + 1)
      * factorial((double_J + double_j1 - double_j2) / 2)
      * factorial((double_J - double_j1 + double_j2) / 2)
      * factorial((-double_J + double_j1 + double_j2) / 2)
      / factorial((double_J + double_j1 + double_j2) / 2 + 1)
    ) * sqrt(
      factorial((double_J + double_M) / 2)
      * factorial((double_J - double_M) / 2)
      * factorial((double_j1 + double_m1) / 2)
      * factorial((double_j1 - double_m1) / 2)
      * factorial((double_j2 + double_m2) / 2)
      * factorial((double_j2 - double_m2) / 2)
    ) * sum;
};

double Gaunt::calculate(double l1, double m1, double l2, double m2, double L, double M)
{
  if ((l1 < 0) || (l2 < 0) || (L < 0))
    return 0.0;

  return
    sqrt((2 * l1 + 1) * (2 * l2 + 1) / 4.0 / pi / (2 * L + 1))
    * ClebschGordan::calculate(l1, 0, l2, 0, L, 0)
    * ClebschGordan::calculate(l1, m1, l2, m2, L, M);
}
