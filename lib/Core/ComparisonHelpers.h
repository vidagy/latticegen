#ifndef LATTICEGEN_COMPARISONHELPERS_H_
#define LATTICEGEN_COMPARISONHELPERS_H_

#include <math.h>

namespace Core
{
  static const double ABS_TOLERANCE = 1e-15;
  static const double REL_TOLERANCE = 1e-10;

  inline bool equalsWithTolerance(
    const double lhs, const double rhs, 
    const double abs_tolerance = ABS_TOLERANCE, const double rel_tolerance = REL_TOLERANCE)
  {
    double diff = fabs(lhs - rhs);
    if (diff > abs_tolerance)
      return false;
    
    double sum = (lhs + rhs) / 2.0;
    if (diff / sum > rel_tolerance)
      return false;

    return true;
  }

  inline bool nearlyZero(const double x, const double abs_tolerance = ABS_TOLERANCE)
  {
    return x < abs_tolerance && x > -1.0 * abs_tolerance;
  }
}
#endif // LATTICEGEN_COMPARISONHELPERS_H_
 
