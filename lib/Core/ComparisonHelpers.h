#ifndef LATTICEGEN_COMPARISONHELPERS_H_
#define LATTICEGEN_COMPARISONHELPERS_H_

#include <math.h>

namespace Core
{
  static const double default_abs_tolerance = 1e-14;
  static const double default_rel_tolerance = 1e-10;

  inline bool equalsWithTolerance(
    double lhs, double rhs,
    double abs_tolerance = default_abs_tolerance, double rel_tolerance = default_rel_tolerance)
  {
    double diff = fabs(lhs - rhs);
    if (diff > abs_tolerance)
      return false;
    
    double sum = (lhs + rhs) / 2.0;
    return sum < abs_tolerance || diff / sum <= rel_tolerance;
      
  }

  inline bool nearlyZero(double x, double abs_tolerance = default_abs_tolerance)
  {
    return x < abs_tolerance && x > -1.0 * abs_tolerance;
  }

  inline bool strictlyGreater(double x, double y, double abs_tolerance = default_abs_tolerance)
  {
    return x > y + abs_tolerance;
  }

  inline bool strictlyLess(double x, double y, double abs_tolerance = default_abs_tolerance)
  {
    return strictlyGreater(y, x, abs_tolerance);
  }

  inline bool greaterEqualsWithTolerance(double x, double y, double abs_tolerance = default_abs_tolerance)
  {
    return x > y - abs_tolerance;
  }

  inline bool lessEqualsWithTolerance(double x, double y, double abs_tolerance = default_abs_tolerance)
  {
    return greaterEqualsWithTolerance(y, x, abs_tolerance);
  }

  inline bool strictlyPositive(double x, double abs_tolerance = default_abs_tolerance)
  {
    return x > abs_tolerance;
  }

  inline bool strictlyNegative(double x, double abs_tolerance = default_abs_tolerance)
  {
    return x < -abs_tolerance;
  }

  inline bool positiveWithTolerance(double x, double abs_tolerance = default_abs_tolerance)
  {
    return ! strictlyNegative(x, abs_tolerance);
  }

  inline bool negativeWithTolerance(double x, double abs_tolerance = default_abs_tolerance)
  {
    return ! strictlyPositive(x, abs_tolerance);
  }
}
#endif // LATTICEGEN_COMPARISONHELPERS_H_
 