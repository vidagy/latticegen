#include "Neighbors.h"

#include <math.h>
#include <stdexcept>
#include <string>

#include <Core/ComparisonHelpers.h>
using namespace Core;

namespace Geometry
{
  CutoffCube::CutoffCube(const double a_)
    : a(a_) 
  {
    if (! strictlyPositive(a))
      throw std::invalid_argument("In CutoffCube::ctor: a must be non-negative but a = " + std::to_string(a));
  }

  bool CutoffCube::is_included(const Point3D& point) const
  {
    return lessEqualsWithTolerance(fabs(point.x), a) && 
      lessEqualsWithTolerance(fabs(point.y), a) && 
      lessEqualsWithTolerance(fabs(point.z), a);
  }

  CutoffSphere::CutoffSphere(const double r_)
    : r(r_) 
  {
    if (! strictlyPositive(r))
      throw std::invalid_argument("In CutoffSphere::ctor: a must be non-negative but r = " + std::to_string(r));
  }

  bool CutoffSphere::is_included(const Point3D& point) const
  {
    return lessEqualsWithTolerance(point.getLength(), r);
  }
}
