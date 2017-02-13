#include "IrreducibleWedge.h"

using namespace Geometry;

std::vector<Point3D>
IrreducibleWedge::reduce_by_symmetries(const std::vector<Point3D> &points,
                                       const SymmetryTransformationFactory::Transformations &transformations)
{
  std::vector<Point3D> result;
  result.reserve(points.size());
  std::vector<bool> already_covered(points.size(), false);
  for (size_t i = 0; i < points.size(); ++i)
  {
    if (! already_covered[i])
    {
      result.push_back(points[i]);

      for (size_t j = i+1; j < points.size(); ++j)
      {
        if (! already_covered[j])
        {
          for (auto transformation: transformations)
          {
            if (points[j] == transformation*points[i])
            {
              already_covered[j] = true;
              break;
            }
          }
        }
      }
    }
    already_covered[i] = true;
  }
  return result;
}

//std::vector<Point3D>
//IrreducibleWedge::get_irreducible_wedge(const UnitCell3D &unit_cell, const double sample_width)
//{
//  return std::vector<Point3D>{};
//}
