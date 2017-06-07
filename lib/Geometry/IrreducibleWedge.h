#ifndef LATTICEGEN_IRREDUCIBLEWEDGE_H
#define LATTICEGEN_IRREDUCIBLEWEDGE_H

#include <vector>
#include <Core/Point3D.h>
#include "SymmetryTransformationFactory.h"
#include "UnitCell3D.h"

using namespace Core;

namespace Geometry
{
  class IrreducibleWedge
  {
  public:
    static std::vector<Point3D> reduce_by_symmetries(
      std::vector<Point3D> points,
      const SymmetryTransformationFactory::Transformations &transformations,
      double abs_tolerance = default_abs_tolerance);

    static std::vector<Point3D> get_irreducible_wedge(const Cell3D &cell, size_t sample, bool centered = false);

    static double get_tolerance(const Cell3D &cell);

    template<typename T>
    static std::vector<std::pair<Point3D, T>> replicate(
      const std::vector<std::pair<Point3D, T>> &irreducible_data,
      const SymmetryTransformationFactory::Transformations &transformations,
      double abs_tolerance = default_abs_tolerance
    )
    {
      auto comparator = Point3DComparator(abs_tolerance, 1.0);
      auto res = std::vector<std::pair<Point3D, T>>();
      for (auto point: irreducible_data) {
        auto sub_res = std::vector<std::pair<Point3D, T>>();
        for (auto transformation: transformations) {
          auto transformed = transformation * point.first;
          auto already_there = false;
          for (auto replicated: sub_res) {
            if (comparator.isEqual(replicated.first, transformed)) {
              if (!(replicated.second == point.second)) {
                THROW_LOGIC_ERROR("in replicate: inconsistent data at " + std::to_string(point.first) +
                                  " and replicated data: " + std::to_string(replicated.first));
              }
              already_there = true;
            }
          }
          if (!already_there) {
            sub_res.push_back(std::make_pair(transformed, point.second));
          }
        }
        std::move(sub_res.begin(), sub_res.end(), std::back_inserter(res));
      }
      return res;
    }
  };
}

#endif //LATTICEGEN_IRREDUCIBLEWEDGE_H
