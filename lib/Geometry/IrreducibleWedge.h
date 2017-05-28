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
      const std::vector<Point3D>& points,
      const SymmetryTransformationFactory::Transformations& transformations);

    static std::vector<Point3D> get_irreducible_wedge(const Cell3D &cell, size_t sample);

    template<typename T>
    static std::vector<std::pair<Point3D, T>> replicate(
      const std::vector<std::pair<Point3D, T>> &irreducible_data,
      const SymmetryTransformationFactory::Transformations &transformations
    )
    {
      auto res = std::vector<std::pair<Point3D, T>>();
      for (auto point: irreducible_data) {
        for (auto transformation: transformations) {
          auto transformed = transformation * point.first;
          auto already_there = false;
          for (auto replicated: res) {
            if (replicated.first == transformed) {
              //std::cout << "===> replicated = " << std::to_string(replicated) << " point = " << std::to_string(replicated) << << std::endl;
              if (!(replicated.second == point.second)) {
                THROW_LOGIC_ERROR("in replicate: inconsistent data at " + std::to_string(point.first) +
                                  " and replicated data: " + std::to_string(replicated.first));
              }
              already_there = true;
            }
          }
          if (!already_there) {
            res.push_back(std::make_pair(transformed, point.second));
            //std::cout << "=> new point = " << std::to_string(transformed) << std::endl;
          }
        }
      }
      return res;
    }
  };
}

#endif //LATTICEGEN_IRREDUCIBLEWEDGE_H
