#ifndef LATTICEGEN_SHELL_H
#define LATTICEGEN_SHELL_H

#include <Core/Exceptions.h>
#include <Core/Point3D.h>
#include <utility>
#include <vector>
#include "Mesh.h"
#include "SymmetryTransformationFactory.h"

using namespace Core;

namespace Geometry
{
  class Shell
  {
  public:
    explicit Shell(std::vector<Point3D> points_) : points(std::move(points_))
    {
      if (points.empty())
        THROW_INVALID_ARGUMENT("Cannot create shell with empty input");
      auto r = points[0].norm();
      for (const auto &p: points) {
        if (!equalsWithTolerance(r, p.norm()))
          THROW_INVALID_ARGUMENT(
            "Shell with different distances: r = " + std::to_string(r) +
            " p.length() = " + std::to_string(p.norm())
          );
      }
    }

    const std::vector<Point3D> points;

    double r() const
    {
      return points[0].norm();
    }

    static std::vector<Shell> get_shells(
      SymmetryTransformationFactory::Transformations &transformations, const std::vector<Point3D> &mesh);
  };
}
#endif //LATTICEGEN_SHELL_H
