#include <algorithm>
#include "Shell.h"

using namespace Geometry;

namespace
{
  bool is_same_shell(
    Point3D &lhs, Point3D &rhs, SymmetryTransformationFactory::Transformations &transformations
  )
  {
    for (auto &transformation: transformations) {
      if (lhs.isApprox(transformation * rhs))
        return true;
    }
    return false;
  }

  std::vector<Shell> get_shells(
    const std::vector<Point3D>::iterator &begin,
    const std::vector<Point3D>::iterator &end,
    SymmetryTransformationFactory::Transformations &transformations
  )
  {
    if (begin == end)
      THROW_INVALID_ARGUMENT("empty range");

    auto collector = std::vector<std::vector<Point3D>>{std::vector<Point3D>{*begin}};

    for (auto it = begin + 1; it != end; ++it) {
      auto ok = false;
      for (auto &shell: collector) {
        if (is_same_shell(*it, shell.front(), transformations)) {
          shell.push_back(*it);
          ok = true;
          break;
        }
      }
      if (!ok) {
        collector.push_back({*it});
      }
    }

    auto res = std::vector<Shell>();
    std::transform(collector.begin(), collector.end(), std::back_inserter(res),
                   [](const std::vector<Point3D> &points)
                   {
                     return Shell(points);
                   });
    return res;
  }
}

std::vector<Shell> Shell::get_shells(
  SymmetryTransformationFactory::Transformations &transformations, const std::vector<Point3D> &mesh)
{
  auto res = std::vector<Shell>();

  if (mesh.empty())
    return res;

  auto sorted_mesh = mesh;
  std::sort(sorted_mesh.begin(), sorted_mesh.end(),
            [](const Point3D &lhs, const Point3D &rhs)
            {
              return lhs.norm() < rhs.norm();
            });

  auto pre_length = sorted_mesh[0].norm();
  auto begin = sorted_mesh.begin();
  for (auto it = sorted_mesh.begin(); it != sorted_mesh.end(); ++it) {
    if (!equalsWithTolerance(pre_length, it->norm())) {
      auto current_shells = ::get_shells(begin, it, transformations);
      std::copy(current_shells.begin(), current_shells.end(), std::back_inserter(res));
      pre_length = it->norm();
      begin = it;
    }
  }
  auto current_shells = ::get_shells(begin, sorted_mesh.end(), transformations);
  std::copy(current_shells.begin(), current_shells.end(), std::back_inserter(res));

  return res;
}
