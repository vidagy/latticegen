#ifndef LATTICEGEN_TEST_TEST_UTILS_BASE_H
#define LATTICEGEN_TEST_TEST_UTILS_BASE_H

#include <gtest/gtest.h>
#include <gmock/gmock.h>
#include <cmath>

#include <Core/Point3D.h>
#include <Core/Matrix3D.h>

MATCHER_P(NearWithTolerance, tolerance, "") {
  *result_listener << "the difference is " << std::setprecision(18) << std::get<0>(arg) - std::get<1>(arg);
  return std::fabs(std::get<0>(arg) - std::get<1>(arg)) <= tolerance;
}

namespace Testing {
  struct Point3DWrap : public Core::Point3D {
    Point3DWrap(double a1, double a2, double a3) : Core::Point3D(a1, a2, a3) {}

    Point3DWrap(const Core::Point3D &point) : Point3DWrap(point[0], point[1], point[2]) {}

    bool operator==(const Point3DWrap &other) const {
      if (Core::nearlyZero(this->norm()) || Core::nearlyZero(other.norm()))
        return Core::nearlyZero((*this - other).norm());
      return this->isApprox(other);
    }
  };

  inline std::ostream &operator<<(std::ostream &o, const Point3DWrap &point) {
    return o << static_cast<Core::Point3D>(point) << "\n";
  }

  inline std::vector<Point3DWrap> wrap(const std::vector<Core::Point3D> &points) {
    auto res = std::vector<Point3DWrap>();
    res.reserve(points.size());
    for (const auto &p : points) {
      res.emplace_back(Point3DWrap(p));
    }
    return res;
  }
}

#endif //LATTICEGEN_TEST_TEST_UTILS_BASE_H
