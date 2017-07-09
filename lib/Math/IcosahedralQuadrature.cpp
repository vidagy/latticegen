#include <algorithm>
#include "IcosahedralQuadrature.h"

using namespace Math;

namespace
{
  struct Triangle
  {
    Triangle(const Point3D &a_, const Point3D &b_, const Point3D &c_) : a(a_), b(b_), c(c_) {}

    void remap()
    {
      a /= a.length();
      b /= b.length();
      c /= c.length();
    };

    Point3D a;
    Point3D b;
    Point3D c;
  };

  void remap(std::vector<Triangle> &faces)
  {
    std::for_each(
      faces.begin(), faces.end(),
      [](Triangle &face)
      {
        return face.remap();
      }
    );
  }

  std::vector<Triangle> icosahedronFaces()
  {
    static const auto p = (1.0 + sqrt(5.0)) / 2.0;
    // a, b and c corresponds to the three rectangles
    static const auto a1 = Point3D{0, p, 1};
    static const auto a2 = Point3D{0, p, -1};
    static const auto a3 = Point3D{0, -p, 1};
    static const auto a4 = Point3D{0, -p, -1};
    static const auto b1 = Point3D{1, 0, p};
    static const auto b2 = Point3D{-1, 0, p};
    static const auto b3 = Point3D{1, 0, -p};
    static const auto b4 = Point3D{-1, 0, -p};
    static const auto c1 = Point3D{p, 1, 0};
    static const auto c2 = Point3D{p, -1, 0};
    static const auto c3 = Point3D{-p, 1, 0};
    static const auto c4 = Point3D{-p, -1, 0};

    auto faces = std::vector<Triangle>{
      {a1, a2, c1},
      {a1, a2, c3},
      {a3, a4, c2},
      {a3, a4, c4},
      {b1, b2, a1},
      {b1, b2, a3},
      {b3, b4, a2},
      {b3, b4, a4},
      {c1, c2, b1},
      {c1, c2, b3},
      {c3, c4, b2},
      {c3, c4, b4},
      {c2, a4, b3},
      {c1, a2, b3},
      {c4, a4, b4},
      {c3, a2, b4},
      {c2, a3, b1},
      {c1, a1, b1},
      {c4, a3, b2},
      {c3, a1, b2}
    };

    remap(faces);
    return faces;
  }

  std::vector<Triangle> refine(std::vector<Triangle> &faces)
  {
    auto res = std::vector<Triangle>();
    res.reserve(faces.size() * 4);
    std::for_each(
      faces.begin(), faces.end(),
      [&res](const Triangle &face)
      {
        auto ab = (face.a + face.b) / 2.0;
        auto ac = (face.a + face.c) / 2.0;
        auto bc = (face.b + face.c) / 2.0;

        res.push_back(Triangle{face.a, ab, ac});
        res.push_back(Triangle{face.b, ab, bc});
        res.push_back(Triangle{face.c, ac, bc});
        res.push_back(Triangle{ab, ac, bc});
      }
    );
    return res;
  }

  std::vector<std::pair<Point3D, double>> centers_and_weights(const std::vector<Triangle> &faces)
  {
    auto res = std::vector<std::pair<Point3D, double>>();
    res.reserve(faces.size());

    std::transform(
      faces.begin(), faces.end(), std::back_inserter(res),
      [](const Triangle &face)
      {
        auto p = (face.a + face.b + face.c) / 3.0;
        auto w = fabs((face.a - face.b) * (face.a - face.c) / 2.0);
        return std::make_pair(p / p.length(), w);
      }
    );

    return res;
  };

  void scale_weights(std::vector<std::pair<Point3D, double>> &centers_and_weights)
  {
    auto sum_weights = std::accumulate(
      centers_and_weights.begin(), centers_and_weights.end(), 0.0,
      [](const double sum, const std::pair<Point3D, double> &rhs)
      {
        return sum + rhs.second;
      }
    );

    std::for_each(
      centers_and_weights.begin(), centers_and_weights.end(),
      [sum_weights](std::pair<Point3D, double> &point)
      {
        point.second /= sum_weights;
      }
    );
  }
}

std::vector<std::pair<Point3D, double>> IcosahedralQuadrature::generate(unsigned int n)
{
  auto faces = icosahedronFaces();
  for (auto i = 0u; i < n; ++i) {
    faces = refine(faces);
    remap(faces);
  }
  auto res = centers_and_weights(faces);
  scale_weights(res);
  return res;
}
