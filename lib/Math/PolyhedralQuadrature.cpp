#include "PolyhedralQuadrature.h"
#include <numeric>
#include <algorithm>

using namespace Math;

namespace
{
  template<typename FaceT>
  struct Face
  {
    virtual void remap() = 0;

    virtual void refine(std::vector<FaceT> &range) const = 0;

    virtual std::pair<Point3D, double> center_and_weight() const = 0;
  };

  struct Triangle : public Face<Triangle>
  {
    Triangle(const Point3D &a_, const Point3D &b_, const Point3D &c_) : a(a_), b(b_), c(c_) {}

    virtual void remap() override final
    {
      a /= a.length();
      b /= b.length();
      c /= c.length();
    };

    virtual void refine(std::vector<Triangle> &range) const override
    {
      auto ab = (a + b) / 2.0;
      auto ac = (a + c) / 2.0;
      auto bc = (b + c) / 2.0;

      range.push_back(Triangle{a, ab, ac});
      range.push_back(Triangle{b, ab, bc});
      range.push_back(Triangle{c, ac, bc});
      range.push_back(Triangle{ab, ac, bc});
    }

    virtual std::pair<Point3D, double> center_and_weight() const override
    {
      auto p = (a + b + c) / 3.0;
      auto w = cross_product(b - a, c - a).length() / 2.0;
      return std::make_pair(p / p.length(), w);
    }

    Point3D a;
    Point3D b;
    Point3D c;
  };

  struct Rectangle : public Face<Rectangle>
  {
    Rectangle(const Point3D &a_, const Point3D &b_, const Point3D &c_, const Point3D &d_)
      : a(a_), b(b_), c(c_), d(d_) {}

    virtual void remap() override final
    {
      a /= a.length();
      b /= b.length();
      c /= c.length();
      d /= d.length();
    };

    virtual void refine(std::vector<Rectangle> &range) const override
    {
      auto ab = (a + b) / 2.0;
      auto bc = (b + c) / 2.0;
      auto cd = (c + d) / 2.0;
      auto da = (d + a) / 2.0;

      auto abcd = (a + b + c + d) / 4.0;

      range.push_back(Rectangle{a, ab, abcd, da});
      range.push_back(Rectangle{ab, b, bc, abcd});
      range.push_back(Rectangle{abcd, bc, c, cd});
      range.push_back(Rectangle{da, abcd, cd, d});
    }

    virtual std::pair<Point3D, double> center_and_weight() const override
    {
      auto p = (a + b + c + d) / 4.0;
      auto w = cross_product(c - a, b - d).length() / 2.0;
      return std::make_pair(p / p.length(), w);
    }

    Point3D a;
    Point3D b;
    Point3D c;
    Point3D d;
  };

  template<typename FaceT>
  void remap(std::vector<FaceT> &faces)
  {
    std::for_each(
      faces.begin(), faces.end(),
      [](FaceT &face)
      {
        return face.remap();
      }
    );
  }

  template<typename FaceT>
  std::vector<FaceT> refine(std::vector<FaceT> &faces)
  {
    auto res = std::vector<FaceT>();
    res.reserve(faces.size() * 4);
    std::for_each(
      faces.begin(), faces.end(),
      [&res](const FaceT &face)
      {
        face.refine(res);
      }
    );
    return res;
  }

  template<typename FaceT>
  Quadrature centers_and_weights(const std::vector<FaceT> &faces)
  {
    auto res = Quadrature();
    res.reserve(faces.size());

    std::transform(
      faces.begin(), faces.end(), std::back_inserter(res),
      [](const FaceT &face)
      {
        return face.center_and_weight();
      }
    );

    return res;
  };

  void scale_weights(Quadrature &centers_and_weights)
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

namespace
{
  std::vector<Triangle> icosahedron_faces()
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
}

Quadrature IcosahedralQuadrature::generate(unsigned int n)
{
  auto faces = icosahedron_faces();
  for (auto i = 0u; i < n; ++i) {
    faces = refine(faces);
    remap(faces);
  }
  auto res = centers_and_weights(faces);
  scale_weights(res);
  return res;
}

namespace
{
  std::vector<Rectangle> cubic_faces()
  {
    static const auto a1 = Point3D{1, 1, 1};
    static const auto a2 = Point3D{-1, 1, 1};
    static const auto a3 = Point3D{-1, -1, 1};
    static const auto a4 = Point3D{1, -1, 1};
    static const auto b1 = Point3D{1, 1, -1};
    static const auto b2 = Point3D{-1, 1, -1};
    static const auto b3 = Point3D{-1, -1, -1};
    static const auto b4 = Point3D{1, -1, -1};

    auto faces = std::vector<Rectangle>{
      {a1, a2, a3, a4},

      {a1, a2, b2, b1},
      {a2, a3, b3, b2},
      {a3, a4, b4, b3},
      {a4, a1, b1, b4},

      {b1, b2, b3, b4}
    };

    remap(faces);
    return faces;
  }
}

Quadrature OctahedralQuadrature::generate(unsigned int n)
{
  auto faces = cubic_faces();
  for (auto i = 0u; i < n; ++i) {
    faces = refine(faces);
    remap(faces);
  }
  auto res = centers_and_weights(faces);
  scale_weights(res);
  return res;
}
