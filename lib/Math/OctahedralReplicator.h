#ifndef LATTICEGEN_OCTAHEDRALREPLICATOR_H
#define LATTICEGEN_OCTAHEDRALREPLICATOR_H

#include <Core/Point3D.h>
#include <vector>

namespace Math
{
  using namespace Core;

  class OctahedralReplicator
  {
  public:
    typedef std::vector<std::pair<Point3D, double>> Quadrature;

    ///  a1: {1, 0, 0} # = 6
    static void gen_a1(double v, Quadrature &res)
    {
      static const double a = 1.0;
      res.push_back(std::make_pair(Point3D{a, 0.0, 0.0}, v));
      res.push_back(std::make_pair(Point3D{-a, 0.0, 0.0}, v));
      res.push_back(std::make_pair(Point3D{0.0, a, 0.0}, v));
      res.push_back(std::make_pair(Point3D{0.0, -a, 0.0}, v));
      res.push_back(std::make_pair(Point3D{0.0, 0.0, a}, v));
      res.push_back(std::make_pair(Point3D{0.0, 0.0, -a}, v));
    }

    ///  a2: 1/sqrt(2) * {1, 1, 0} # = 12
    static void gen_a2(double v, Quadrature &res)
    {
      static const double a = sqrt(0.5);
      res.push_back(std::make_pair(Point3D{a, a, 0.0}, v));
      res.push_back(std::make_pair(Point3D{-a, a, 0.0}, v));
      res.push_back(std::make_pair(Point3D{a, -a, 0.0}, v));
      res.push_back(std::make_pair(Point3D{-a, -a, 0.0}, v));
      res.push_back(std::make_pair(Point3D{a, 0.0, a}, v));
      res.push_back(std::make_pair(Point3D{-a, 0.0, a}, v));
      res.push_back(std::make_pair(Point3D{a, 0.0, -a}, v));
      res.push_back(std::make_pair(Point3D{-a, 0.0, -a}, v));
      res.push_back(std::make_pair(Point3D{0.0, a, a}, v));
      res.push_back(std::make_pair(Point3D{0.0, -a, a}, v));
      res.push_back(std::make_pair(Point3D{0.0, a, -a}, v));
      res.push_back(std::make_pair(Point3D{0.0, -a, -a}, v));
    }

    ///  a3: 1/sqrt(3) * {1, 1, 1} # = 8
    static void gen_a3(double v, Quadrature &res)
    {
      const static double a = sqrt(1.0 / 3.0);
      res.push_back(std::make_pair(Point3D{a, a, a}, v));
      res.push_back(std::make_pair(Point3D{-a, a, a}, v));
      res.push_back(std::make_pair(Point3D{a, -a, a}, v));
      res.push_back(std::make_pair(Point3D{a, a, -a}, v));
      res.push_back(std::make_pair(Point3D{-a, -a, a}, v));
      res.push_back(std::make_pair(Point3D{-a, a, -a}, v));
      res.push_back(std::make_pair(Point3D{a, -a, -a}, v));
      res.push_back(std::make_pair(Point3D{-a, -a, -a}, v));
    }

    ///  b: {l, l, m} where 2*l^2+m^2 = 1 # = 24
    static void gen_b(double a, double v, Quadrature &res)
    {
      const double b = sqrt(1.0 - 2.0 * a * a);
      res.push_back(std::make_pair(Point3D{b, a, a}, v));
      res.push_back(std::make_pair(Point3D{-b, a, a}, v));
      res.push_back(std::make_pair(Point3D{b, -a, a}, v));
      res.push_back(std::make_pair(Point3D{b, a, -a}, v));
      res.push_back(std::make_pair(Point3D{-b, -a, a}, v));
      res.push_back(std::make_pair(Point3D{-b, a, -a}, v));
      res.push_back(std::make_pair(Point3D{b, -a, -a}, v));
      res.push_back(std::make_pair(Point3D{-b, -a, -a}, v));
      res.push_back(std::make_pair(Point3D{a, b, a}, v));
      res.push_back(std::make_pair(Point3D{-a, b, a}, v));
      res.push_back(std::make_pair(Point3D{a, -b, a}, v));
      res.push_back(std::make_pair(Point3D{a, b, -a}, v));
      res.push_back(std::make_pair(Point3D{-a, -b, a}, v));
      res.push_back(std::make_pair(Point3D{-a, b, -a}, v));
      res.push_back(std::make_pair(Point3D{a, -b, -a}, v));
      res.push_back(std::make_pair(Point3D{-a, -b, -a}, v));
      res.push_back(std::make_pair(Point3D{a, a, b}, v));
      res.push_back(std::make_pair(Point3D{-a, a, b}, v));
      res.push_back(std::make_pair(Point3D{a, -a, b}, v));
      res.push_back(std::make_pair(Point3D{a, a, -b}, v));
      res.push_back(std::make_pair(Point3D{-a, -a, b}, v));
      res.push_back(std::make_pair(Point3D{-a, a, -b}, v));
      res.push_back(std::make_pair(Point3D{a, -a, -b}, v));
      res.push_back(std::make_pair(Point3D{-a, -a, -b}, v));
    }

    ///  c: {p, q, 0} where p^2+q^2 = 1 # = 24
    static void gen_c(double a, double v, Quadrature &res)
    {
      const double b = sqrt(1.0 - a * a);
      res.push_back(std::make_pair(Point3D{a, b, 0.0}, v));
      res.push_back(std::make_pair(Point3D{-a, b, 0.0}, v));
      res.push_back(std::make_pair(Point3D{a, -b, 0.0}, v));
      res.push_back(std::make_pair(Point3D{-a, -b, 0.0}, v));
      res.push_back(std::make_pair(Point3D{b, a, 0.0}, v));
      res.push_back(std::make_pair(Point3D{-b, a, 0.0}, v));
      res.push_back(std::make_pair(Point3D{b, -a, 0.0}, v));
      res.push_back(std::make_pair(Point3D{-b, -a, 0.0}, v));
      res.push_back(std::make_pair(Point3D{a, 0.0, b}, v));
      res.push_back(std::make_pair(Point3D{-a, 0.0, b}, v));
      res.push_back(std::make_pair(Point3D{a, 0.0, -b}, v));
      res.push_back(std::make_pair(Point3D{-a, 0.0, -b}, v));
      res.push_back(std::make_pair(Point3D{b, 0.0, a}, v));
      res.push_back(std::make_pair(Point3D{-b, 0.0, a}, v));
      res.push_back(std::make_pair(Point3D{b, 0.0, -a}, v));
      res.push_back(std::make_pair(Point3D{-b, 0.0, -a}, v));
      res.push_back(std::make_pair(Point3D{0.0, a, b}, v));
      res.push_back(std::make_pair(Point3D{0.0, -a, b}, v));
      res.push_back(std::make_pair(Point3D{0.0, a, -b}, v));
      res.push_back(std::make_pair(Point3D{0.0, -a, -b}, v));
      res.push_back(std::make_pair(Point3D{0.0, b, a}, v));
      res.push_back(std::make_pair(Point3D{0.0, -b, a}, v));
      res.push_back(std::make_pair(Point3D{0.0, b, -a}, v));
      res.push_back(std::make_pair(Point3D{0.0, -b, -a}, v));
    }

    ///  d: {r, S, W} where r^2+S^2+W^2 = 1 # = 48
    static void gen_d(double a, double b, double v, Quadrature &res)
    {
      const double c = sqrt(1.0 - a * a - b * b);
      res.push_back(std::make_pair(Point3D{a, b, c}, v));
      res.push_back(std::make_pair(Point3D{-a, b, c}, v));
      res.push_back(std::make_pair(Point3D{a, -b, c}, v));
      res.push_back(std::make_pair(Point3D{a, b, -c}, v));
      res.push_back(std::make_pair(Point3D{-a, -b, c}, v));
      res.push_back(std::make_pair(Point3D{-a, b, -c}, v));
      res.push_back(std::make_pair(Point3D{a, -b, -c}, v));
      res.push_back(std::make_pair(Point3D{-a, -b, -c}, v));
      res.push_back(std::make_pair(Point3D{a, c, b}, v));
      res.push_back(std::make_pair(Point3D{-a, c, b}, v));
      res.push_back(std::make_pair(Point3D{a, -c, b}, v));
      res.push_back(std::make_pair(Point3D{a, c, -b}, v));
      res.push_back(std::make_pair(Point3D{-a, -c, b}, v));
      res.push_back(std::make_pair(Point3D{-a, c, -b}, v));
      res.push_back(std::make_pair(Point3D{a, -c, -b}, v));
      res.push_back(std::make_pair(Point3D{-a, -c, -b}, v));
      res.push_back(std::make_pair(Point3D{b, a, c}, v));
      res.push_back(std::make_pair(Point3D{-b, a, c}, v));
      res.push_back(std::make_pair(Point3D{b, -a, c}, v));
      res.push_back(std::make_pair(Point3D{b, a, -c}, v));
      res.push_back(std::make_pair(Point3D{-b, -a, c}, v));
      res.push_back(std::make_pair(Point3D{-b, a, -c}, v));
      res.push_back(std::make_pair(Point3D{b, -a, -c}, v));
      res.push_back(std::make_pair(Point3D{-b, -a, -c}, v));
      res.push_back(std::make_pair(Point3D{b, c, a}, v));
      res.push_back(std::make_pair(Point3D{-b, c, a}, v));
      res.push_back(std::make_pair(Point3D{b, -c, a}, v));
      res.push_back(std::make_pair(Point3D{b, c, -a}, v));
      res.push_back(std::make_pair(Point3D{-b, -c, a}, v));
      res.push_back(std::make_pair(Point3D{-b, c, -a}, v));
      res.push_back(std::make_pair(Point3D{b, -c, -a}, v));
      res.push_back(std::make_pair(Point3D{-b, -c, -a}, v));
      res.push_back(std::make_pair(Point3D{c, a, b}, v));
      res.push_back(std::make_pair(Point3D{-c, a, b}, v));
      res.push_back(std::make_pair(Point3D{c, -a, b}, v));
      res.push_back(std::make_pair(Point3D{c, a, -b}, v));
      res.push_back(std::make_pair(Point3D{-c, -a, b}, v));
      res.push_back(std::make_pair(Point3D{-c, a, -b}, v));
      res.push_back(std::make_pair(Point3D{c, -a, -b}, v));
      res.push_back(std::make_pair(Point3D{-c, -a, -b}, v));
      res.push_back(std::make_pair(Point3D{c, b, a}, v));
      res.push_back(std::make_pair(Point3D{-c, b, a}, v));
      res.push_back(std::make_pair(Point3D{c, -b, a}, v));
      res.push_back(std::make_pair(Point3D{c, b, -a}, v));
      res.push_back(std::make_pair(Point3D{-c, -b, a}, v));
      res.push_back(std::make_pair(Point3D{-c, b, -a}, v));
      res.push_back(std::make_pair(Point3D{c, -b, -a}, v));
      res.push_back(std::make_pair(Point3D{-c, -b, -a}, v));
    }
  };
}

#endif //LATTICEGEN_OCTAHEDRALREPLICATOR_H
