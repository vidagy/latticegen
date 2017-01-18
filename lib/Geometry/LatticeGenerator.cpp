#include "LatticeGenerator.h"

namespace Geometry
{
  namespace
  {
    size_t get_max_steps(const Point3D& v, const Point3D& perp1, const Point3D& perp2, const double max_distance)
    {
      Point3D perp = cross_product(perp1, perp2);
      const double perp_length = perp.length();
      perp = 1.0 / perp_length * perp;

      double perp_component = v * perp;
      return (size_t) (max_distance / perp_component + 1);
    }
  }

  CutoffCube::CutoffCube(const UnitCell3D& unit_cell_, const double a_)
    : Cutoff(unit_cell_), a(a_)
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

  std::tuple<size_t, size_t, size_t> CutoffCube::steps_to_cover() const
  {
    return std::make_tuple(
      get_max_steps(unit_cell.a, unit_cell.b, unit_cell.c, a*sqrt(3)),
      get_max_steps(unit_cell.b, unit_cell.c, unit_cell.a, a*sqrt(3)),
      get_max_steps(unit_cell.c, unit_cell.a, unit_cell.b, a*sqrt(3))
    );
  }

  CutoffSphere::CutoffSphere(const UnitCell3D& unit_cell_, const double r_)
    : Cutoff(unit_cell_), r(r_)
  {
    if (! strictlyPositive(r))
      throw std::invalid_argument("In CutoffSphere::ctor: a must be non-negative but r = " + std::to_string(r));
  }

  bool CutoffSphere::is_included(const Point3D& point) const
  {
    return lessEqualsWithTolerance(point.length(), r);
  }

  std::tuple<size_t, size_t, size_t> CutoffSphere::steps_to_cover() const
  {
    return std::make_tuple(
      get_max_steps(unit_cell.a, unit_cell.b, unit_cell.c, r),
      get_max_steps(unit_cell.b, unit_cell.c, unit_cell.a, r),
      get_max_steps(unit_cell.c, unit_cell.a, unit_cell.b, r)
    );
  }

  CutoffUnitVectors::CutoffUnitVectors(
    const UnitCell3D& unit_cell_,
    const size_t a_max_, const size_t b_max_, const size_t c_max_)
    : Cutoff(unit_cell_)
    , a_max(a_max_), b_max(b_max_), c_max(c_max_)
  {
  }

  bool CutoffUnitVectors::is_included(const Point3D& point) const
  {
    long nx, ny, nz;
    std::tie(nx, ny, nz) = unit_cell.get_offsets(point);
    return (size_t)labs(nx) <= a_max &&
           (size_t)labs(ny) <= b_max &&
           (size_t)labs(nz) <= c_max;
  }

  std::tuple<size_t, size_t, size_t> CutoffUnitVectors::steps_to_cover() const
  {
    return std::make_tuple(a_max, b_max, c_max);
  }

  std::vector<Point3D> LatticeGenerator::generate(const Cutoff& cutoff, const bool positive_only)
  {
    long max_x(0), max_y(0), max_z(0);
    std::tie(max_x, max_y, max_z) = cutoff.steps_to_cover();

    std::vector<Point3D> lattice;
    long min_x(0), min_y(0), min_z(0);
    if (positive_only)
    {
      lattice.reserve((unsigned long) ((max_x + 1) * (max_y + 1) * (max_z + 1)));
    }
    else
    {
      min_x = -max_x;
      min_y = -max_y;
      min_z = -max_z;
      lattice.reserve((unsigned long) ((2 * max_x + 1) * (2 * max_y + 1) * (2 * max_z + 1)));
    }

    const UnitCell3D& unit_cell = cutoff.unit_cell;
    
    for (long n_z = min_z; n_z <= max_z; ++n_z)
    {
      auto c = n_z * unit_cell.c;
      for (long n_y = min_y; n_y <= max_y; ++n_y)
      {
        auto b = n_y * unit_cell.b;
        for (long n_x = min_x; n_x <= max_x; ++n_x)
        {
          auto a = n_x * unit_cell.a;
          auto v = a + b + c;
          if (cutoff.is_included(v))
            lattice.push_back(v);
        } 
      } 
    }
    lattice.shrink_to_fit();
    return lattice;
  }
}
