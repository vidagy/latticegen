#include "LatticeGenerator.h"

namespace Geometry
{
  LatticeGenerator::LatticeGenerator(std::shared_ptr<Cutoff> cutoff_)
    : cutoff(cutoff_)
  {
  } 

  std::vector<Point3D> LatticeGenerator::generate(const bool positive_only)
  {
    long max_x(0), max_y(0), max_z(0);
    std::tie(max_x, max_y, max_z) = cutoff->steps_to_cover();

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

    const UnitCell3D& unit_cell = cutoff->unit_cell;
    
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
          if (cutoff->is_included(v))
            lattice.push_back(v);
        } 
      } 
    }
    lattice.shrink_to_fit();
    return lattice;
  }
}