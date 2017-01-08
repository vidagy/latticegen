#include "LatticeGenerator.h"

#include <tuple>

namespace Geometry
{
  LatticeGenerator::LatticeGenerator(std::unique_ptr<Cutoff> cutoff_)
    : cutoff(std::move(cutoff_))
  {
  } 

  std::vector<Point3D> LatticeGenerator::generate(const bool positive_only)
  {
    long max_x, max_y, max_z;
    std::tie(max_x, max_y, max_z) = cutoff->steps_to_cover();

    std::vector<Point3D> lattice;
    long n_x(0), n_y(0), n_z(0);
    if (positive_only)
    {
      lattice.reserve( (max_x+1) *(max_y+1) * (max_z+1));
    }
    else
    {
      n_x = -max_x; 
      n_y = -max_y;
      n_z = -max_z;
      lattice.reserve( (2*max_x+1) *(2*max_y+1) * (2*max_z+1));
    }

    const UnitCell3D& unit_cell = cutoff->unit_cell;
    
    for (; n_z <= max_z; ++n_z)
    {
      auto c = n_z * unit_cell.c;
      for (; n_y <= max_y; ++n_y)
      {
        auto b = n_y * unit_cell.b;
        for (; n_x < max_x; ++n_x)
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