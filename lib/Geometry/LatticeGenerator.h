#ifndef LATTICEGEN_LATTICEGENERATOR_H_
#define LATTICEGEN_LATTICEGENERATOR_H_

#include "Point3D.h"
#include "Cutoff.h"

#include <vector>
#include <memory>

namespace Geometry
{
  class LatticeGenerator
  {
  public:
    LatticeGenerator(std::unique_ptr<Cutoff> cutoff);

    std::vector<Point3D> generate(const bool positive_only = false);

    std::unique_ptr<Cutoff> cutoff;
  };
}

#endif // LATTICEGEN_LATTICEGENERATOR_H_
