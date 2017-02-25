#ifndef LATTICEGEN_INTEGRATOR_H
#define LATTICEGEN_INTEGRATOR_H

#include <vector>

namespace Math
{
  class IntegratorEquidistant
  {
  public:
    static double simpson(const std::vector<double> &f, double dx);
  };
}

#endif //LATTICEGEN_INTEGRATOR_H
