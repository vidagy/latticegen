#ifndef LATTICEGEN_INTEGRATOR_H
#define LATTICEGEN_INTEGRATOR_H

#include <vector>

namespace Math
{
  class IntegratorEquidistant
  {
  public:
    static double simpson(const std::vector<double> &f, double dx);
    static double simpson_3_8(const std::vector<double> &f, double dx);

    static double simpson_alt(const std::vector<double> &f, double dx);
  };
}

#endif //LATTICEGEN_INTEGRATOR_H
