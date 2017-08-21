#ifndef LATTICEGEN_INTEGRATOR_H
#define LATTICEGEN_INTEGRATOR_H

#include <vector>
#include <Core/RadialMesh.h>

namespace Math
{
  class IntegratorEquidistant
  {
  public:
    static double simpson(const std::vector<double> &f, double dx);
    static double simpson_3_8(const std::vector<double> &f, double dx);
    static double simpson_alt(const std::vector<double> &f, double dx);

    static double trapezoidal(const std::vector<double> &f, double dx, unsigned int quadrature = 7);
  };

  class IntegratorExponential
  {
  public:
    static double simpson(const std::vector<double> &f, const Core::ExponentialMesh &x);
    static double simpson_alt(const std::vector<double> &f, const Core::ExponentialMesh &x);
  };

  class IntegratorGeneric
  {
  public:
    static double integrate(const std::vector<double> &f, const std::vector<double> &x);
  };
}

#endif //LATTICEGEN_INTEGRATOR_H
