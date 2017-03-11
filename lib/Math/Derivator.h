#ifndef LATTICEGEN_DERIVATOR_H
#define LATTICEGEN_DERIVATOR_H

#include <vector>

namespace Math
{
  class Derivator
  {
  public:
    /// Generates the nth order Lagrangian quadrature matrix [m]
    /// (f')_i = sum(j=0..n-1) ( m_ij * f_j )
    /// where f_i = f(x_i) and (f')_i = f'(x_i)
    /// and x_i is an equidistant grid
    static std::vector<std::vector<double>> lagrange_quadrature(int n);
  };
}

#endif //LATTICEGEN_DERIVATOR_H
