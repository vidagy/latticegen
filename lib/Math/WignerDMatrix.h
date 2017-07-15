#ifndef LATTICEGEN_WIGNERDMATRIX_H
#define LATTICEGEN_WIGNERDMATRIX_H

#include <complex>
#include <Eigen/Dense>

namespace Math
{
  class WignerDMatrix
  {
  public:
    /// @brief as defined on Wikipedia.
    /// Using the z-y-z convention, right-handed frame, right-hand screw rule, active interpretation
    /// Y_{lm'}(r) = \sum_{m=-l}^{l} d_{m'm} * Y_{lm}(R*r)
    static Eigen::MatrixXcd calculate(double alpha, double beta, double gamma, unsigned int l);

    /// @brief as defined on Wikipedia.
    /// Using the z-y-z convention, right-handed frame, right-hand screw rule, active interpretation
    /// Y_{lm'}(r) = \sum_{m=-l}^{l} d_{m'm} * Y_{lm}(R*r)
    static Eigen::MatrixXd calculate_small(double beta, unsigned int l);
  };
}


#endif //LATTICEGEN_WIGNERDMATRIX_H
