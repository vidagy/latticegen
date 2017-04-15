#ifndef LATTICEGEN_LAPACKWRAPPER_H
#define LATTICEGEN_LAPACKWRAPPER_H

#include <vector>

namespace Math
{
  class LapackWrapper
  {
  public:
    enum Majority
    {
      Column,
      Row
    };

    /// @brief inverts column ordered matrix in place
    static void invert_matrix(std::vector<double> &matrix, Majority majority = Majority::Column);
  };
}
#endif //LATTICEGEN_LAPACKWRAPPER_H
