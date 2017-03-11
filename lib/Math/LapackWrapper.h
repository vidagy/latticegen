#ifndef LATTICEGEN_LAPACKWRAPPER_H
#define LATTICEGEN_LAPACKWRAPPER_H

#include <vector>
#include <lapacke.h>

namespace Math
{
  class LapackWrapper
  {
  public:
    /// @brief inverts column ordered matrix in place
    static void invert_matrix(std::vector<double> &matrix, int majority = LAPACK_COL_MAJOR);
  };
}
#endif //LATTICEGEN_LAPACKWRAPPER_H
