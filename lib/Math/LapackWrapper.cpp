#include <cmath>
#include <Core/Exceptions.h>
#include "LapackWrapper.h"
#include <lapacke.h>

using namespace Math;

namespace
{
  int get_majority(LapackWrapper::Majority majority)
  {
    switch (majority) {
      case LapackWrapper::Majority::Column:
        return LAPACK_COL_MAJOR;
      case LapackWrapper::Majority::Row:
        return LAPACK_ROW_MAJOR;
      default:
        THROW_LOGIC_ERROR("get_majority not implemented for " + std::to_string(majority));
    }
  }
}

void LapackWrapper::invert_matrix(std::vector<double> &matrix, Majority majority)
{
  auto lapack_majority = get_majority(majority);
  unsigned int n = static_cast<unsigned int>(lround(sqrt(matrix.size())));
  if (n * n != matrix.size())
    THROW_INVALID_ARGUMENT("matrix not square in invert_matrix length = " + std::to_string(matrix.size()));

  std::vector<int> pivotArray;
  pivotArray.reserve(n);

  int errorHandler = LAPACKE_dgetrf(lapack_majority, n, n, &(matrix[0]), n, &(pivotArray[0]));
  if (errorHandler)
    std::logic_error("Error in LAPACKE_dgetrf. Error code = " + std::to_string(errorHandler));

  std::vector<double> lapackWorkspace;
  lapackWorkspace.reserve(n * n);
  errorHandler = LAPACKE_dgetri(lapack_majority, n, &(matrix[0]), n, &(pivotArray[0]));
  if (errorHandler)
    std::logic_error("Error in LAPACKE_dgetri. Error code = " + std::to_string(errorHandler));
}