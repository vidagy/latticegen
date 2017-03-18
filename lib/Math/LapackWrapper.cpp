#include <cmath>
#include <Core/Exceptions.h>
#include "LapackWrapper.h"

using namespace Math;

void LapackWrapper::invert_matrix(std::vector<double> &matrix, int majority)
{
  unsigned int n = static_cast<unsigned int>(lround(sqrt(matrix.size())));
  if (n * n != matrix.size())
    THROW_INVALID_ARGUMENT("matrix not square in invert_matrix length = " + std::to_string(matrix.size()));

  std::vector<int> pivotArray;
  pivotArray.reserve(n);

  int errorHandler = LAPACKE_dgetrf(majority, n, n, &(matrix[0]), n, &(pivotArray[0]));
  if (errorHandler)
    std::logic_error("Error in LAPACKE_dgetrf. Error code = " + std::to_string(errorHandler));

  std::vector<double> lapackWorkspace;
  lapackWorkspace.reserve(n * n);
  errorHandler = LAPACKE_dgetri(majority, n, &(matrix[0]), n, &(pivotArray[0]));
  if (errorHandler)
    std::logic_error("Error in LAPACKE_dgetri. Error code = " + std::to_string(errorHandler));
}