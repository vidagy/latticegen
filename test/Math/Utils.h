#ifndef LATTICEGEN_UTILS_H
#define LATTICEGEN_UTILS_H

#include <vector>

namespace Utils
{
  void print_square_matrix(const std::vector<double> &matrix, bool is_row_ordered);

  void print_matrix(const std::vector<std::vector<double>> &matrix);
}
#endif //LATTICEGEN_UTILS_H
