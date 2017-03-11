#include <cmath>
#include <stdexcept>
#include <iostream>
#include <iomanip>
#include "Utils.h"

namespace Utils
{
  void print_square_matrix(const std::vector<double> &matrix, bool is_row_ordered)
  {
    unsigned int n = static_cast<unsigned int>(lround(sqrt(matrix.size())));
    if (n * n != matrix.size())
      throw std::invalid_argument("matrix not square in print_square_matrix size = " + std::to_string(matrix.size()));
    for (auto i = 0u; i < n; ++i) {
      for (auto j = 0u; j < n; ++j) {
        int index;
        if (is_row_ordered)
          index = i * n + j;
        else
          index = i + j * n;
        std::cout << std::setw(20) << std::setprecision(17) << std::fixed << matrix[index] << "  ";
      }
      std::cout << '\n';
    }
    std::cout << '\n';
  }

  void print_matrix(const std::vector<std::vector<double>> &matrix)
  {
    for (auto i = 0u; i < matrix.size(); ++i) {
      for (auto j = 0u; j < matrix[i].size(); ++j) {
        std::cout << std::setw(20) << std::setprecision(17) << std::fixed << matrix[i][j] << "  ";
      }
      std::cout << '\n';
    }
    std::cout << '\n';
  }
}