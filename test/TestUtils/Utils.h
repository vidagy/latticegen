#ifndef LATTICEGEN_UTILS_H
#define LATTICEGEN_UTILS_H

#include <vector>
#include <Core/Point3D.h>

namespace Utils
{
  void print_square_matrix(const std::vector<double> &matrix, bool is_row_ordered);

  void print_matrix(const std::vector<std::vector<double>> &matrix);

  void log(const std::vector<double> &R, const std::string &filename);

  void log(const std::vector<Core::Point3D> &R, const std::string &filename);
}

#endif //LATTICEGEN_UTILS_H
