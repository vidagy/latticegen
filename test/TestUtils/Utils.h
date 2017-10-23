#ifndef LATTICEGEN_UTILS_H
#define LATTICEGEN_UTILS_H

#include <vector>
#include <complex>
#include <Core/Point3D.h>

namespace Utils
{
  // TODO: add [[maybe_unused]] when on c++17
  void print_square_matrix(const std::vector<double> &matrix, bool is_row_ordered);

  // TODO: add [[maybe_unused]] when on c++17
  void print_matrix(const std::vector<std::vector<double>> &matrix);

  void log(const std::vector<double> &R, const std::string &filename);

  void log(const std::vector<std::complex<double>> &R, const std::string &filename);

  void log(const std::vector<Core::Point3D> &R, const std::string &filename);
}

#endif //LATTICEGEN_UTILS_H
