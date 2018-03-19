#include <cmath>
#include <stdexcept>
#include <iostream>
#include <iomanip>
#include <fstream>
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

  void log(const std::vector<double> &R, const std::string &filename)
  {
    std::ofstream out_R;
    out_R.open(filename + ".dat");
    for (auto r : R)
      out_R << std::setw(20) << std::setprecision(17) << std::fixed << r << "\n";
    out_R.close();
  }

  void log(const std::vector<std::complex<double>> &R, const std::string &filename)
  {
    std::ofstream out_R;
    out_R.open(filename + ".dat");
    for (auto r : R)
      out_R << std::setw(20) << std::setprecision(17) << std::fixed << r.real() << "\t" << r.imag() << "\n";
    out_R.close();
  }

  void log(const std::vector<Core::Point3D> &R, const std::string &filename)
  {
    std::ofstream out_R;
    out_R.open(filename + ".dat");
    for (auto p : R)
      out_R << std::setfill(' ') << std::setw(20) << std::setprecision(17) << std::fixed
            << p(0) << '\t' << p(1) << '\t' << p(2) << "\n";
    out_R.close();
  }
}