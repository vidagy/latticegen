#include "Derivator.h"

#include <Math/Factorial.h>
#include <stdexcept>

using namespace Math;

namespace
{
  std::vector<double> generate_A(int n)
  {
    std::vector<double> A;
    A.reserve(n);

    A[0] = factorial(n - 1);
    for (auto i = 1; i < n; ++i) {
      A[i] = A[i - 1] * i / (n - i);
    }
    return A;
  }
}

std::vector<std::vector<double>> Derivator::lagrange_quadrature(int n)
{
  if (n < 1)
    throw std::invalid_argument("in Derivator::lagrange_quadrature n = " + std::to_string(n));

  std::vector<std::vector<double>> result(n, std::vector<double>(n, 0.0));

  auto A = generate_A(n);
  for (auto i = 0; i < (n + 1) / 2; ++i) {
    auto sum = 0.0;
    for (auto j = 0; j < n; ++j) {
      if (i != j) {
        auto result_ij = ((i - j + 1) % 2 ? -1.0 : 1.0) / (j - i) * A[i] / A[j];
        result[i][j] = result_ij;
        sum += result_ij;
      }
    }
    result[i][i] = -sum;
  }

  for (auto i = (n + 1) / 2; i < n; ++i) {
    for (auto j = 0; j < n; ++j) {
      result[i][j] = -result[n - 1 - i][n - 1 - j];
    }
  }
  return result;
}