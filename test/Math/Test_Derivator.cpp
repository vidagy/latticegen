#include <TestUtils/base.h>

#include <Math/Derivator.h>

using namespace Math;

//namespace
//{
//  void print_matrix(const std::vector<std::vector<double>>& matrix) {
//    for (auto i = 0u; i < matrix.size(); ++i) {
//      for (auto j = 0u; j < matrix[i].size(); ++j) {
//        std::cout << std::setw(20) << std::setprecision(17) << std::fixed << matrix[i][j] << "  ";
//      }
//      std::cout << '\n';
//    }
//    std::cout << '\n';
//  }
//}

TEST(TestDerivator, 4)
{
  auto res = Derivator::lagrange_quadrature(4);

  auto reference = std::vector<std::vector<double>>{
    {-1.83333333333333326, 3.00000000000000000,  -1.50000000000000000, 0.33333333333333331},
    {-0.33333333333333331, -0.50000000000000011, 1.00000000000000000,  -0.16666666666666666},
    {0.16666666666666666,  -1.00000000000000000, 0.50000000000000011,  0.33333333333333331},
    {-0.33333333333333331, 1.50000000000000000,  -3.00000000000000000, 1.83333333333333326}
  };

  EXPECT_THAT(res, ::testing::UnorderedElementsAreArray(reference));
}
