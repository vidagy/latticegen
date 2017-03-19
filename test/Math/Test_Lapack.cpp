#include <TestUtils/base.h>

#include <Math/LapackWrapper.h>

using namespace Math;

TEST(TestLapack, InvertMatrixColumnMajor)
{
  std::vector<double> matrix = {
    1, 4, 7,
    2, 1, 8,
    3, 6, 1
  };

  std::vector<double> reference = {
    -47.0 / 104.0, 19.0 / 52.0, 25.0 / 104.0,
    11.0 / 52.0, -5.0 / 26.0, 3.0 / 52.0,
    9.0 / 104.0, 3.0 / 52.0, -7.0 / 104.0
  };

  LapackWrapper::invert_matrix(matrix);
  EXPECT_THAT(matrix, ::testing::Pointwise(NearWithTolerance(std::numeric_limits<double>::epsilon()), reference));
}

TEST(TestLapack, InvertMatrixRowMajor)
{
  std::vector<double> matrix = {
    1, 2, 3,
    4, 1, 6,
    7, 8, 1
  };

  std::vector<double> reference = {
    -47.0 / 104.0, 11.0 / 52.0, 9.0 / 104.0,
    19.0 / 52.0, -5.0 / 26.0, 3.0 / 52.0,
    25.0 / 104.0, 3.0 / 52.0, -7.0 / 104.0
  };

  LapackWrapper::invert_matrix(matrix, LAPACK_ROW_MAJOR);
  EXPECT_THAT(matrix, ::testing::Pointwise(NearWithTolerance(std::numeric_limits<double>::epsilon()), reference));
}