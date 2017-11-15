#include <TestUtils/base.h>
#include <Core/Mute.h>

LATTICEGEN_MUTE_BEGIN
LATTICEGEN_MUTE_EIGEN
#include <Eigen/Dense>

LATTICEGEN_MUTE_END

TEST(TestEigen, First)
{
  Eigen::MatrixXd m(3, 3);
  m <<
    1, 4, 7,
    2, 1, 8,
    3, 6, 1;

  std::cout << m.size();

  Eigen::MatrixXd reference(3, 3);
  reference <<
            -47.0 / 104.0, 19.0 / 52.0, 25.0 / 104.0,
    11.0 / 52.0, -5.0 / 26.0, 3.0 / 52.0,
    9.0 / 104.0, 3.0 / 52.0, -7.0 / 104.0;

  auto inverse = m.inverse().eval();

  EXPECT_DOUBLE_EQ(reference(0, 0), inverse(0, 0));
  EXPECT_DOUBLE_EQ(reference(0, 1), inverse(0, 1));
  EXPECT_DOUBLE_EQ(reference(0, 2), inverse(0, 2));
  EXPECT_DOUBLE_EQ(reference(1, 0), inverse(1, 0));
  EXPECT_DOUBLE_EQ(reference(1, 1), inverse(1, 1));
  EXPECT_DOUBLE_EQ(reference(1, 2), inverse(1, 2));
  EXPECT_DOUBLE_EQ(reference(2, 0), inverse(2, 0));
  EXPECT_DOUBLE_EQ(reference(2, 1), inverse(2, 1));
  EXPECT_DOUBLE_EQ(reference(2, 2), inverse(2, 2));
}
