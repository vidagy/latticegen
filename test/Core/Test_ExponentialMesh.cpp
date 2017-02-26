#include <gtest/gtest.h>
#include <gmock/gmock.h>

#include <Core/ExponentialMesh.h>

using namespace Core;

TEST(TestExponentialMesh, Simple)
{
  auto mesh = ExponentialMesh(0.01, 1.0, 5);
  auto &points = mesh.points;

  // for(auto p: points)
  //   std::cout << std::setprecision(18) << p << std::endl;

  const std::vector<double> reference{
    0.0100000000000000002,
    0.0316227766016837983,
    0.100000000000000019,
    0.316227766016838052,
    1.00000000000000044};

  EXPECT_THAT(points, ::testing::ContainerEq(reference));
}
