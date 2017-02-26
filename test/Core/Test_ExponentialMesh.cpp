#include <gtest/gtest.h>
#include <gmock/gmock.h>

#include <Core/ExponentialMesh.h>

using namespace Core;

TEST(TestExponentialMesh, Simple)
{
  auto mesh = ExponentialMesh(1.0, 0.01, 5);
  auto &points = mesh.points;

  const std::vector<double> reference{
    0,
    0.246262593225747922,
    0.49500016665999963,
    0.746237594267377125,
    1};

  EXPECT_THAT(points, ::testing::ContainerEq(reference));
}
