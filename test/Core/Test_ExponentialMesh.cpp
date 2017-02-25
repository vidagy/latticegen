#include <gtest/gtest.h>
#include <gmock/gmock.h>

#include <Core/ExponentialMesh.h>

using namespace Core;

TEST(TestExponentialMesh, Simple)
{
  auto mesh = ExponentialMesh(1.0, 0.01, 5);
  auto &points = mesh.points;

  const std::vector<double> reference{0.1960201320927419388,
                                      0.3940102992648801816,
                                      0.5939903006981216427,
                                      0.7959801345592550925,
                                      1};

  EXPECT_THAT(points, ::testing::ContainerEq(reference));
}
