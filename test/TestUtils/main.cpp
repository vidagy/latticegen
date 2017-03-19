#include <gtest/gtest.h>
#include <gmock/gmock.h>
#include <Core/App.h>

int main(int argc, char **argv)
{
  // initialize Google Test framework
  ::testing::InitGoogleTest(&argc, argv);

  // initialize App
  Core::App::initialize(argc, static_cast<char const *const *>(argv));

  return RUN_ALL_TESTS();
}
