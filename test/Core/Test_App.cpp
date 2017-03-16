#include <gtest/gtest.h>
#include <gmock/gmock.h>

#include <Core/App.h>

using namespace Core;

TEST(TestApp, OnlyOnce)
{
  const char *a = "a";
  const char *b = "--log_level=debug";
  const char *const ab[2] = {a, b};

  App::Create(2, ab);

  auto args = App::get_args();
  EXPECT_EQ(args["LOG_LEVEL"], "debug");

  auto environ = App::get_environment_variables();
  EXPECT_TRUE(environ.find("PATH") != environ.end());
}