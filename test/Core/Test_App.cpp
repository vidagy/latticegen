#include <TestUtils/base.h>

#include <Core/App.h>

using namespace Core;

const char *a = "a";
const char *b = "--log_level=debug";
const char *const ab[2] = {a, b};

TEST(TestApp, OnlyOnce)
{
  App::initialize(2, ab);
  auto args = App::get_args();
  EXPECT_EQ(args["LOG_LEVEL"], "info");

  auto environ = App::get_environment_variables();
  EXPECT_TRUE(environ.find("PATH") != environ.end());
}