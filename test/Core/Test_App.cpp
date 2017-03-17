#include <gtest/gtest.h>
#include <gmock/gmock.h>

#include <Core/App.h>

using namespace Core;

const char *a = "a";
const char *b = "--log_level=debug";
const char *const ab[2] = {a, b};

TEST(TestApp, OnlyOnce)
{
  App::Create(2, ab);

  auto args = App::get_args();
  EXPECT_EQ(args["LOG_LEVEL"], "debug");

  auto environ = App::get_environment_variables();
  EXPECT_TRUE(environ.find("PATH") != environ.end());
}

TEST(TestApp, Logger)
{
  // redundant, just to make sure it is working if run withour only once
  App::Create(2, ab);

  App::get_logger().log(Logger::Severity::Error, "this is an error message");
  App::get_logger().log(Logger::Severity::Warning, "this is a warning message");
  App::get_logger().log(Logger::Severity::Info, "this is an info message");
  App::get_logger().log(Logger::Severity::Debug, "this is a debug message");
  App::get_logger().flush();
}