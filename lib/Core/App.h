#ifndef LATTICEGEN_APP_H
#define LATTICEGEN_APP_H

#include <mutex>
#include <map>
#include <memory>
#include <string>
#include "Logger.h"

namespace Core
{
  class Logger;

  class App
  {
  public:
    static void Create(int argc_, char const *const *argv_);

    static const std::map<std::string, std::string> &get_args() { return args; }

    static const std::map<std::string, std::string> &get_environment_variables() { return environment_variables; }

    static const Logger &get_logger() { return logger; }

    App() = delete;

    App(App const &) = delete;

    App(App &&) = delete;

    App &operator=(App const &) = delete;

    App &operator=(App &&) = delete;

  private:
    static void create_impl(int argc, char const *const *argv);

    static std::once_flag onceFlag;
    static std::map<std::string, std::string> args;
    static std::map<std::string, std::string> environment_variables;
    static Logger logger;
  };
}

#endif //LATTICEGEN_APP_H