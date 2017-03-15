#ifndef LATTICEGEN_APP_H
#define LATTICEGEN_APP_H

#include <mutex>
#include <map>

namespace Core
{
  class App
  {
  public:
    static void Create(int argc_, char const *const *argv_);

    App() = delete;

    App(App const &) = delete;

    App(App &&) = delete;

    App &operator=(App const &) = delete;

    App &operator=(App &&) = delete;

  private:
    static void set_up_args(int argc, char const *const *argv);

    static std::once_flag onceFlag;
    static std::map<std::string, std::string> args;
  };
}

#endif //LATTICEGEN_APP_H
