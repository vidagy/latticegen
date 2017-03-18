#ifndef LATTICEGEN_LOGGER_H
#define LATTICEGEN_LOGGER_H

#include <memory>
#include <stdexcept>

namespace Core
{
  class Logger
  {
  public:
    enum Severity
    {
      Debug = 0,
      Info = 1,
      Warning = 2,
      Error = 3
    };

    static void log(Severity severity, const char *message);

    static void
    log(Severity severity, const char *filename, int line, const char *function_name, const std::string &message);

    static void flush();

    Logger() = delete;

    Logger(const Logger &) = delete;

    Logger(Logger &&) = delete;

    Logger &operator=(const Logger &) = delete;

    Logger &operator=(const Logger &&) = delete;

  private:
    static void initialize();
    friend class App;
  };
}

#endif //LATTICEGEN_LOGGER_H
