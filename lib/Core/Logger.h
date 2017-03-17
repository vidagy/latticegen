#ifndef LATTICEGEN_LOGGER_H
#define LATTICEGEN_LOGGER_H

#include <memory>

namespace Core
{
  class LoggerImpl;

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

    void log(Severity severity, const char *message) const;

    void flush() const;

  private:
    Logger();

    ~Logger();

    void initialize();

    std::unique_ptr<LoggerImpl> impl;

    friend class App;
  };
}

#endif //LATTICEGEN_LOGGER_H
