#include "Logger.h"
#include "App.h"

#include <boost/log/common.hpp>
#include <boost/log/sinks.hpp>
#include <boost/core/null_deleter.hpp>

#include <fstream>

using namespace Core;
using namespace boost::log;

namespace
{
  class LoggerImpl
  {
  public:
    LoggerImpl(Logger::Severity severity, const std::string &target_string)
    {
      typedef sinks::asynchronous_sink<sinks::text_ostream_backend> text_sink;
      boost::shared_ptr<text_sink> sink = boost::make_shared<text_sink>();

      boost::shared_ptr<std::ostream> stream{&std::cout, boost::null_deleter{}};
      sink->locked_backend()->add_stream(stream);
      sink->set_filter(
        [severity](const attribute_value_set &set)
        {
          return set["Severity"].extract<Logger::Severity>() >= severity;
        });

      core::get()->add_sink(sink);
    }

    void log(Logger::Severity severity, const char *message)
    {
      BOOST_LOG_SEV(lg, severity) << message;
    }

    void flush() const
    {
      core::get()->flush();
    }

    sources::severity_logger<Logger::Severity> lg;
  };

  static std::unique_ptr<LoggerImpl> loggerImpl;
}

namespace
{
  Logger::Severity parse_severity(const std::string &s)
  {
    if (s == "DEBUG")
      return Logger::Severity::Debug;
    if (s == "INFO")
      return Logger::Severity::Info;
    if (s == "WARNING")
      return Logger::Severity::Warning;
    if (s == "ERROR")
      return Logger::Severity::Error;

    return Logger::Severity::Info;
  }

//  std::ostream parse_target(const std::string& s) {
//    if (s == "COUT")
//      return std::cout;
//    else {
//      std::ofstream f(s,std::fstream::out | std::fstream::app);
//      if (!f.is_open()) {
//        throw std::invalid_argument("Could not open file " + s + " as log output");
//      }
//      return f;
//    }
//  }

  std::string get_setting(
    const std::map<std::string, std::string> &args,
    const std::map<std::string, std::string> &env,
    const std::string &key
  )
  {
    auto arg = args.find(key);
    auto envv = env.find(key);
    if (arg != args.end()) {
      auto value = arg->second;
      std::transform(value.begin(), value.end(), value.begin(), ::toupper);
      return value;
    } else if (envv != env.end()) {
      auto value = envv->second;
      std::transform(value.begin(), value.end(), value.begin(), ::toupper);
      return value;
    } else
      return std::string();
  }
}

void Logger::initialize()
{
  auto args = App::get_args();
  auto env = App::get_environment_variables();
  auto severity = parse_severity(get_setting(args, env, "LOG_LEVEL"));
  auto target_string = get_setting(args, env, "LOG_TARGET");

  ::loggerImpl = std::make_unique<LoggerImpl>(severity, target_string);
}

void Logger::log(Logger::Severity severity, const char *message)
{
  ::loggerImpl->log(severity, message);
}

void Logger::flush()
{
  ::loggerImpl->flush();
}
