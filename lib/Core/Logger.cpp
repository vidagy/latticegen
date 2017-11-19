#include "Logger.h"
#include "App.h"

#include <boost/log/common.hpp>
#include <boost/log/sinks.hpp>
#include <boost/log/expressions.hpp>
#include <boost/log/utility/setup/file.hpp>
#include <boost/log/utility/setup/common_attributes.hpp>
#include <boost/log/support/date_time.hpp>
#include <boost/core/null_deleter.hpp>

#include <fstream>

using namespace Core;

namespace
{
  using namespace boost::log;
  const std::string default_filename = "latticegen.log";

  class LoggerImpl
  {
  public:
    LoggerImpl(Logger::Level severity, Logger::Target target, const std::string &target_string)
    {
      if (target == Logger::Target::cout) {
        typedef sinks::asynchronous_sink<sinks::text_ostream_backend> text_sink;
        boost::shared_ptr<text_sink> sink = boost::make_shared<text_sink>();

        boost::shared_ptr<std::ostream> stream{&std::cout, boost::null_deleter{}};
        sink->locked_backend()->add_stream(stream);
        sink->set_filter(
          [severity](const attribute_value_set &set)
          {
            return set["Severity"].extract<Logger::Level>() >= severity;
          });

        core::get()->add_sink(sink);
        add_common_attributes();
      } else {
        std::string filename;
        if (target == Logger::Target::default_file)
          filename = default_filename;
        else
          filename = target_string;

        add_file_log
          (
            keywords::file_name = filename,
            keywords::rotation_size = 10 * 1024 * 1024,
            keywords::format =
              (
                expressions::stream
                  << '['
                  << expressions::format_date_time<boost::posix_time::ptime>("TimeStamp", "%Y%m%d %H:%M:%S")
                  << "] - "
                  << expressions::smessage
              )
          );

        add_common_attributes();
      }
    }

    void flush() const
    {
      core::get()->flush();
    }

    sources::severity_logger<Logger::Level> lg;
  };

  std::unique_ptr<LoggerImpl> loggerImpl;
}

namespace
{
  Logger::Level parse_severity(const std::string &s)
  {
    if (s == "DEBUG")
      return Logger::Level::Debug;
    if (s == "INFO")
      return Logger::Level::Info;
    if (s == "WARNING")
      return Logger::Level::Warning;
    if (s == "ERROR")
      return Logger::Level::Error;

    return Logger::Level::Info;
  }

  Logger::Target parse_target_type(const std::string &s)
  {
    auto caps = s;
    std::transform(caps.begin(), caps.end(), caps.begin(), ::toupper);
    if (caps == "COUT")
      return Logger::Target::cout;
    if (caps == "DEFAULT")
      return Logger::Target::default_file;
    else
      return Logger::Target::file;
  }

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
  auto target_type = parse_target_type(target_string);

  ::loggerImpl = std::make_unique<LoggerImpl>(severity, target_type, target_string);
}

void Logger::log(Logger::Level severity, const std::string &message)
{
  BOOST_LOG_SEV(loggerImpl->lg, severity) << message;
}

void Logger::log(Core::Logger::Level severity, const char *filename, int line, const char *function_name,
                 const std::string &message)
{
  BOOST_LOG_SEV(loggerImpl->lg, severity) << filename << ':' << line << " in function " << function_name << " : "
                                          << message;
}

void Logger::flush()
{
  ::loggerImpl->flush();
}
