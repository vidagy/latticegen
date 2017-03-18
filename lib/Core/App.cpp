#include "App.h"
#include "Logger.h"

#include <boost/program_options.hpp>
#include <iostream>

using namespace Core;

std::once_flag App::onceFlag;
std::map<std::string, std::string> App::args;
std::map<std::string, std::string> App::environment_variables;

void App::initialize(int argc, char const *const *argv)
{
  std::call_once(onceFlag, initialize_impl, argc, argv);
}

void App::initialize_impl(int argc, char const *const *argv)
{
  using namespace boost::program_options;
  try {
    options_description desc{"Options"};
    desc.add_options()
      ("help,h", "Help")
      ("log_level,ll", value<std::string>()->default_value("info"), "Log level")
      ("log_target,lt", value<std::string>()->default_value("cout"), "Log target: cout or filename");

    variables_map vm;
    store(
      command_line_parser(argc, argv)
        .options(desc)
        .style(command_line_style::unix_style | command_line_style::case_insensitive)
        .run(), vm);

    notify(vm);

    if (vm.count("help")) {
      std::cout << desc << '\n';
      exit(0);
    }
    if (vm.count("log_level")) {
      args["LOG_LEVEL"] = vm["log_level"].as<std::string>();
    }
    if (vm.count("log_target")) {
      args["LOG_TARGET"] = vm["log_target"].as<std::string>();
    }
  }
  catch (const error &ex) {
    std::cerr << ex.what() << std::endl;
  }

  for (char **it = environ; *it != nullptr; ++it) {
    std::string entry(*it);
    auto sep = entry.find('=');
    if (sep != std::string::npos) {
      auto key = entry.substr(0, sep);
      std::transform(key.begin(), key.end(), key.begin(), ::toupper);
      environment_variables[key] = entry.substr(sep + 1);
    }
  }

  Logger::initialize();
}
