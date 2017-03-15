#include "App.h"

#include <boost/program_options.hpp>
#include <iostream>

std::once_flag Core::App::onceFlag;
std::map<std::string, std::string> Core::App::args;

void Core::App::Create(int argc, char const *const *argv)
{
  std::call_once(onceFlag, set_up_args, argc, argv);
}

void Core::App::set_up_args(int argc, char const *const *argv)
{
  using namespace boost::program_options;
  try {
    options_description desc{"Options"};
    desc.add_options()
      ("help,h", "Help")
      ("log_level", value<std::string>()->default_value("info"), "Log level");

    variables_map vm;
    store(parse_command_line(argc, argv, desc), vm);
    notify(vm);

    if (vm.count("help"))
      std::cout << desc << '\n';
    else if (vm.count("log_level")) {
      // std::cout << "log_level: " << vm["log_level"].as<std::string>() << '\n';
      args["log_level"] = vm["log_level"].as<std::string>();
    }

  }
  catch (const error &ex) {
    std::cerr << ex.what() << std::endl;
  }
}
