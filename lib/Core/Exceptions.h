#ifndef LATTICEGEN_EXCEPTIONS_H
#define LATTICEGEN_EXCEPTIONS_H

#include "Logger.h"

#define MULTILINE_MACRO_BEGIN do { \

#define MULTILINE_MACRO_END } while (false)

#define THROW_INVALID_ARGUMENT(msg__) \
MULTILINE_MACRO_BEGIN \
  Core::Logger::log(Core::Logger::Severity::Error, __FILE__, __LINE__, __func__, std::string("std::invalid_argument : ") + msg__); \
  throw std::invalid_argument(msg__); \
MULTILINE_MACRO_END

#define THROW_LOGIC_ERROR(msg__) \
MULTILINE_MACRO_BEGIN \
  Core::Logger::log(Core::Logger::Severity::Error, __FILE__, __LINE__, __func__, std::string("std::logic_error : ") + msg__); \
  throw std::logic_error(msg__); \
MULTILINE_MACRO_END


#endif //LATTICEGEN_EXCEPTIONS_H
