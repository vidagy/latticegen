cmake_minimum_required(VERSION 3.2 FATAL_ERROR)

# project name
project (LatticeGen CXX)

include(ExternalProject)
include(include/include.cmake)

set(BOOST_COMPONENTS program_options log log_setup system)
find_package(Boost 1.58.0 REQUIRED COMPONENTS ${BOOST_COMPONENTS})
set(Boost_DEBUG OFF)
set(Boost_USE_STATIC_LIBS OFF)
set(Boost_USE_MULTITHREADED ON)
set(Boost_USE_STATIC_RUNTIME OFF)
set(BOOST_ALL_DYN_LINK ON)
add_definitions(-DBOOST_LOG_DYN_LINK)

# Setting compiler flags for debug and release versions
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS_DEBUG} -std=c++14")
set(CMAKE_CXX_FLAGS_DEBUG   "${CMAKE_CXX_FLAGS_DEBUG} -fprofile-arcs -ftest-coverage -g -pg -O0 -Wall -Werror")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O2 -g -pg -Wall -Werror")

# install directory
set(CMAKE_INSTALL_PREFIX "${CMAKE_SOURCE_DIR}")

# subdirectories
add_subdirectory(lib/)
include_directories("lib")

# set up testing
enable_testing()

add_subdirectory(test)

##########################################################################
#set( CTEST_MEMORYCHECK_COMMAND "/usr/bin/valgrind")
#set( CTEST_MEMORYCHECK_COMMAND_OPTIONS "--tool=memcheck" ) # --xml=yes --xml-file=memtest.xml
#include(CTest)
