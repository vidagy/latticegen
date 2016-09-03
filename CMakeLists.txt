cmake_minimum_required(VERSION 3.5 FATAL_ERROR)

# project name
project (LatticeGen CXX)

set(MULTITHREADED_BUILD "4")
set(CMAKE_MAKE_PROGRAM "${CMAKE_MAKE_PROGRAM} -j${MULTITHREADED_BUILD}")
message(STATUS "Added arguments to CMAKE_BUILD_TOOL: ${CMAKE_MAKE_PROGRAM}")

# Setting compipler
set (CMAKE_CXX_COMPILER g++)

# Setting compiler flags for debug and release versions
set(CMAKE_CXX_FLAGS         "-std=c++11 ${CMAKE_CXX_FLAGS}")
set(CMAKE_CXX_FLAGS_DEBUG   "-std=c++11 ${CMAKE_CXX_FLAGS_DEBUG}   -g -pg -O0 -Wall")
set(CMAKE_CXX_FLAGS_RELEASE "-std=c++11 ${CMAKE_CXX_FLAGS_RELEASE} -O2 ")

# install directory
set(CMAKE_INSTALL_PREFIX "${CMAKE_SOURCE_DIR}")

# subdirectories
add_subdirectory(lib/)
include_directories("lib")

# set up testing
enable_testing()

set(GTEST_LIBRARY_DEBUG "/usr/src/gtest/build/libgtest.a")
set(GTEST_MAIN_LIBRARY_DEBUG "/usr/src/gtest/build/libgtest_main.a")
include_directories("/usr/include/")

add_subdirectory(test)

##########################################################################
#set( CTEST_MEMORYCHECK_COMMAND "/usr/bin/valgrind")
#set( CTEST_MEMORYCHECK_COMMAND_OPTIONS "--tool=memcheck" ) # --xml=yes --xml-file=memtest.xml
#include(CTest)
