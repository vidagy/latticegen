cmake_minimum_required(VERSION 3.5 FATAL_ERROR)

# project name
project (LatticeGen CXX)

include(ExternalProject)
include(include/include.cmake)

find_package(Boost REQUIRED)
add_definitions(-DBOOST_LOG_DYN_LINK)

# Setting compipler
set (CMAKE_CXX_COMPILER g++)

# Setting compiler flags for debug and release versions
set(CMAKE_CXX_FLAGS         "-std=c++14 ${CMAKE_CXX_FLAGS}")
set(CMAKE_CXX_FLAGS_DEBUG   "${CMAKE_CXX_FLAGS_DEBUG} -fprofile-arcs -ftest-coverage -g -pg -O0 -Wall -Werror")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O2 ")

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
