cmake_minimum_required(VERSION 3.2 FATAL_ERROR)

# project name
project (LatticeGen CXX)

# DEPENDENCIES

# googletest
include(ExternalProject)
include(include/include.cmake)

# boost
set(BOOST_COMPONENTS program_options log log_setup system)
find_package(Boost 1.58.0 REQUIRED COMPONENTS ${BOOST_COMPONENTS})
set(Boost_DEBUG OFF)
set(Boost_USE_STATIC_LIBS OFF)
set(Boost_USE_MULTITHREADED ON)
set(Boost_USE_STATIC_RUNTIME OFF)
set(BOOST_ALL_DYN_LINK ON)
add_definitions(-DBOOST_LOG_DYN_LINK)

# eigen
find_package(Eigen3 REQUIRED NO_MODULE)
# TODO generalize this somehow, it looks extremely lame...
include_directories(/usr/include/eigen3/)

# Setting compiler flags for debug and release versions
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS_DEBUG} -std=c++14")
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -fprofile-arcs -ftest-coverage -g -O0 -Wall -Werror")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O2 -g -Wall -Werror")

# install directory
set(CMAKE_INSTALL_PREFIX "${CMAKE_SOURCE_DIR}")

# subdirectories
add_subdirectory(lib/)
include_directories("lib")

# set up testing
enable_testing()

add_subdirectory(test)
