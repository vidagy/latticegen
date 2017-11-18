cmake_minimum_required(VERSION 3.2 FATAL_ERROR)

# project name
project (LatticeGen CXX)

# DEPENDENCIES

# googletest
include(ExternalProject)
include(include/include.cmake)

# boost
set(BOOST_COMPONENTS program_options log log_setup system)
find_package(Boost 1.63.0 REQUIRED COMPONENTS ${BOOST_COMPONENTS})
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
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++17")
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

# Doxygen documentation
find_package(Doxygen)
option(BUILD_DOCUMENTATION "Create and install HTML based API documentation (requires Doxygen)" ${DOXYGEN_FOUND})
if (BUILD_DOCUMENTATION)
    if (NOT DOXYGEN_FOUND)
        message(FATAL_ERROR "Doxygen is needed to build the documentation.")
    endif ()

    set(doxyfile_in ${CMAKE_CURRENT_SOURCE_DIR}/Doxyfile)
    set(doxyfile ${CMAKE_CURRENT_BINARY_DIR}/doxyfile)

    configure_file(${doxyfile_in} ${doxyfile} @ONLY)
    add_custom_target(
            doc
            COMMAND ${DOXYGEN_EXECUTABLE} ${doxyfile_in}
            WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
            COMMENT "Generating API documentation with Doxygen"
            VERBATIM)
endif ()
