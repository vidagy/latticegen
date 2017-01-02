# Build gtest from existing sources
ExternalProject_Add(
    project_gtest
    DOWNLOAD_COMMAND "" # No download required
    SOURCE_DIR "${CMAKE_SOURCE_DIR}/include/googletest/googletest" 
    BINARY_DIR "${CMAKE_BINARY_DIR}/include/googletest/googletest" 
    INSTALL_COMMAND ""
    LOG_CONFIGURE ON
    LOG_BUILD ON
)

ExternalProject_Get_Property(project_gtest source_dir)
ExternalProject_Get_Property(project_gtest binary_dir)

include_directories(${source_dir}/include)

add_library(gtest STATIC IMPORTED GLOBAL)
set_property(TARGET gtest PROPERTY IMPORTED_LOCATION ${binary_dir}/libgtest.a)
add_dependencies(gtest project_gtest)

add_library(gtest_main STATIC IMPORTED GLOBAL)
set_property(TARGET gtest_main PROPERTY IMPORTED_LOCATION ${binary_dir}/libgtest_main.a)
add_dependencies(gtest_main project_gtest)

# Build gmock from existing sources
ExternalProject_Add(
    project_gmock
    DOWNLOAD_COMMAND "" 
    SOURCE_DIR "${CMAKE_SOURCE_DIR}/include/googletest/googlemock" 
    BINARY_DIR "${CMAKE_BINARY_DIR}/include/googletest/googlemock"
    INSTALL_COMMAND ""
    LOG_CONFIGURE ON
    LOG_BUILD ON
)

ExternalProject_Get_Property(project_gmock source_dir)
ExternalProject_Get_Property(project_gmock binary_dir)

include_directories(${source_dir}/include)

add_library(gmock STATIC IMPORTED)
set_property(TARGET gmock PROPERTY IMPORTED_LOCATION ${binary_dir}/libgmock.a)
add_dependencies(gmock project_gmock)

add_library(gmock_main STATIC IMPORTED)
set_property(TARGET gmock_main PROPERTY IMPORTED_LOCATION ${binary_dir}/libgmock_main.a)
add_dependencies(gmock_main project_gmock)

