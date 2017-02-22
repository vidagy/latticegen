# Core

file(GLOB all_core_src "Core/*.cpp")
file(GLOB all_core_hdr "Core/*.h")
add_library(Core SHARED ${all_core_src} ${all_core_hdr})
target_link_libraries(Core)
target_include_directories(Core PUBLIC .)

# Geometry

file(GLOB all_geometry_src "Geometry/*.cpp")
file(GLOB all_geometry_hdr "Geometry/*.h")
add_library(Geometry SHARED ${all_geometry_src} ${all_geometry_hdr})
target_link_libraries(Geometry Core)
target_include_directories(Geometry PUBLIC .)

