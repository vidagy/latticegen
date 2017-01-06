# Core

file(GLOB all_core_hdr "Core/*.h")
add_library(Core INTERFACE)
target_sources(Core INTERFACE ${all_core_hdr})
target_include_directories(Core INTERFACE .)

# Geometry

file(GLOB all_geometry_src "Geometry/*.cpp")
file(GLOB all_geometry_hdr "Geometry/*.h")
add_library(Geometry SHARED ${all_geometry_src} ${all_geometry_hdr})
target_link_libraries(Geometry Core)

