# Core

file(GLOB all_core_hdr "Core/*.h")
add_library(Core INTERFACE)
target_sources(Core INTERFACE ${all_core_hdr})
target_include_directories(Core INTERFACE .)

# Geometry

file(GLOB all_geometry_src "Geometry/*.cpp")
file(GLOB all_geometry_hdr "Geometry/*.h")
add_library(Geometry SHARED ${all_geometry_src} ${all_geometry_hdr})
# target_link_libraries(Geometry )
target_include_directories(Geometry PUBLIC .)

# Math

file(GLOB all_math_src "Math/*.cpp")
file(GLOB all_math_hdr "Math/*.h")
add_library(Math SHARED ${all_math_src} ${all_math_hdr})
# target_link_libraries(Math Core)
target_include_directories(Math PUBLIC .)

# Physics

file(GLOB all_physics_src "Physics/CoreElectrons/*.cpp")
file(GLOB all_physics_hdr "Physics/CoreElectrons/*.h")
add_library(Physics SHARED ${all_physics_src} ${all_physics_hdr})
target_link_libraries(Physics Math)
target_include_directories(Physics PUBLIC .)

