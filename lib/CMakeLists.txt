# Core

file(GLOB all_core_hdr "Core/*.h")
file(GLOB all_core_src "Core/*.cpp")
add_library(Core SHARED ${all_core_src} ${all_core_hdr})
target_link_libraries(Core boost_program_options boost_log boost_log_setup boost_system)
target_include_directories(Core PUBLIC .)

# Geometry

file(GLOB all_geometry_src "Geometry/*.cpp")
file(GLOB all_geometry_hdr "Geometry/*.h")
add_library(Geometry SHARED ${all_geometry_src} ${all_geometry_hdr})
target_link_libraries(Geometry Core)
target_include_directories(Geometry PUBLIC .)

# Math

file(GLOB all_math_src "Math/*.cpp")
file(GLOB all_math_hdr "Math/*.h")
add_library(Math SHARED ${all_math_src} ${all_math_hdr})
# TODO something like find_package(LAPACKE REQUIRED) should be used instead just linking here
target_link_libraries(Math Core ${BLAS_LIBRARIES} lapacke)
target_include_directories(Math PUBLIC .)

# Physics

file(GLOB all_physics_src "Physics/CoreElectrons/*.cpp")
file(GLOB all_physics_hdr "Physics/CoreElectrons/*.h")
add_library(Physics SHARED ${all_physics_src} ${all_physics_hdr})
target_link_libraries(Physics Core Math)
target_include_directories(Physics PUBLIC .)

