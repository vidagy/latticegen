# Core
file(GLOB all_core_src "Core/*.cpp")
file(GLOB all_core_hdr "Core/*.h")

add_library(LatticeGenCore SHARED ${all_core_src} ${all_core_hdr})

