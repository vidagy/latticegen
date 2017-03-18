# Core

file(GLOB all_core_test_src "main.cpp" "Core/*.cpp")
add_executable(Test_Core ${all_core_test_src})
target_link_libraries(Test_Core Core gmock pthread)
install(TARGETS Test_Core DESTINATION tests/)
add_test(NAME Test_Core COMMAND Test_Core --gtest_output=xml:testresults_core.xml)

# Geometry

file(GLOB all_geometry_test_src "main.cpp" "Geometry/*.cpp")
add_executable(Test_Geometry ${all_geometry_test_src})
target_link_libraries(Test_Geometry Geometry gmock pthread)
install(TARGETS Test_Geometry DESTINATION tests/)
add_test(NAME Test_Geometry COMMAND Test_Geometry --gtest_output=xml:testresults_geometry.xml)

# Math

file(GLOB all_math_test_src "main.cpp" "Math/*.cpp")
add_executable(Test_Math ${all_math_test_src})
target_link_libraries(Test_Math Math gmock pthread)
install(TARGETS Test_Math DESTINATION tests/)
add_test(NAME Test_Math COMMAND Test_Math --gtest_output=xml:testresults_math.xml)

# Physics

file(GLOB all_physics_test_src "main.cpp" "Physics/CoreElectrons/*.cpp")
add_executable(Test_Physics ${all_physics_test_src})
target_link_libraries(Test_Physics Physics gmock pthread)
install(TARGETS Test_Physics DESTINATION tests/)
add_test(NAME Test_Physics COMMAND Test_Physics --gtest_output=xml:testresults_physics.xml)
