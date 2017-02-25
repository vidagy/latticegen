# Core

file(GLOB all_core_test_src "Core/*.cpp")
add_executable(Test_Core ${all_core_test_src})
target_link_libraries(Test_Core gmock_main pthread)
install(TARGETS Test_Core DESTINATION tests/)
add_test(NAME Test_Core COMMAND Test_Core --gtest_output=xml:testresults_core.xml)

# Geometry

file(GLOB all_geometry_test_src "Geometry/*.cpp")
add_executable(Test_Geometry ${all_geometry_test_src})
target_link_libraries(Test_Geometry Geometry gmock_main pthread)
install(TARGETS Test_Geometry DESTINATION tests/)
add_test(NAME Test_Geometry COMMAND Test_Geometry --gtest_output=xml:testresults_geometry.xml)

# Math

file(GLOB all_math_test_src "Math/*.cpp")
add_executable(Test_Math ${all_math_test_src})
target_link_libraries(Test_Math Math gmock_main pthread)
install(TARGETS Test_Math DESTINATION tests/)
add_test(NAME Test_Math COMMAND Test_Math --gtest_output=xml:testresults_math.xml)
