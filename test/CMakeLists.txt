# TestUtils

file(GLOB all_test_utils_hdr "TestUtils/*.h")
file(GLOB all_test_utils_src "TestUtils/*.cpp")
add_library(TestUtils STATIC ${all_test_utils_src} ${all_test_utils_hdr})
target_link_libraries(TestUtils gmock pthread)
target_include_directories(TestUtils PUBLIC .)

# Core

file(GLOB all_core_test_src "Core/*.cpp")
add_executable(Test_Core ${all_core_test_src} ${test_utils})
target_link_libraries(Test_Core Core TestUtils)
install(TARGETS Test_Core DESTINATION tests/)
add_test(NAME Test_Core COMMAND Test_Core --gtest_output=xml:testresults_core.xml)

# Geometry

file(GLOB all_geometry_test_src "Geometry/*.cpp")
add_executable(Test_Geometry ${all_geometry_test_src} ${test_utils})
target_link_libraries(Test_Geometry Geometry TestUtils)
install(TARGETS Test_Geometry DESTINATION tests/)
add_test(NAME Test_Geometry COMMAND Test_Geometry --gtest_output=xml:testresults_geometry.xml)

# Math

file(GLOB all_math_test_src "Math/*.cpp")
add_executable(Test_Math ${all_math_test_src})
target_link_libraries(Test_Math Math TestUtils)
install(TARGETS Test_Math DESTINATION tests/)
add_test(NAME Test_Math COMMAND Test_Math --gtest_output=xml:testresults_math.xml)

# Physics

file(GLOB all_physics_test_src
        "Physics/NonRelativistic/CoreElectrons/*.cpp"
        "Physics/NonRelativistic/*.cpp")
add_executable(Test_Physics ${all_physics_test_src})
target_link_libraries(Test_Physics Physics TestUtils)
install(TARGETS Test_Physics DESTINATION tests/)
add_test(NAME Test_Physics COMMAND Test_Physics --gtest_output=xml:testresults_physics.xml)
