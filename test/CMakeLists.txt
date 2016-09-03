# Core
file(GLOB all_core_test_src "Core/*.cpp")
add_executable(Test_Core ${all_core_test_src})
target_link_libraries(Test_Core ${GTEST_LIBRARY_DEBUG} ${GTEST_MAIN_LIBRARY_DEBUG} pthread LatticeGenCore)
install(TARGETS Test_Core DESTINATION tests/)
add_test(NAME Test_Core COMMAND Test_Core)