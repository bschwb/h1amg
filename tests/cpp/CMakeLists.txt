add_executable(test_h1amg
  test_main.cpp
  test_h1_helpers.cpp
  test_h1amg_helpers.cpp)

set_property(TARGET test_h1amg PROPERTY CXX_STANDARD 17)
set_property(TARGET test_h1amg PROPERTY CXX_STANDARD_REQUIRED ON)

add_dependencies(test_h1amg h1amg)
target_include_directories(test_h1amg
  PRIVATE "${CMAKE_SOURCE_DIR}/src")
target_link_libraries(test_h1amg h1amg ngsolve)

add_test(NAME cpp_tests COMMAND test_h1amg)
