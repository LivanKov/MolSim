# Google Test dependency
include(FetchContent)
FetchContent_Declare(
  googletest
  URL https://github.com/google/googletest/archive/b514bdc898e2951020cbdca1304b75f5950d1f59.zip
)
FetchContent_MakeAvailable(googletest)
enable_testing()
add_executable(
  hello_test
  tests/hello_test.cc
)
target_link_libraries(
  hello_test
  GTest::gtest_main
)
include(GoogleTest)
gtest_discover_tests(hello_test)
