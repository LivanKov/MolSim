# Google Test dependency
include(FetchContent)
FetchContent_Declare(
  googletest
  URL https://github.com/google/googletest/archive/b514bdc898e2951020cbdca1304b75f5950d1f59.zip
)
FetchContent_MakeAvailable(googletest)
enable_testing()
add_executable(basic_test basic_test.cc ../src/ParticleContainer.cpp ../src/Particle.cpp)

target_include_directories(basic_test PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/../src)

target_link_libraries(
  basic_test
  GTest::gtest_main
)
target_compile_features(basic_test
        PRIVATE
            cxx_std_20
)

include(GoogleTest)
gtest_discover_tests(basic_test)
