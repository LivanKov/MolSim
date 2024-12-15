# Google Test dependency
include(FetchContent)
FetchContent_Declare(
  googletest
  URL https://github.com/google/googletest/archive/b514bdc898e2951020cbdca1304b75f5950d1f59.zip
)
FetchContent_MakeAvailable(googletest)
enable_testing()

add_executable(basic_test basic_test.cpp generator_test.cpp linked_cell.cpp boundary_conditions_test.cpp ../src/simulator/particle/container/DirectSumContainer.cpp ../src/simulator/particle/Particle.cpp ../src/utils/logger/Logger.cpp
        ../src/simulator/particle/ParticleGenerator.cpp ../src/simulator/particle/container/LinkedCellContainer.cpp ../src/simulator/particle/container/DirectSumContainer.cpp ../src/simulator/Thermostat.cpp
        )

target_include_directories(basic_test PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/../src)

target_link_libraries(
  basic_test
  GTest::gtest_main
  spdlog::spdlog
)
target_compile_features(basic_test
        PRIVATE
            cxx_std_20
)

include(GoogleTest)
gtest_discover_tests(basic_test)
