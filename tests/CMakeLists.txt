# Google Test dependency
include(FetchContent)
FetchContent_Declare(
  googletest
  URL https://github.com/google/googletest/archive/b514bdc898e2951020cbdca1304b75f5950d1f59.zip
)
FetchContent_MakeAvailable(googletest)
enable_testing()

# Find OpenMP
find_package(OpenMP REQUIRED)

# add_executable(basic_test  boundary_conditions_test.cpp ../src/simulator/particle/container/DirectSumContainer.cpp ../src/simulator/particle/Particle.cpp ../src/utils/logger/Logger.cpp
#         ../src/simulator/particle/ParticleGenerator.cpp ../src/simulator/particle/container/LinkedCellContainer.cpp ../src/simulator/particle/container/DirectSumContainer.cpp ../src/simulator/Thermostat.cpp 
#         ../src/simulator/calculations/Calculation.h ../src/simulator/calculations/BoundaryConditions.cpp ../src/simulator/calculations/Position.cpp ../src/io/input/cli/SimParams.cpp ../src/simulator/calculations/Force.h ../src/simulator/calculations/Force.cpp)

add_executable(basic_test basic_test.cpp generator_test.cpp linked_cell.cpp boundary_conditions_test.cpp thermostat_test.cpp periodic_conditions_test.cpp ../src/simulator/particle/container/DirectSumContainer.cpp ../src/simulator/particle/Particle.cpp ../src/utils/logger/Logger.cpp
        ../src/simulator/particle/ParticleGenerator.cpp ../src/simulator/particle/container/LinkedCellContainer.cpp ../src/simulator/particle/container/DirectSumContainer.cpp ../src/simulator/Thermostat.cpp 
        ../src/simulator/calculations/Calculation.h ../src/simulator/calculations/BoundaryConditions.cpp ../src/simulator/calculations/Position.cpp ../src/io/input/cli/SimParams.cpp ../src/simulator/calculations/Force.h ../src/simulator/calculations/Force.cpp)   

target_include_directories(basic_test PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/../src)

target_link_libraries(
  basic_test
  GTest::gtest_main
  spdlog::spdlog
  OpenMP::OpenMP_CXX
)
target_compile_features(basic_test
        PRIVATE
            cxx_std_20
)

# Add OpenMP-specific compiler options
if(OpenMP_CXX_FOUND)
    target_compile_options(MolSim PRIVATE ${OpenMP_CXX_FLAGS})
endif()

include(GoogleTest)
gtest_discover_tests(basic_test)
