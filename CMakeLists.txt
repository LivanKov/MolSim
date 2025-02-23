cmake_minimum_required(VERSION 3.10)

# define project name, version
project(PSEMolDyn_GroupH VERSION 0.0.1)

# let ccmake and cmake-gui offer the default build type options
set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug;Release;RelWithDebInfo;MinSizeRel")

# set Release as the default build type if it is not yet set.
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE "Release" CACHE STRING
            "Choose the type of build, options are: Debug Release RelWithDebInfo MinSizeRel." FORCE)
endif(NOT CMAKE_BUILD_TYPE)


# Include the FetchContent module
include(FetchContent)

# Fetch the latest version of spdlog
FetchContent_Declare(
        spdlog
        GIT_REPOSITORY https://github.com/gabime/spdlog.git
        GIT_TAG v1.x
)

FetchContent_MakeAvailable(spdlog)

# Find OpenMP
find_package(OpenMP REQUIRED)

# collect all cpp files
file(GLOB_RECURSE MY_SRC
    "${CMAKE_CURRENT_SOURCE_DIR}/src/*.cpp"
    # header don't need to be included but this might be necessary for some IDEs
    "${CMAKE_CURRENT_SOURCE_DIR}/src/*.h"
)

# Add the generated .cxx files from CodeSynthesis XSD
set(GENERATED_SRC
    ${CMAKE_CURRENT_SOURCE_DIR}/src/io/input/MolSim.cxx
)   

# create make target
add_executable(MolSim ${MY_SRC} ${GENERATED_SRC})

# set cxx standard. You may raise this if you want.
target_compile_features(MolSim
        PRIVATE
            cxx_std_20
)

target_include_directories(MolSim
        PUBLIC
            ${CMAKE_CURRENT_SOURCE_DIR}/libs/libxsd
        PRIVATE
            ${CMAKE_CURRENT_SOURCE_DIR}/src
)

target_link_libraries(MolSim
        # stuff that is used in headers and source files
        PUBLIC
        xerces-c
        OpenMP::OpenMP_CXX # Link OpenMP library
        PRIVATE
        spdlog::spdlog
)

option(OPEN_MP "Build with OpenMP" OFF)
if(OPEN_MP)
    add_compile_definitions(PARALLEL)
endif()

# activate all compiler warnings. Clean up your code :P
# depending on the compiler different flags are used
target_compile_options(MolSim
        PRIVATE
            $<$<CXX_COMPILER_ID:GNU>:-Wall>
            # disable some spam warnings for icpc...
            $<$<CXX_COMPILER_ID:Intel>:-w3 -wd383,981,1418,1572,2259>
)

# Add OpenMP-specific compiler options
if(OpenMP_CXX_FOUND)
    target_compile_options(MolSim PRIVATE ${OpenMP_CXX_FLAGS})
endif()

add_subdirectory(tests ${CMAKE_BINARY_DIR}/tests)

enable_testing()
add_test(NAME BasicTest COMMAND tests/basic_test)

list(APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake/modules)

include(doxygen)