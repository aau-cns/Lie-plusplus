cmake_minimum_required(VERSION 3.13)

project(
    lieplusplus
    VERSION 1.0
    DESCRIPTION "Lie++ A header-only Eigen-based C++ library for Lie group operations"
    LANGUAGES CXX
)

## Options
if (NOT BUILD_TESTS)
    set(BUILD_TESTS ON)
endif()

## Include and set up external libraries
include(FetchContent)

# Googletest
if(${BUILD_TESTS})
    FetchContent_Declare(
        googletest
        GIT_REPOSITORY  https://github.com/google/googletest.git
        GIT_TAG         release-1.12.1
        GIT_SHALLOW     TRUE
        GIT_PROGRESS    TRUE
    )
    set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
    list(APPEND external googletest)
    list(APPEND include_dirs ${googletest_INCLUDE_DIR})
    list(APPEND libs ${googletest_LIBRARIES})
endif()

# Eigen
FetchContent_Declare(
    Eigen
    GIT_REPOSITORY  https://gitlab.com/libeigen/eigen.git
    GIT_TAG         3.4.0
    GIT_SHALLOW     TRUE
    GIT_PROGRESS    TRUE
)
set(CMAKE_POLICY_DEFAULT_CMP0077 NEW)
set(EIGEN_BUILD_DOC OFF)
set(EIGEN_BUILD_PKGCONFIG OFF)
list(APPEND external Eigen)
list(APPEND include_dirs ${EIGEN3_INCLUDE_DIR})
list(APPEND libs Eigen3::Eigen ${EIGEN3_LIBRARIES})

FetchContent_MakeAvailable(${external})

message(STATUS "Eigen:    " ${EIGEN3_VERSION})

## Try to compile with C++ 17 (make sure GCC>=7, clang>=5, MSVC>=19.12)
include(CheckCXXCompilerFlag)
CHECK_CXX_COMPILER_FLAG("-std=c++17" COMPILER_SUPPORTS_CXX17)
if(COMPILER_SUPPORTS_CXX17)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++17")
    set(CMAKE_CXX_STANDARD 17)
else()
    message(STATUS "The compiler ${CMAKE_CXX_COMPILER} has no c++17 support. Please use a different C++ compiler.")
endif()
set(CMAKE_CXX_STANDARD_REQUIRED ON)
message(STATUS "Set compiler ${CMAKE_CXX_COMPILER}.")

## Define build type, Release by default
if (NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE "Release")
endif()
message(STATUS "CMAKE_BUILD_TYPE ${CMAKE_BUILD_TYPE}.")

# Set CMAKE compile flags
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS} -DNDEBUG -O3 -fsee -fomit-frame-pointer -fno-signed-zeros -fno-math-errno -funroll-loops -fno-unsafe-math-optimizations")
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS} -O0 -g3 -Wall -Wextra -Werror -Wuninitialized -Wmaybe-uninitialized -pedantic -fno-omit-frame-pointer -fsanitize=address -fsanitize-recover=address")
# set(CMAKE_CXX_FLAGS_DEBUG "-O0 -g3 -Wall -Wextra -Werror -Wuninitialized -Wmaybe-uninitialized -pedantic -fno-omit-frame-pointer -fsanitize=undefined")
message(STATUS "CMAKE_CXX_FLAGS ${CMAKE_CXX_FLAGS}.")

## Define includes
list(APPEND include_dirs ${CMAKE_CURRENT_SOURCE_DIR}/include)
include_directories(${include_dirs})

## Declare C++ tests
if(${BUILD_TESTS})
    add_executable(lieplusplus_tests tests/tests.cpp)
    target_include_directories(lieplusplus_tests PRIVATE ${include_dirs})
    target_link_libraries(lieplusplus_tests GTest::gtest_main)
    include(GoogleTest)
    gtest_discover_tests(lieplusplus_tests)
endif()
