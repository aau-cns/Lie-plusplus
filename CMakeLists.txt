cmake_minimum_required(VERSION 3.13)

project(
    LiePlusPlus
    VERSION 1.0
    DESCRIPTION "Lie++ A header-only Eigen-based C++ library for Lie group operations"
    LANGUAGES CXX
)

# Set build type, Release by default
if (NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE "Release")
endif()

# Set options
option(LIEPLUSPLUS_TESTS "Build Lie++ tests" OFF)
option(ENABLE_ADDRESS_SANITIZER "Enable address sanitizer" OFF)
option(ENABLE_UNDEFINED_SANITIZER "Enable undefined behavior sanitizer" ON)

## Include and set up external libraries
include(FetchContent)

# Googletest
if(LIEPLUSPLUS_TESTS)
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
set(EIGEN_TEST_CXX11 ON)
set(EIGEN_HAS_CXX11_MATH ON)
list(APPEND external Eigen)
list(APPEND libs Eigen3::Eigen)

FetchContent_MakeAvailable(${external})

include(CheckCXXCompilerFlag)

# Set compiler flags
set(RELEASE_FLAGS "-DNDEBUG" "-O3" "-fsee" "-fomit-frame-pointer" "-fno-signed-zeros" "-fno-math-errno" "-funroll-loops" "-fno-unsafe-math-optimizations" "-flto" "-march=native")
set(DEBUG_FLAGS "-O0" "-g3" "-Wall" "-Wextra" "-Werror" "-Wuninitialized" "-Wmaybe-uninitialized" "-pedantic" "-fno-omit-frame-pointer")

CHECK_CXX_COMPILER_FLAG("-std=c++17" COMPILER_SUPPORTS_CXX17)
if(COMPILER_SUPPORTS_CXX17)
    set(CMAKE_CXX_STANDARD 17)
    set(CMAKE_CXX_STANDARD_REQUIRED ON)
else()
    message(STATUS "The compiler ${CMAKE_CXX_COMPILER} has no c++17 support. Please use a different C++ compiler.")
endif()
if (CMAKE_BUILD_TYPE STREQUAL "Debug")
    foreach(FLAG ${DEBUG_FLAGS})
        string(REPLACE "-" "" FLAG_NAME ${FLAG})
        string(TOUPPER ${FLAG_NAME} FLAG_VAR_NAME)
        CHECK_CXX_COMPILER_FLAG("${FLAG}" COMPILER_SUPPORTS_${FLAG_VAR_NAME})
        if (COMPILER_SUPPORTS_${FLAG_VAR_NAME})
            set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} ${FLAG}")
        else()
            message("${FLAG} was requested but is not supported.")
        endif()
    endforeach()
    if (ENABLE_ADDRESS_SANITIZER)
        CHECK_CXX_COMPILER_FLAG("-fsanitize=address" COMPILER_SUPPORTS_ADDRESS_SANITIZER)
        if (COMPILER_SUPPORTS_ADDRESS_SANITIZER)
            set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -fsanitize=address -fsanitize-recover=address")
        else()
            message("-fsanitize=address -fsanitize-recover=address were requested but is not supported.")
        endif()
    elseif (ENABLE_UNDEFINED_SANITIZER)
        CHECK_CXX_COMPILER_FLAG("-fsanitize=undefined" COMPILER_SUPPORTS_UNDEFINED_SANITIZER)
        if (COMPILER_SUPPORTS_UNDEFINED_SANITIZER)
            set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -fsanitize=undefined")
        else()
            message("-fsanitize=undefined was requested but is not supported.")
        endif()
    endif()
else()
    foreach(FLAG ${RELEASE_FLAGS})
        string(REPLACE "-" "" FLAG_NAME ${FLAG})
        string(TOUPPER ${FLAG_NAME} FLAG_VAR_NAME)
        CHECK_CXX_COMPILER_FLAG("${FLAG}" COMPILER_SUPPORTS_${FLAG_VAR_NAME})
        if (COMPILER_SUPPORTS_${FLAG_VAR_NAME})
            set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} ${FLAG}")
        else()
            message("${FLAG} was requested but is not supported.")
        endif()
    endforeach()
endif()

message(STATUS "COMPILER: " ${CMAKE_CXX_COMPILER})
message(STATUS "CMAKE_BUILD_TYPE: " ${CMAKE_BUILD_TYPE})

## Define includes
list(APPEND include_dirs ${CMAKE_CURRENT_SOURCE_DIR}/include)
include_directories(${include_dirs})

# Add the library
add_library(LiePlusPlus INTERFACE)
target_include_directories(LiePlusPlus INTERFACE ${include_dirs})

if (LIEPLUSPLUS_TESTS)
    message(STATUS "Building Lie++ tests")
    enable_testing()
    add_executable(lieplusplus_tests tests/tests.cpp)
    target_link_libraries(lieplusplus_tests gtest_main LiePlusPlus ${libs})
    include(GoogleTest)
    gtest_discover_tests(lieplusplus_tests)
endif()