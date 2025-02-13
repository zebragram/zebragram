cmake_minimum_required(VERSION 3.10)
project(GRAM LANGUAGES CXX)

# Set C++ Standard to C++11 (you can change this if needed)
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# option(BUILD_AS_SUBMODULE "Set to ON if building as a submodule of emp-tool" OFF)

# # Detect if included as a submodule
# if (NOT BUILD_AS_SUBMODULE AND NOT "${CMAKE_CURRENT_SOURCE_DIR}" STREQUAL "${CMAKE_SOURCE_DIR}")
#     set(BUILD_AS_SUBMODULE ON)
# endif()
# find_program(EMP_TOOL_PATH emp-tool)
find_package(emp-tool CONFIG)

# Find OpenSSL package
find_package(OpenSSL REQUIRED)
find_package(OpenMP REQUIRED)
# Find the GMP library

set(COMMON_INCLUDE_DIRS 
. crypto gadgets storage tests types utils scripts emp ${OPENSSL_INCLUDE_DIR}
)

if(EMP-TOOL_FOUND)
    message(STATUS "emp-tool is installed and found at: ${EMP-TOOL_INCLUDE_DIR}")
    list(APPEND COMMON_INCLUDE_DIRS ${EMP-TOOL_INCLUDE_DIR}/emp-tool/io)
else()
    message(STATUS "emp-tool not found: use replicated headers")
    list(APPEND COMMON_INCLUDE_DIRS emp-headers)
endif()

# Include the OpenSSL headers
include_directories(${COMMON_INCLUDE_DIRS})

# Enable testing
enable_testing()

# Download and configure GoogleTest if it's not installed
# This assumes the GTest source is in the 'external/googletest' folder
add_subdirectory(googletest)

file(GLOB TEST_SOURCES tests/*.cpp)
file(GLOB GADGET_SOURCES gadgets/*.cpp)
file(GLOB TYPE_SOURCES types/*.cpp)
file(GLOB UTIL_SOURCES utils/*.cpp)
file(GLOB STORAGE_SOURCES storage/*.cpp)
file(GLOB CRYPTO_SOURCES crypto/*.cpp)
file(GLOB EMP_SOURCE emp/*.cpp)
# combine all sources
set(COMMON_SOURCES ${GADGET_SOURCES} ${TYPE_SOURCES} ${UTIL_SOURCES} ${STORAGE_SOURCES} ${CRYPTO_SOURCES} ${EMP_SOURCE})
# Link GoogleTest to the test executable
add_executable(test_basic ${COMMON_SOURCES} ${TEST_SOURCES})
target_link_libraries(test_basic gtest gtest_main OpenSSL::SSL OpenSSL::Crypto OpenMP::OpenMP_CXX)
target_compile_options(test_basic PRIVATE -maes)

add_library(gramimpl ${COMMON_SOURCES})
target_link_libraries(gramimpl OpenSSL::SSL OpenSSL::Crypto OpenMP::OpenMP_CXX)
target_compile_options(gramimpl PRIVATE -maes)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -march=native -Wall -Wextra -DEMP_CHANNEL_PATH")

# Release build specific flags
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -DNDEBUG")

# Debug build specific flags
set(CMAKE_CXX_FLAGS_DEBUG "-g -O0")

# Add tests to CTest
add_test(NAME Test_Basic COMMAND test_basic)

# Cppcheck configuration
find_program(CPPCHECK_EXECUTABLE cppcheck)

if(CPPCHECK_EXECUTABLE)
    add_custom_target(
        cppcheck
        COMMAND ${CPPCHECK_EXECUTABLE} --enable=all --inconclusive --std=c++11 ${COMMON_SOURCES}
        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
        COMMENT "Running Cppcheck"
        VERBATIM
    )
else()
    message(WARNING "Cppcheck not found. Please install it.")
endif()