cmake_minimum_required(VERSION 3.10)
project(GRAM LANGUAGES CXX)

# Set C++ Standard to C++11 (you can change this if needed)
set(CMAKE_CXX_STANDARD 20)
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
# Find the IPCL library
find_package(IPCL 2.0.0 REQUIRED)
# Find the GMP library

set(COMMON_INCLUDE_DIRS 
. crypto gadgets storage tests types utils scripts emp ${OPENSSL_INCLUDE_DIR}
)

# Try to locate FLINT (Fast Library for Number Theory)
# You can point CMake to a custom FLINT installation by exporting FLINT_ROOT
#   export FLINT_ROOT=/path/to/flint/prefix (containing include/ and lib/)
set(_FLINT_HINTS)
if(DEFINED ENV{FLINT_ROOT})
    list(APPEND _FLINT_HINTS $ENV{FLINT_ROOT})
endif()

find_path(FLINT_INCLUDE_DIR
    NAMES flint/flint.h
    HINTS ${_FLINT_HINTS}
    PATH_SUFFIXES include
    PATHS
        /usr/local/include
        /usr/include
        /opt/local/include
)

find_library(FLINT_LIBRARY
    NAMES flint
    HINTS ${_FLINT_HINTS}
    PATH_SUFFIXES lib lib64
    PATHS
        /usr/local/lib
        /usr/local/lib64
        /usr/lib
        /usr/lib64
        /usr/lib/x86_64-linux-gnu
        /opt/local/lib
)

if(FLINT_LIBRARY)
    get_filename_component(FLINT_LIBRARY_DIR ${FLINT_LIBRARY} DIRECTORY)
    message(STATUS "Found FLINT library: ${FLINT_LIBRARY}")
    if(FLINT_INCLUDE_DIR)
        message(STATUS "Found FLINT include dir: ${FLINT_INCLUDE_DIR}")
        list(APPEND COMMON_INCLUDE_DIRS ${FLINT_INCLUDE_DIR})
    endif()
else()
    message(WARNING "FLINT library not found in standard locations. Will try to link with '-lflint'. Set FLINT_ROOT to help CMake find it.")
endif()

if(FLINT_LIBRARY)
    set(FLINT_LIB_TO_LINK ${FLINT_LIBRARY})
else()
    set(FLINT_LIB_TO_LINK flint)
endif()

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
target_link_libraries(test_basic gtest gtest_main OpenSSL::SSL OpenSSL::Crypto OpenMP::OpenMP_CXX gmp gmpxx ${FLINT_LIB_TO_LINK} IPCL::ipcl)
target_compile_options(test_basic PRIVATE -maes)
target_compile_definitions(test_basic PRIVATE TESTING=1 IPCL_AVAILABLE=1)
target_include_directories(test_basic PRIVATE /opt/intel/ipcl/include/ipcl/ippcrypto /opt/intel/ipcl/include/ipcl/)

add_library(gramimpl ${COMMON_SOURCES})
target_link_libraries(gramimpl OpenSSL::SSL OpenSSL::Crypto OpenMP::OpenMP_CXX gmp gmpxx ${FLINT_LIB_TO_LINK})
target_compile_options(gramimpl PRIVATE -maes)
target_compile_definitions(gramimpl PRIVATE IPCL_AVAILABLE=1)
target_include_directories(gramimpl PRIVATE /opt/intel/ipcl/include/ipcl/ippcrypto /opt/intel/ipcl/include/ipcl/)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -march=native -Wall -Wextra -Wno-deprecated-declarations -DEMP_CHANNEL_PATH")

# Release build specific flags
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -DNDEBUG")

# Debug build specific flags
set(CMAKE_CXX_FLAGS_DEBUG "-g -O0")

# Add tests to CTest
add_test(NAME Test_Basic COMMAND test_basic)

# Embed RPATH so executables and shared libs can find FLINT at runtime
if(FLINT_LIBRARY)
    # Respect typical install paths and the FLINT lib dir
    set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)
    list(APPEND CMAKE_INSTALL_RPATH "${FLINT_LIBRARY_DIR}")
    list(REMOVE_DUPLICATES CMAKE_INSTALL_RPATH)

    # Also set build-tree rpath for the binary in build/ so running from build works
    set(_build_rpath "${FLINT_LIBRARY_DIR}")
    foreach(tgt test_basic gramimpl)
        if(TARGET ${tgt})
            set_target_properties(${tgt} PROPERTIES
                BUILD_RPATH "${_build_rpath}"
                INSTALL_RPATH "${CMAKE_INSTALL_RPATH}"
            )
        endif()
    endforeach()
endif()