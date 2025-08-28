cmake_minimum_required(VERSION 3.20)
project(noLZSS LANGUAGES CXX)

# Prefer modern FindPython/pybind11 behavior and avoid warnings about
# deprecated FindPythonInterp/FindPythonLibs usage in older modules.
cmake_policy(SET CMP0148 NEW)
set(PYBIND11_FINDPYTHON ON)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

# This CMake configuration only supports vendoring the SDSL v3 (xxsds fork).
# It does not attempt to detect or use a system-installed sdsl.

include(FetchContent)

# Allow overriding the sdsl tag for advanced users, but default to a stable release.
set(SDSL_LITE_GIT_TAG "v3.0.3" CACHE STRING "Tag/branch to fetch for sdsl-lite (xxsds fork)")
message(STATUS "Vendoring sdsl-lite (xxsds) tag ${SDSL_LITE_GIT_TAG}")
FetchContent_Declare(sdsl-lite
  GIT_REPOSITORY https://github.com/xxsds/sdsl-lite.git
  GIT_TAG        ${SDSL_LITE_GIT_TAG}
)
# Use modern FetchContent workflow to make the declared content available.
FetchContent_MakeAvailable(sdsl-lite)
# Do NOT execute the sdsl-lite subproject CMake. Treat SDSL v3 (xxsds fork)
# as header-only for downstream users by exposing an INTERFACE target that
# points at the fetched headers. This avoids running third-party CMake during
# pip/PEP builds and avoids CMake policy/version incompatibilities.
if (NOT DEFINED sdsl-lite_SOURCE_DIR)
  message(FATAL_ERROR "sdsl-lite source directory not available after FetchContent")
endif()

# Create an INTERFACE target that exposes the sdsl include dirs from the
# fetched source tree. Also include the bundled libdivsufsort headers if
# present; this matches how consumers expect to include SDSL v3 headers.
add_library(sdsl_vendor INTERFACE)
target_include_directories(
  sdsl_vendor
  INTERFACE
    $<BUILD_INTERFACE:${sdsl-lite_SOURCE_DIR}/include>
    $<BUILD_INTERFACE:${sdsl-lite_SOURCE_DIR}/external/libdivsufsort/include>
)
set(SDSL_TARGET sdsl_vendor)
message(STATUS "Created sdsl_vendor INTERFACE target exposing sdsl-lite headers from ${sdsl-lite_SOURCE_DIR}/include")

# Vendoring pybind11 if not available
find_package(Python COMPONENTS Interpreter Development QUIET)
# Vendoring pybind11 if not available; ensure Python is discovered first so
# pybind11 uses the modern FindPython implementation.
if(NOT Python_FOUND)
  # If Python wasn't found, require it explicitly (most environments will
  # have Python available via pybind11's own discovery but setting this
  # ensures the policy CMP0148 doesn't trigger warnings).
  find_package(Python COMPONENTS Interpreter Development REQUIRED)
endif()

find_package(pybind11 CONFIG QUIET)
if (NOT pybind11_FOUND)
  message(STATUS "Vendoring pybind11 via FetchContent")
  FetchContent_Declare(pybind11
    GIT_REPOSITORY https://github.com/pybind/pybind11.git
    GIT_TAG        v2.10.4
  )
  FetchContent_MakeAvailable(pybind11)
endif()

add_library(noLZSS_core STATIC src/noLZSS/factorizer.cpp)
target_link_libraries(noLZSS_core PUBLIC ${SDSL_TARGET})

pybind11_add_module(_noLZSS src/noLZSS/bindings.cpp)
add_dependencies(_noLZSS noLZSS_core)
target_link_libraries(_noLZSS PRIVATE noLZSS_core)

install(TARGETS _noLZSS
  LIBRARY DESTINATION noLZSS
  RUNTIME DESTINATION noLZSS
  ARCHIVE DESTINATION noLZSS
)
