cmake_minimum_required(VERSION 3.20)

# Get version from scikit-build-core variable
if(DEFINED SKBUILD_PROJECT_VERSION)
    set(PROJECT_VER ${SKBUILD_PROJECT_VERSION})
else()
    set(PROJECT_VER 0.0.0)  # Fallback for manual builds
endif()

project(noLZSS VERSION ${PROJECT_VER} LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

# This CMake configuration only supports vendoring the SDSL v3 (xxsds fork).
# It does not attempt to detect or use a system-installed sdsl.

include(FetchContent)
# Force FetchContent to use the build directory for downloaded/populated subprojects.
# When building in an out-of-source or copied source tree (for example conda-build),
# existing _skbuild/_deps present in the original source can cause FetchContent to
# reuse CMakeCache files from the wrong directory and fail with a CMakeCache.txt
# mismatch. Pointing FETCHCONTENT_BASE_DIR to the current CMake binary dir ensures
# subbuilds are created inside the active build tree and avoids the mismatch.
set(FETCHCONTENT_BASE_DIR "${CMAKE_BINARY_DIR}/_deps")

# Allow overriding the sdsl tag for advanced users, but default to a stable release.
set(SDSL_LITE_GIT_TAG "v3.0.3" CACHE STRING "Tag/branch to fetch for sdsl-lite (xxsds fork)")
message(STATUS "Vendoring sdsl-lite (xxsds) tag ${SDSL_LITE_GIT_TAG}")
FetchContent_Declare(sdsl-lite
  GIT_REPOSITORY https://github.com/xxsds/sdsl-lite.git
  GIT_TAG        ${SDSL_LITE_GIT_TAG}
)
## Populate sdsl-lite source only (do not run its CMake). The upstream sdsl
## CMake enforces older policies/versions and may fail when configured inside
## containerized build environments (conda-build, scikit-build). We only need
## the headers for this project, so avoid invoking sdsl-lite's CMake.
FetchContent_GetProperties(sdsl-lite)
if (NOT sdsl-lite_POPULATED)
  message(STATUS "Populating sdsl-lite source (no configure)")
  FetchContent_Populate(sdsl-lite)
endif()

# Expose sdsl-lite's headers via sdsl-lite_SOURCE_DIR (FetchContent_Populate
# sets this variable when the content is downloaded).
if (NOT DEFINED sdsl-lite_SOURCE_DIR)
  message(FATAL_ERROR "sdsl-lite source directory not available after FetchContent_Populate")
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
find_package(pybind11 CONFIG QUIET)
if (NOT pybind11_FOUND)
  message(STATUS "Vendoring pybind11 via FetchContent")
  FetchContent_Declare(pybind11
    GIT_REPOSITORY https://github.com/pybind/pybind11.git
    GIT_TAG        v2.10.4
  )
  FetchContent_MakeAvailable(pybind11)
endif()

add_library(noLZSS_core STATIC src/cpp/factorizer.cpp src/cpp/fasta_processor.cpp)
target_include_directories(noLZSS_core PUBLIC src/cpp)
target_link_libraries(noLZSS_core PUBLIC ${SDSL_TARGET})

pybind11_add_module(_noLZSS src/cpp/bindings.cpp)
target_include_directories(_noLZSS PRIVATE src/cpp)
add_dependencies(_noLZSS noLZSS_core)
target_link_libraries(_noLZSS PRIVATE noLZSS_core)

# Add version definition for C++ code
target_compile_definitions(_noLZSS PRIVATE NOLZSS_VERSION="${PROJECT_VERSION}")

# Configure Doxyfile from template to ensure version consistency
configure_file(Doxyfile.in Doxyfile @ONLY)

install(TARGETS _noLZSS
  LIBRARY DESTINATION noLZSS
  RUNTIME DESTINATION noLZSS
  ARCHIVE DESTINATION noLZSS
)

# Install Python source files
install(FILES 
  src/noLZSS/__init__.py
  src/noLZSS/core.py
  src/noLZSS/utils.py
  DESTINATION noLZSS
)

# Install genomics subpackage
install(FILES 
  src/noLZSS/genomics/__init__.py
  src/noLZSS/genomics/fasta.py
  src/noLZSS/genomics/sequences.py
  DESTINATION noLZSS/genomics
)
