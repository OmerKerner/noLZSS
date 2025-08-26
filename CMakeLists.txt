cmake_minimum_required(VERSION 3.20)
project(noLZSS LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

option(NOLZSS_FORCE_VENDOR_SDSL "Always fetch (vendor) sdsl-lite even if system copy exists" OFF)
option(NOLZSS_FORCE_SYSTEM_SDSL "Require system sdsl-lite (do not fetch)" OFF)
set(SDSL_ROOT "" CACHE PATH "Prefix of a pre-installed sdsl-lite (optional)")

set(_sdsl_found FALSE)
if (NOT NOLZSS_FORCE_VENDOR_SDSL)
  if (SDSL_ROOT)
    list(APPEND CMAKE_PREFIX_PATH "${SDSL_ROOT}")
  endif()
  find_path(SDSL_INCLUDE_DIR NAMES sdsl/suffix_trees.hpp)
  find_library(SDSL_LIBRARY NAMES sdsl)
  if (SDSL_INCLUDE_DIR AND SDSL_LIBRARY)
    set(_sdsl_found TRUE)
    message(STATUS "Using system sdsl-lite: inc=${SDSL_INCLUDE_DIR} lib=${SDSL_LIBRARY}")
  endif()
endif()

if (NOT _sdsl_found)
  if (NOLZSS_FORCE_SYSTEM_SDSL)
    message(FATAL_ERROR "NOLZSS_FORCE_SYSTEM_SDSL=ON but system sdsl-lite not found. Provide SDSL_ROOT or install the package.")
  endif()
  include(FetchContent)
  # By default fetch the main branch of sdsl-lite. For reproducible builds you
  # can set SDSL_LITE_GIT_TAG to a specific commit hash or tag (for example
  # via -DSDSL_LITE_GIT_TAG=<commit> when configuring). The previous hard-
  # coded commit used here was not available and caused fetch failures.
  # Pin to a recent stable tag to ensure reproducible and reliable fetches.
  # Update this value when you intentionally upgrade sdsl-lite.
  set(SDSL_LITE_GIT_TAG "v2.1.1")
  message(STATUS "Fetching vendored sdsl-lite tag ${SDSL_LITE_GIT_TAG}")
  FetchContent_Declare(sdsl-lite
    GIT_REPOSITORY https://github.com/simongog/sdsl-lite.git
    GIT_TAG        ${SDSL_LITE_GIT_TAG}
  )
  set(BUILD_EXAMPLES OFF CACHE BOOL "" FORCE)
  set(BUILD_TESTING OFF CACHE BOOL "" FORCE)
  FetchContent_MakeAvailable(sdsl-lite)
  if (NOT TARGET sdsl)
    find_library(SDSL_LIBRARY_FALLBACK NAMES sdsl PATHS ${sdsl-lite_SOURCE_DIR}/lib ${sdsl-lite_BINARY_DIR}/lib)
    add_library(sdsl STATIC IMPORTED)
    set_target_properties(sdsl PROPERTIES
      IMPORTED_LOCATION "${SDSL_LIBRARY_FALLBACK}"
      INTERFACE_INCLUDE_DIRECTORIES "${sdsl-lite_SOURCE_DIR}/include")
  endif()
endif()

find_package(pybind11 CONFIG REQUIRED)

add_library(noLZSS_core STATIC src/noLZSS/factorizer.cpp)

if (TARGET sdsl)
  target_link_libraries(noLZSS_core PUBLIC sdsl)
elseif(_sdsl_found)
  target_include_directories(noLZSS_core PUBLIC ${SDSL_INCLUDE_DIR})
  target_link_libraries(noLZSS_core PUBLIC ${SDSL_LIBRARY})
else()
  message(FATAL_ERROR "sdsl target not found and system library missing")
endif()

pybind11_add_module(_noLZSS src/noLZSS/bindings.cpp)
add_dependencies(_noLZSS noLZSS_core)
target_link_libraries(_noLZSS PRIVATE noLZSS_core)

install(TARGETS _noLZSS DESTINATION .)
