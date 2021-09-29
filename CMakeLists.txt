cmake_minimum_required(VERSION 3.12...4.0)

project(
  s2let
  VERSION "2.2.4"
  DESCRIPTION "Fast wavelets on the sphere"
  HOMEPAGE_URL "http://astro-informatics.github.io/s2let/"
  LANGUAGES C)

option(cfitsio "Links to cfitsio" OFF)
option(conan_deps "Download dependencies using conan" OFF)

if(conan_deps AND NOT cfitsio)
  message(
    WARNING "cfitsio is always included when managing dependencies via conan")
endif()

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "Debug")
endif()

include(CTest)
if(conan_deps OR CONAN_EDITABLE_MODE)
  include("${PROJECT_SOURCE_DIR}/cmake/conan_dependencies.cmake")
endif()
if(EXISTS "${PROJECT_BINARY_DIR}/conan_paths.cmake")
  include("${PROJECT_BINARY_DIR}/conan_paths.cmake")
elseif(EXISTS "${PROJECT_BINARY_DIR}/FindFFTW3.cmake")
  list(APPEND CMAKE_MODULE_PATH "${PROJECT_BINARY_DIR}")
else()
  list(APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake")
endif()

find_package(astro-informatics-so3 REQUIRED)
find_library(MATH_LIBRARY m)
if(cfitsio)
  find_package(cfitsio REQUIRED)
endif()

add_subdirectory(src/main/c)

include(CTest)
if(BUILD_TESTING)
  include("${PROJECT_SOURCE_DIR}/cmake/fetch_cmocka.cmake")
  add_subdirectory(src/test/c)
endif()

if(SKBUILD)
  find_package(PythonExtensions REQUIRED)
  find_package(Cython REQUIRED)
  find_package(NumPy REQUIRED)
  add_subdirectory(src/main/pys2let)
elseif(NOT CONAN_EXPORTED)
  include("${PROJECT_SOURCE_DIR}/cmake/exporting.cmake")
endif()

# only run documentation if this is not a sub-project
if(CMAKE_SOURCE_DIR STREQUAL CMAKE_CURRENT_SOURCE_DIR)
  find_package(Doxygen)
  if(DOXYGEN_FOUND)
    set(DOXYGEN_PROJECT_LOGO ${PROJECT_SOURCE_DIR}/docs/images/wav.png)
    set(DOXYGEN_OPTIMIZE_OUTPUT_FOR_C "YES")
    set(DOXYGEN_EXTRACT_ALL "YES")
    doxygen_add_docs(docs src/main/c src/test/c include)
  endif()
endif()
