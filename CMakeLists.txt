cmake_minimum_required(VERSION 3.24)

if(NOT SKBUILD_PROJECT_NAME)
  set(project_name s2let)
else()
  set(project_name ${SKBUILD_PROJECT_NAME})
endif()

project(
  ${project_name}
  VERSION "2.2.7"
  DESCRIPTION "Fast wavelets on the sphere"
  HOMEPAGE_URL "http://astro-informatics.github.io/s2let/"
  LANGUAGES C)

option(cfitsio "Links to cfitsio" OFF)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "Debug")
endif()

include(CTest)
if(EXISTS "${PROJECT_BINARY_DIR}/FindFFTW3.cmake")
  list(APPEND CMAKE_MODULE_PATH "${PROJECT_BINARY_DIR}")
else()
  list(APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake")
endif()

find_package(astro-informatics-so3 REQUIRED)
find_library(MATH_LIBRARY m)
if(cfitsio)
  find_package(cfitsio REQUIRED)
endif()

add_subdirectory(src/c)

include(CTest)
if(BUILD_TESTING)
  include("${PROJECT_SOURCE_DIR}/cmake/fetch_cmocka.cmake")
  add_subdirectory(test/c)
endif()

if(SKBUILD)
  find_package(
    Python
    COMPONENTS Interpreter Development.Module NumPy
    REQUIRED)
  add_subdirectory(src/pys2let)
else()
  include("${PROJECT_SOURCE_DIR}/cmake/exporting.cmake")
endif()

# only run documentation if this is not a sub-project
if(CMAKE_SOURCE_DIR STREQUAL CMAKE_CURRENT_SOURCE_DIR AND NOT SKBUILD)
  find_package(Doxygen)
  if(DOXYGEN_FOUND)
    set(DOXYGEN_PROJECT_LOGO ${PROJECT_SOURCE_DIR}/docs/images/wav.png)
    set(DOXYGEN_OPTIMIZE_OUTPUT_FOR_C "YES")
    set(DOXYGEN_EXTRACT_ALL "YES")
    doxygen_add_docs(docs src/main/c src/test/c include)
  endif()
endif()
