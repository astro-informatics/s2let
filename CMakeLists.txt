cmake_minimum_required(VERSION 3.12)

project(
  s2let
  VERSION "2.2.0"
  DESCRIPTION "Fast wavelets on the sphere"
  HOMEPAGE_URL "http://astro-informatics.github.io/s2let/"
  LANGUAGES C)

option(tests "Enable testing" ON)
option(python "Creates python package only" OFF)
option(cfitsio "Links to cfitsio" OFF)
option(fPIC "Compile with fPIC" ON)
option(conan_deps "Download dependencies using conan" OFF)

if(WIN32)
  set(fpic OFF)
endif()
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "Debug")
endif()

list(APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)
include("sanitizers")
if(conan_deps OR CONAN_EDITABLE_MODE)
  include("conan_dependencies")
  if(NOT CONAN_EDITABLE_MODE)
    set(SO3_TARGET CONAN_PKG::so3)
  endif()
  if(cfitsio)
    set(CFITSIO_TARGET CONAN_PKG::cfitsio)
  endif()
endif()
if(NOT conan_deps)
  find_package(So3 REQUIRED)
  set(SO3_TARGET so3)
  if(cfitsio AND NOT CONAN_EDITABLE_MODE)
    find_package(CFITSIO REQUIRED)
    set(CFITSIO_TARGET cfitsio::cfitsio)
  endif()
endif()

add_subdirectory(src/main/c)

if(tests)
  enable_testing()
  include("fetch_cmocka")
  add_subdirectory(src/test/c)
endif()

if(NOT python)
  include("exporting")
else()
  include("python-setup")
  find_python()
  setup_skbuild()
  find_package(PythonExtensions REQUIRED)
  find_package(Cython REQUIRED)
  find_package(NumPy REQUIRED)
  add_subdirectory(src/main/pys2let)
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
