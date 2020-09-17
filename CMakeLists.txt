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

list(APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)
if(conan_deps)
  include("conan_dependencies")
  set(SO3_TARGET CONAN_PKG::so3)
  if(cfitsio)
    set(CFITSIO_TARGET CONAN_PKG::cfitsio)
  endif()
else()
  find_package(So3 REQUIRED)
  if(cfitsio)
    find_package(CFITSIO REQUIRED)
    set(CFITSIO_TARGET cfitsio::cfitsio)
  endif()
  set(SO3_TARGET So3)
endif()

add_subdirectory(src/main/c)

if(tests)
  enable_testing()
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
