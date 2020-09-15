cmake_minimum_required(VERSION 3.12)

project(s2let
  VERSION "2.2.0"
  DESCRIPTION "Fast wavelets on the sphere"
  HOMEPAGE_URL "http://astro-informatics.github.io/s2let/"
  LANGUAGES C)

option(tests "Enable testing" ON)
option(python "Creates python package" OFF)
option(fPIC "Compile with fPIC" ON)
option(conan_deps "Download dependencies using conan" OFF)

list(APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)
if(conan_deps)
  include("conan_dependencies")
  set(SO3_TARGET CONAN_PKG::so3)
else()
  find_package(So3 REQUIRED)
  set(SO3_TARGET So3)
endif()

add_subdirectory(src/main/c)

if(tests)
    enable_testing()
    add_subdirectory(src/test/c)
endif()

include("exporting")
