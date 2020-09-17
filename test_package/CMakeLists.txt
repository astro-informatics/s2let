cmake_minimum_required(VERSION 3.12)
project(PackageTest C)

include(${CMAKE_BINARY_DIR}/conanbuildinfo.cmake)
conan_basic_setup(TARGETS)

add_executable(example example.c)
target_include_directories(example PUBLIC ${PROJECT_SOURCE_DIR})
target_link_libraries(example PUBLIC CONAN_PKG::s2let)

add_test(
  NAME example
  WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/bin
  COMMAND example)
