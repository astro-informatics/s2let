if(NOT EXISTS "${CMAKE_BINARY_DIR}/conan.cmake")
  message(
    STATUS
      "Downloading conan.cmake from https://github.com/conan-io/cmake-conan")
  file(DOWNLOAD
       "https://github.com/conan-io/cmake-conan/raw/v0.16.1/conan.cmake"
       "${CMAKE_BINARY_DIR}/conan.cmake" TLS_VERIFY ON)
endif()
include(${CMAKE_BINARY_DIR}/conan.cmake)

conan_cmake_configure(CONANFILE ${PROJECT_SOURCE_DIR}/conanfile.txt GENERATORS
                      cmake_find_package)
conan_cmake_autodetect(settings)
conan_cmake_install(
  PATH_OR_REFERENCE
  ${PROJECT_SOURCE_DIR}/conanfile.txt
  GENERATOR
  cmake_find_package
  cmake_paths
  BUILD
  missing
  REMOTE
  conancenter
  SETTINGS
  ${settings})
