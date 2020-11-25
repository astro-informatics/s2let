if(NOT EXISTS "${CMAKE_BINARY_DIR}/conan.cmake")
  message(
    STATUS
      "Downloading conan.cmake from https://github.com/conan-io/cmake-conan")
  file(DOWNLOAD "https://github.com/conan-io/cmake-conan/raw/v0.15/conan.cmake"
       "${CMAKE_BINARY_DIR}/conan.cmake" TLS_VERIFY ON)
endif()
include(${CMAKE_BINARY_DIR}/conan.cmake)

if(fPIC AND NOT CONAN_OPTIONS)
  set(fpic_value "True")
elseif(NOT CONAN_OPTIONS)
  set(fpic_value "False")
endif()
if(NOT CONAN_EDITABLE_MODE)
  if(NOT CONAN_DEPS)
    set(CONAN_DEPS "so3/1.3.1@astro-informatics/stable")
  endif()
  list(APPEND CONAN_OPTIONS "so3:fPIC=${fpic_value}")
endif()
if(cfitsio)
  list(APPEND CONAN_DEPS "cfitsio/3.490")
  list(APPEND CONAN_OPTIONS "cfitsio:shared=False" "cfitsio:fPIC=${fpic_value}")
endif()
if(NOT CONAN_BUILD)
  set(CONAN_BUILD "missing")
endif()

conan_check(REQUIRED)
conan_add_remote(
  NAME astro-informatics URL
  https://api.bintray.com/conan/astro-informatics/astro-informatics)
conan_cmake_run(
  REQUIRES
  ${CONAN_DEPS}
  BASIC_SETUP
  OPTIONS
  "${CONAN_OPTIONS}"
  KEEP_RPATHS
  CMAKE_TARGETS
  NO_OUTPUT_DIRS
  BUILD
  ${CONAN_BUILD})
