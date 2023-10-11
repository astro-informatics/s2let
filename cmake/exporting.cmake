# Exports s2let so other packages can access it
export(
  TARGETS s2let
  FILE "${PROJECT_BINARY_DIR}/s2letTargets.cmake"
  NAMESPACE s2let::)

# Avoids creating an entry in the cmake registry.
if(NOT NOEXPORT)
  export(PACKAGE s2let)
endif()

# First in binary dir
set(INCLUDE_INSTALL_DIR include/)
include(CMakePackageConfigHelpers)
configure_package_config_file(
  cmake/S2letConfig.in.cmake "${PROJECT_BINARY_DIR}/s2letConfig.cmake"
  INSTALL_DESTINATION lib/cmake/s2let
  PATH_VARS INCLUDE_INSTALL_DIR)
write_basic_package_version_file(
  s2letConfigVersion.cmake
  VERSION ${PROJECT_VERSION}
  COMPATIBILITY SameMajorVersion)

if(NOT CONAN_EXPORTED)
  install(FILES "${PROJECT_BINARY_DIR}/s2letConfig.cmake"
                "${PROJECT_BINARY_DIR}/s2letConfigVersion.cmake"
          DESTINATION lib/cmake/s2let)
endif()

install(
  EXPORT s2letTargets
  DESTINATION lib/cmake/s2let
  NAMESPACE s2let::)
