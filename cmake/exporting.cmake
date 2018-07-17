# Exports s2let so other packages can access it
export(TARGETS s2let FILE "${PROJECT_BINARY_DIR}/S2letTargets.cmake")

# Avoids creating an entry in the cmake registry.
if(NOT NOEXPORT)
  export(PACKAGE s2let)
endif()

# First in binary dir
set(ALL_INCLUDE_DIRS "${PROJECT_SOURCE_DIR}")
configure_File(cmake/S2letConfig.in.cmake
  "${PROJECT_BINARY_DIR}/S2letConfig.cmake" @ONLY
)
configure_File(cmake/S2letConfigVersion.in.cmake
  "${PROJECT_BINARY_DIR}/S2letConfigVersion.cmake" @ONLY
)

# Then for installation tree
file(RELATIVE_PATH REL_INCLUDE_DIR
    "${CMAKE_INSTALL_PREFIX}/share/cmake/s2let"
    "${CMAKE_INSTALL_PREFIX}/include/s2let"
)
set(ALL_INCLUDE_DIRS "\${s2let_CMAKE_DIR}/${REL_INCLUDE_DIR}")
configure_file(cmake/S2letConfig.in.cmake
  "${PROJECT_BINARY_DIR}/CMakeFiles/S2letConfig.cmake" @ONLY
)

# Finally install all files
install(FILES
  "${PROJECT_BINARY_DIR}/CMakeFiles/S2letConfig.cmake"
  "${PROJECT_BINARY_DIR}/S2letConfigVersion.cmake"
  DESTINATION share/cmake/s2let
    COMPONENT dev
)

install(EXPORT S2letTargets DESTINATION share/cmake/s2let COMPONENT dev)
