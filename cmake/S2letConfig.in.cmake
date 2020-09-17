get_filename_component(S2let_CMAKE_DIR "${CMAKE_CURRENT_LIST_FILE}" PATH)
message(STATUS "Linking to S2let package in ${S2let_CMAKE_DIR}")
if(NOT TARGET S2let AND EXISTS "${S2let_CMAKE_DIR}/S2letTargets.cmake")
  include("${S2let_CMAKE_DIR}/S2letTargets.cmake")
endif()

set(S2let_INCLUDE_DIRS "@ALL_INCLUDE_DIRS@")
set(S2let_LIBRARIES s2let)
