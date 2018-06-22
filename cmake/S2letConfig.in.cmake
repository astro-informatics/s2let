get_filename_component(So3_CMAKE_DIR "${CMAKE_CURRENT_LIST_FILE}" PATH)
message(STATUS "Linking to so3 package in ${So3_CMAKE_DIR}")
if(NOT TARGET so3 AND EXISTS "${So3_CMAKE_DIR}/So3Targets.cmake")
  include("${So3_CMAKE_DIR}/So3Targets.cmake")
endif()

set(So3_INCLUDE_DIRS "@ALL_INCLUDE_DIRS@")
set(So3_LIBRARIES so3)
