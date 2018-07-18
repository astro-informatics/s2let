find_path(So3_INCLUDE_DIR NAMES so3.h PATH_SUFFIXES so3)
find_library(So3_LIBRARY NAMES so3)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(
  So3  DEFAULT_MSG So3_INCLUDE_DIR So3_LIBRARY)
MESSAGE(STATUS "Using so3 include ${So3_INCLUDE_DIR}")
MESSAGE(STATUS "Using so3 lib ${So3_LIBRARY}")
