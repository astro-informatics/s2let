find_path(So3_INCLUDE_DIR NAMES ssht.h PATH_SUFFIXES ssht)
find_library(So3_LIBRARY NAMES ssht)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(
  So3  DEFAULT_MSG So3_INCLUDE_DIR So3_LIBRARY)
MESSAGE(STATUS "Using ssht include ${So3_INCLUDE_DIR}")
MESSAGE(STATUS "Using ssht lib ${So3_LIBRARY}")
