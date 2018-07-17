find_path(Ssht_INCLUDE_DIR NAMES ssht.h PATH_SUFFIXES ssht)
find_library(Ssht_LIBRARY NAMES ssht)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(
  Ssht  DEFAULT_MSG Ssht_INCLUDE_DIR Ssht_LIBRARY)
MESSAGE(STATUS "Using ssht include ${Ssht_INCLUDE_DIR}")
MESSAGE(STATUS "Using ssht lib ${Ssht_LIBRARY}")
