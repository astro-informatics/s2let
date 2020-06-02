find_path(Ssht_INCLUDE_DIR NAMES ssht.h HINTS /usr/local/opt/src_ssht/include/c)
find_library(Ssht_LIBRARY NAMES libssht.a HINTS /usr/local/opt/src_ssht/build/src/c)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(
  Ssht  DEFAULT_MSG Ssht_INCLUDE_DIR Ssht_LIBRARY)
MESSAGE(STATUS "Using ssht include ${Ssht_INCLUDE_DIR}")
MESSAGE(STATUS "Using ssht lib ${Ssht_LIBRARY}")
