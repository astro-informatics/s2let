cmake_minimum_required(VERSION 3.7)
option(openmp         "Enable OpenMP"                                   on)
option(fPIC           "Enable fPIC (only for gcc)"                      on)
option(tests          "Enable testing"                                  on)

project(s2let C)


set(S2let_VERSION "1.1b1")
set(S2let_BUILD "CMake")
set(CMAKE_C_STANDARD 99)

list(APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)

find_package(FFTW3 REQUIRED COMPONENT SERIAL DOUBLE)
find_package(OpenMP QUIET)
find_package(Ssht REQUIRED)
find_package(So3 REQUIRED)

list(APPEND C_SRC src/main/c/s2let_transform_axisym_lm.c
  src/main/c/s2let_transform_axisym_mw.c src/main/c/s2let_math.c
  src/main/c/s2let_mw.c src/main/c/s2let_tiling.c
  src/main/c/s2let_synthesis.c
  src/main/c/s2let_synthesis_adjoint.c
  src/main/c/s2let_lm.c
  src/main/c/s2let_alloc.c
  src/main/c/s2let_helper.c
  src/main/c/s2let_analysis.c
  src/main/c/s2let_analysis_adjoint.c
  )

add_library(s2let STATIC ${C_SRC})
target_include_directories(s2let PUBLIC ${FFTW3_INCLUDE_DIRS})
target_include_directories(s2let PUBLIC ${Ssht_INCLUDE_DIR})
target_include_directories(s2let PUBLIC ${So3_INCLUDE_DIR})
target_include_directories(s2let PRIVATE ${PROJECT_SOURCE_DIR}/include/)
if(OpenMP_C_FOUND AND openmp)
  target_link_libraries(s2let OpenMP::OpenMP_C)
endif()
target_link_libraries(s2let ${Ssht_LIBRARY})
target_link_libraries(s2let ${So3_LIBRARY})
target_link_libraries(s2let ${FFTW3_DOUBLE_SERIAL_LIBRARY})

if(fPIC)
  set_target_properties(s2let PROPERTIES COMPILE_FLAGS "${FFTW3_DEFINITIONS} -fPIC")
else()
  set_target_properties(s2let PROPERTIES COMPILE_FLAGS "${FFTW3_DEFINITIONS}")
endif()

install(TARGETS s2let EXPORT S2letTargets ARCHIVE DESTINATION lib PUBLIC_HEADER)
install(FILES
  include/s2let.h                      include/s2let_error.h                include/s2let_hpx.h                  include/s2let_lm.h                   include/s2let_mw.h                   include/s2let_tiling.h               include/s2let_transform_axisym_mw.h
  include/s2let_alloc.h                include/s2let_fits.h                 include/s2let_idl_hpx.h              include/s2let_math.h                 include/s2let_so3.h                  include/s2let_transform_axisym_hpx.h include/s2let_types.h
  include/s2let_analysis.h             include/s2let_helper.h               include/s2let_idl_mw.h               include/s2let_mex.h                  include/s2let_synthesis.h            include/s2let_transform_axisym_lm.h
  include/s2let_analysis_adjoint.h     include/s2let_synthesis_adjoint.h
  DESTINATION include)
#building about and test executables
add_executable(s2let_about src/main/c/s2let_about.c)
target_compile_definitions(s2let_about PRIVATE S2LET_VERSION="${S2let_VERSION}" S2LET_BUILD="${S2let_BUILD}")
set_target_properties(s2let_about PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR})

include("exporting")

if(tests)
  add_executable(s2let_test src/test/c/s2let_test.c)
  target_include_directories(s2let_test PRIVATE ${PROJECT_SOURCE_DIR}/include/)
  target_link_libraries(s2let_test s2let)
endif()
