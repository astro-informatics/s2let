set(C_SRC
    s2let_transform_axisym_lm.c
    s2let_transform_axisym_mw.c
    s2let_math.c
    s2let_mw.c
    s2let_tiling.c
    s2let_synthesis.c
    s2let_synthesis_adjoint.c
    s2let_lm.c
    s2let_alloc.c
    s2let_helper.c
    s2let_analysis.c
    s2let_analysis_adjoint.c)
set(HEADER_PREFIX ${PROJECT_SOURCE_DIR}/include)
set(C_HEADERS
    s2let.h
    s2let_error.h
    s2let_hpx.h
    s2let_lm.h
    s2let_mw.h
    s2let_tiling.h
    s2let_transform_axisym_mw.h
    s2let_alloc.h
    s2let_idl_hpx.h
    s2let_math.h
    s2let_so3.h
    s2let_transform_axisym_hpx.h
    s2let_types.h
    s2let_analysis.h
    s2let_helper.h
    s2let_idl_mw.h
    s2let_mex.h
    s2let_synthesis.h
    s2let_transform_axisym_lm.h
    s2let_analysis_adjoint.h
    s2let_synthesis_adjoint.h
    s2let_fits.h)

if(cfitsio)
  list(APPEND C_SRC s2let_fits.c)
endif()
list(TRANSFORM C_HEADERS PREPEND "${PROJECT_SOURCE_DIR}/include/")

add_library(s2let STATIC ${C_SRC})
target_link_libraries(s2let PUBLIC ${SO3_TARGET} ${CFITSIO_TARGET})
set_target_properties(s2let PROPERTIES C_STANDARD 99)
target_include_directories(
  s2let PUBLIC $<INSTALL_INTERFACE:include>
               $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>)
if(fPIC)
  set_target_properties(s2let PROPERTIES POSITION_INDEPENDENT_CODE TRUE)
endif()

if(NOT python)
  add_executable(s2let_about s2let_about.c)
  target_compile_definitions(
    s2let_about PRIVATE S2LET_VERSION="${PROJECT_VERSION}" S2LET_BUILD="CMake")
  set_target_properties(s2let_about PROPERTIES RUNTIME_OUTPUT_DIRECTORY
                                               ${PROJECT_BINARY_DIR})

  install(
    TARGETS s2let
    EXPORT S2letTargets
    ARCHIVE DESTINATION lib
    PUBLIC_HEADER)

  install(FILES ${C_HEADERS} DESTINATION include/s2let)
endif()
