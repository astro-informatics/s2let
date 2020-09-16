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

add_library(s2let STATIC ${C_SRC})
target_link_libraries(s2let PUBLIC ${SO3_TARGET})
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

  install(
    FILES ${PROJECT_SOURCE_DIR}/include/s2let.h
          ${PROJECT_SOURCE_DIR}/include/s2let_error.h
          ${PROJECT_SOURCE_DIR}/include/s2let_hpx.h
          ${PROJECT_SOURCE_DIR}/include/s2let_lm.h
          ${PROJECT_SOURCE_DIR}/include/s2let_mw.h
          ${PROJECT_SOURCE_DIR}/include/s2let_tiling.h
          ${PROJECT_SOURCE_DIR}/include/s2let_transform_axisym_mw.h
          ${PROJECT_SOURCE_DIR}/include/s2let_alloc.h
          ${PROJECT_SOURCE_DIR}/include/s2let_fits.h
          ${PROJECT_SOURCE_DIR}/include/s2let_idl_hpx.h
          ${PROJECT_SOURCE_DIR}/include/s2let_math.h
          ${PROJECT_SOURCE_DIR}/include/s2let_so3.h
          ${PROJECT_SOURCE_DIR}/include/s2let_transform_axisym_hpx.h
          ${PROJECT_SOURCE_DIR}/include/s2let_types.h
          ${PROJECT_SOURCE_DIR}/include/s2let_analysis.h
          ${PROJECT_SOURCE_DIR}/include/s2let_helper.h
          ${PROJECT_SOURCE_DIR}/include/s2let_idl_mw.h
          ${PROJECT_SOURCE_DIR}/include/s2let_mex.h
          ${PROJECT_SOURCE_DIR}/include/s2let_synthesis.h
          ${PROJECT_SOURCE_DIR}/include/s2let_transform_axisym_lm.h
          ${PROJECT_SOURCE_DIR}/include/s2let_analysis_adjoint.h
          ${PROJECT_SOURCE_DIR}/include/s2let_synthesis_adjoint.h
    DESTINATION include)
endif()
