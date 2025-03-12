add_custom_command(
  OUTPUT pys2let.c
  MAIN_DEPENDENCY pys2let.pyx
  VERBATIM
  COMMAND
    Python::Interpreter -m cython "${CMAKE_CURRENT_SOURCE_DIR}/pys2let.pyx"
    --output-file "${CMAKE_CURRENT_BINARY_DIR}/pys2let.c")

python_add_library(pys2let MODULE "${CMAKE_CURRENT_BINARY_DIR}/pys2let.c"
                   WITH_SOABI)
if(WIN32)
  target_compile_definitions(pys2let PRIVATE MS_WIN64)
endif()

target_link_libraries(pys2let PUBLIC s2let)
target_include_directories(pys2let PUBLIC ${Python3_NumPy_INCLUDE_DIR})
install(TARGETS pys2let LIBRARY DESTINATION ${SKBUILD_PROJECT_NAME})
