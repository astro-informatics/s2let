add_cython_target(cython_sources pys2let.pyx C PY3)

add_library(pys2let MODULE ${cython_sources})
python_extension_module(pys2let)
if(WIN32)
  target_compile_definitions(pys2let PRIVATE MS_WIN64)
endif()

target_link_libraries(pys2let s2let)
target_include_directories(pys2let PUBLIC "${PROJECT_SOURCE_DIR}/include"
                                          ${NumPy_INCLUDE_DIRS})
install(TARGETS pys2let LIBRARY DESTINATION src/main/pys2let/)
