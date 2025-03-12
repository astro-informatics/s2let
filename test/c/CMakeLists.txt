foreach(name tiling adjoint back_and_forth axisym_vs_directional)
  add_executable(test_${name} test_${name}.c)
  target_link_libraries(test_${name} PUBLIC s2let cmocka)
  add_test(NAME test_${name} COMMAND test_${name})
  set_target_properties(
    test_${name} PROPERTIES C_STANDARD 11 RUNTIME_OUTPUT_DIRECTORY
                                          ${PROJECT_BINARY_DIR}/bin)
endforeach()
