add_executable(s2let_test s2let_test.c)
target_link_libraries(s2let_test PUBLIC s2let)
add_test(NAME s2let_test COMMAND s2let_test)

foreach(name tiling adjoint back_and_forth axisym_vs_directional)
  add_executable(test_${name} test_${name}.c)
  target_link_libraries(test_${name} PUBLIC s2let cmocka)
  add_test(NAME test_${name} COMMAND test_${name})
endforeach()
