add_executable(s2let_test s2let_test.c)
target_link_libraries(s2let_test PUBLIC s2let)
add_test(NAME s2let_test COMMAND s2let_test)
