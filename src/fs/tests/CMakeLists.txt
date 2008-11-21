add_executable (ex1 ex1.c)
target_link_libraries(ex1 ${Dohp_LIBRARIES})
add_executable (ex2 ex2.c)
target_link_libraries(ex2 ${Dohp_LIBRARIES})
add_executable (ex3 ex3.c)
target_link_libraries(ex3 ${Dohp_LIBRARIES})

add_test (fs-mf-projection ex1 -snes_mf -ksp_type minres -const_BDeg 6 -snes_monitor -ksp_converged_reason)
