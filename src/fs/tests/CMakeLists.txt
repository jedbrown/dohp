macro (dohp_link_executable name source)
  add_executable (${name} ${source})
  target_link_libraries (${name} dohp)
endmacro ()

dohp_link_executable (fs-ex1 ex1.c)
dohp_link_executable (ellip ellip.c)
# add_executable (ex2 ex2.c)
# target_link_libraries(ex2 ${Dohp_LIBRARIES})
# add_executable (ex3 ex3.c)
#target_link_libraries(ex3 ${Dohp_LIBRARIES})

add_test (fs-mf-projection fs-ex1 -snes_mf -ksp_type minres -const_BDeg 6 -snes_monitor -ksp_converged_reason)
add_test (fs-mf-op-proj-0  fs-ex1 -snes_mf_operator -ksp_type minres -pc_type jacobi -const_BDeg 10 -snes_monitor -ksp_converged_reason
  -ksp_rtol 1e-10 -q1scale 0.5 -frequency 5,4,3 -require_ptwise 2e-6,3e-6,2.3e-5 -require_grad 7e-7,3e-6,2e-4)
add_test (fs-mf-op-proj-1  fs-ex1 -const_BDeg 7 -snes_mf_operator -ksp_converged_reason -pc_type jacobi -ksp_type minres
  -ksp_rtol 1e-10 -snes_monitor -frequency 2,2,2 -require_ptwise 2e-6,2e-6,9e-6 -require_grad 3e-6,7e-6,1.1e-4)
